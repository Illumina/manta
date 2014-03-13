// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#include "GenerateSVCandidates.hh"
#include "EdgeRetrieverBin.hh"
#include "EdgeRetrieverLocus.hh"
#include "EdgeRuntimeTracker.hh"
#include "GSCOptions.hh"
#include "SVCandidateProcessor.hh"
#include "SVFinder.hh"

#include "blt_util/log.hh"
#include "blt_util/bam_header_util.hh"

#include "common/Exceptions.hh"
#include "manta/ReadGroupStatsSet.hh"
#include "manta/SVCandidateUtil.hh"
#include "manta/SVMultiJunctionCandidate.hh"
#include "truth/TruthTracker.hh"

#include "boost/foreach.hpp"

#include <iostream>
#include <memory>

//#define DEBUG_GSV


/// provide additional edge details, intended for attachment to an in-flight exception:
static
void
dumpEdgeInfo(
    const EdgeInfo& edge,
    const SVLocusSet& set,
    std::ostream& os)
{
    os << edge;
    os << "\tnode1:" << set.getLocus(edge.locusIndex).getNode(edge.nodeIndex1);
    os << "\tnode2:" << set.getLocus(edge.locusIndex).getNode(edge.nodeIndex2);
}



/// we can either traverse all edges in a single locus (disjoint subgraph) of the graph
/// OR
/// traverse all edges in one "bin" -- that is, one out of binCount subsets of the total
/// graph edges. Each bin is designed to be of roughly equal size in terms of total
/// anticipated workload, so that we have good parallel processing performance.
///
static
EdgeRetriever*
edgeRFactory(
    const SVLocusSet& set,
    const EdgeOptions& opt)
{
    if (opt.isLocusIndex)
    {
        return (new EdgeRetrieverLocus(set, opt.graphNodeMaxEdgeCount, opt.locusOpt));
    }
    else
    {
        return (new EdgeRetrieverBin(set, opt.graphNodeMaxEdgeCount, opt.binCount, opt.binIndex));
    }
}



#if 0
/// edge indices+graph evidence counts and regions:
///
/// this is designed to be useful even when the locus graph is not present
struct EhancedEdgeInfo
{

};

/// reduce the full (very-large) graph down to just the information we need during SVCandidate generation:
struct ReducedGraphInfo
{
    ReducedGraphInfo(const GSCOptions& opt)

    bam_header_info header;

    std::vector<EnhancedEdgeInfo> edges;
};
#endif



/// should we filter out the SVCandidate?
///
/// Note this logic belongs in SVFinder and should make its way there once stable:
static
bool
isFilterCandidate(
    const SVCandidate& sv)
{
    // don't consider candidates created from only semi-mapped read pairs:
    if (sv.bp1.isLocalOnly() && sv.bp2.isLocalOnly()) return true;

    // candidates must have a minimum amount of evidence:
    if (isSpanningSV(sv))
    {
        static const unsigned minCandidateSpanningCount(3);
        if (sv.bp1.getSpanningCount() < minCandidateSpanningCount) return true;
    }
    else if (isComplexSV(sv))
    {
        static const unsigned minCandidateComplexCount(2);
        if (sv.bp1.lowresEvidence.getTotal() < minCandidateComplexCount) return true;
    }
    else
    {
        assert(false && "Unknown SV candidate type");
    }

    return false;
}



static
unsigned
getIntervalDist(
    const GenomeInterval& intervalA,
    const GenomeInterval& intervalB)
{
    static const unsigned far(std::numeric_limits<unsigned>::max());

    if (intervalA.tid != intervalB.tid) return far;

    return std::abs(intervalA.range.center_pos() - intervalB.range.center_pos());
}



///
static
bool
isIntervalPairGroupCandidate(
    const GenomeInterval& intervalA,
    const GenomeInterval& intervalB,
    const unsigned minFilterDist)
{
    return (getIntervalDist(intervalA,intervalB) < minFilterDist);
}



/// return  1 if dist(A1,B1) and dist(A2,B2) are both less than dist(A1,B2) and dist(A2,B1)
/// return -1 if dist(A1,B2) and dist(A2,B1) are both less than dist(A1,B1) and dist(A2,B2)
/// return 0 for all other cases
static
int
getJunctionBpAlignment(
    const SVCandidate& svA,
    const SVCandidate& svB)
{
    unsigned dist11(getIntervalDist(svA.bp1.interval,svB.bp1.interval));
    unsigned dist12(getIntervalDist(svA.bp1.interval,svB.bp2.interval));
    unsigned dist21(getIntervalDist(svA.bp2.interval,svB.bp1.interval));
    unsigned dist22(getIntervalDist(svA.bp2.interval,svB.bp2.interval));

    if (((dist11 < dist12) && (dist11 < dist21)) && ((dist22 < dist12) && (dist22 < dist21))) return  1;
    if (((dist12 < dist11) && (dist12 < dist22)) && ((dist21 < dist11) && (dist21 < dist22))) return -1;
    return 0;
}



/// are two breakend pairs candidates for a multi-junction analysis?:
///
static
bool
isBreakendPairGroupCandidate(
    const SVBreakend& bpA,
    const SVBreakend& bpB,
    const unsigned groupRange = 1000)
{
    if (! isOppositeOrientation(bpA.state, bpB.state)) return false;

    return isIntervalPairGroupCandidate(bpA.interval, bpB.interval, groupRange);
}



/// test to see if a breakend can participate in a multi-junction analysis:
///
/// right now our only criteria is to exclude small non-inversions, just because
/// such pairs can spontaneously occur at relatively high rates:
static
bool
isSVMJExcluded(
    const SVCandidate& sv)
{
    static const unsigned minInnieSVSize(5000);

    {
        using namespace SV_TYPE;
        const SV_TYPE::index_t svt(getSVType(sv));
        if ((svt != INDEL) && (svt != TANDUP)) return false;
    }

    return (getIntervalDist(sv.bp1.interval, sv.bp2.interval) < minInnieSVSize);
}



namespace MJ_INTERACTION
{
    enum index_t {
        NONE,
        SAME,
        FLIP,
        CONFLICT
    };

    struct MJState
    {
        MJState() :
            type(NONE),
            partnerId(0)
        {}

        index_t type;
        unsigned partnerId;
    };
}



static
void
findMultiJunctionCandidates(
    const std::vector<SVCandidate>& svs,
    std::vector<SVMultiJunctionCandidate>& mjSVs)
{
    mjSVs.clear();

    std::vector<SVCandidate> complexSVs;
    std::vector<SVCandidate> spanningSVs;

    BOOST_FOREACH(const SVCandidate& candidateSV, svs)
    {
        /// Filter various candidates types:
        if (isFilterCandidate(candidateSV)) continue;

        const bool isComplex(isComplexSV(candidateSV));

        if (isComplex)
        {
            complexSVs.push_back(candidateSV);
        }
        else
        {
            spanningSVs.push_back(candidateSV);
        }
    }

    /// do a brute-force intersection test to see if we can associate candidates:
    ///
    /// intersection rules : breakend region center must be within distance N
    /// intersecting breakend orientation makes it possible for these to be a single event -- ie. pointing away or towards each other
    /// full set of intersections must complete a loop, this is an intentionally conservative requirement to make sure we start into
    ///    this without getting involved in the really difficult stuff
    ///
    /// just for the starting version, the number of SVCandidates which can intersect is limited to 2
    ///

    const unsigned spanCount(spanningSVs.size());
    std::vector<MJ_INTERACTION::MJState> spanPartners(spanCount);
    {
        using namespace MJ_INTERACTION;

        for (unsigned spanIndexA(0); (spanIndexA+1)<spanCount; ++spanIndexA)
        {
            const SVCandidate& spanA(spanningSVs[spanIndexA]);
            if (isSVMJExcluded(spanA)) continue;

            for (unsigned spanIndexB(spanIndexA+1); spanIndexB<spanCount; ++spanIndexB)
            {
                const SVCandidate& spanB(spanningSVs[spanIndexB]);
                if (isSVMJExcluded(spanB)) continue;

                const bool isSameBpGroup(isBreakendPairGroupCandidate(spanA.bp1,spanB.bp1) && isBreakendPairGroupCandidate(spanA.bp2,spanB.bp2));
                const bool isFlipBpGroup(isBreakendPairGroupCandidate(spanA.bp1,spanB.bp2) && isBreakendPairGroupCandidate(spanA.bp2,spanB.bp1));

                bool isGroup(false);
                if (isSameBpGroup || isFlipBpGroup)
                {
                    /// check that this isn't a flipped association as breakpoints get near each other,
                    /// if it is treat the association as independent junctions:
                    if (isSameBpGroup)
                    {
                        isGroup = (getJunctionBpAlignment(spanA, spanB) == 1);
                    }
                    else
                    {
                        isGroup = (getJunctionBpAlignment(spanA, spanB) == -1);
                    }
                }

                if (!isGroup) continue;

                if ((spanPartners[spanIndexA].type == NONE) && (spanPartners[spanIndexB].type == NONE))
                {
                    const index_t newType( isSameBpGroup ? SAME : FLIP );
                    spanPartners[spanIndexA].type = newType;
                    spanPartners[spanIndexA].partnerId = spanIndexB;
                    spanPartners[spanIndexB].type = newType;
                    spanPartners[spanIndexB].partnerId = spanIndexA;
                }
                else
                {
                    spanPartners[spanIndexA].type = CONFLICT;
                    spanPartners[spanIndexB].type = CONFLICT;
                }
            }
        }
    }

    /// complex SVs are translated directly into single partner candidates:
    BOOST_FOREACH(const SVCandidate& candidateSV, complexSVs)
    {
        SVMultiJunctionCandidate mj;
        mj.junction.push_back(candidateSV);
        mjSVs.push_back(mj);
    }

    for (unsigned spanIndex(0); spanIndex<spanCount; ++spanIndex)
    {
        SVMultiJunctionCandidate mj;
        mj.junction.push_back(spanningSVs[spanIndex]);

        using namespace MJ_INTERACTION;
        if ((spanPartners[spanIndex].type == SAME) ||
            (spanPartners[spanIndex].type == FLIP))
        {
            const unsigned partnerId(spanPartners[spanIndex].partnerId);
            assert(partnerId < spanCount);

            // only include the connected pair once:
            if (partnerId < spanIndex) continue;

            mj.junction.push_back(spanningSVs[partnerId]);
        }
        mjSVs.push_back(mj);
    }
}



static
void
runGSC(
    const GSCOptions& opt,
    const char* progName,
    const char* progVersion)
{
#if 0
    {
        // to save memory, load the graph and process/store only the information we need from it:
    }
#endif

    SVFinder svFind(opt);
    const SVLocusSet& cset(svFind.getSet());

    TruthTracker truthTracker(opt.truthVcfFilename, cset);
    EdgeRuntimeTracker edgeTracker(opt.edgeRuntimeFilename);

    SVCandidateProcessor svProcessor(opt, progName, progVersion, cset,  truthTracker, edgeTracker);

    std::auto_ptr<EdgeRetriever> edgerPtr(edgeRFactory(cset, opt.edgeOpt));
    EdgeRetriever& edger(*edgerPtr);

    SVCandidateSetData svData;
    std::vector<SVCandidate> svs;
    std::vector<SVMultiJunctionCandidate> mjSVs;

    static const std::string logtag("runGSC");
    if (opt.isVerbose)
    {
        log_os << logtag << " " << cset.header << "\n";
    }

    const std::map<std::string, int32_t>& chromToIndex(cset.header.chrom_to_index);

    while (edger.next())
    {
        const EdgeInfo& edge(edger.getEdge());

        try
        {
            truthTracker.addEdge(edge);
            edgeTracker.start();

            if (opt.isVerbose)
            {
                log_os << logtag << " starting analysis of edge: ";
                dumpEdgeInfo(edge,cset,log_os);
            }

            // find number, type and breakend range (or better: breakend distro) of SVs on this edge:
            svFind.findCandidateSV(chromToIndex, edge, svData, svs,
                                   truthTracker);

            truthTracker.reportNumCands(svs.size(), edge);

            if (opt.isVerbose)
            {
                log_os << logtag << " Low-resolution candidate generation complete. Candidate count: " << svs.size() << "\n";
            }

            const unsigned svCount(svs.size());
            for (unsigned i(0);i<svCount;++i)
            {
                truthTracker.addCandSV();
            }

            findMultiJunctionCandidates(svs, mjSVs);

            BOOST_FOREACH(const SVMultiJunctionCandidate& mjCandidateSV, mjSVs)
            {
                svProcessor.evaluateCandidate(edge, mjCandidateSV, svData);
            }
        }
        catch (illumina::common::ExceptionData& e)
        {
            std::ostringstream oss;
            dumpEdgeInfo(edge,cset,oss);
            e << illumina::common::ExceptionMsg(oss.str());
            throw;
        }
        catch (...)
        {
            log_os << "Exception caught while processing graph component: ";
            dumpEdgeInfo(edge,cset,log_os);
            throw;
        }

        edgeTracker.stop(edge);
        if (opt.isVerbose)
        {
            log_os << logtag << " Processing this edge took " << edgeTracker.getLastEdgeTime() << " seconds.\n";
        }
    }

    truthTracker.dumpAll();
}



void
GenerateSVCandidates::
runInternal(int argc, char* argv[]) const
{
    GSCOptions opt;

    parseGSCOptions(*this,argc,argv,opt);
#ifdef DEBUG_GSV
    opt.isVerbose=true;
#endif
    runGSC(opt, name(), version());
}
