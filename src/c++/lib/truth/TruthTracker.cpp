// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Richard Shaw
///

#include <cassert>

#include "boost/foreach.hpp"

#include "blt_util/bam_header_util.hh"

#include "common/VcfFile.hh"

#include "truth/TruthTracker.hh"

/*****************************************************************************/
// ObsRecord
/*****************************************************************************/

ObsRecord::ObsRecord()
{
    _svObs = SVObservation();
    _discardReason = KEPT;
}

/*****************************************************************************/

ObsRecord::ObsRecord(const SVObservation& svObs)
{
    _svObs = svObs;
    _discardReason = KEPT;
}

/*****************************************************************************/

bool ObsRecord::discard(DiscardReason discardReason)
{
    if (_discardReason != KEPT)
    {
        return false;
    }

    _discardReason = discardReason;
    return true;
}

/*****************************************************************************/

bool ObsRecord::discarded(DiscardReason& discardReason) const
{
    discardReason = _discardReason;
    return (_discardReason != KEPT);
}

/*****************************************************************************/

bool ObsRecord::operator<(const ObsRecord& obsRecordB) const
{
    if (_svObs.candidateIndex < obsRecordB._svObs.candidateIndex)
    {
        return true;
    }

    if (_svObs.assemblyAlignIndex < obsRecordB._svObs.assemblyAlignIndex)
    {
        return true;
    }

    if (_svObs.assemblySegmentIndex < obsRecordB._svObs.assemblySegmentIndex)
    {
        return true;
    }


    return false;
}

/*****************************************************************************/

bool ObsRecord::operator==(const ObsRecord& obsRecordB) const
{
    return ((_svObs.candidateIndex == obsRecordB._svObs.candidateIndex)
            && (_svObs.assemblyAlignIndex
                == obsRecordB._svObs.assemblyAlignIndex)
            && (_svObs.assemblySegmentIndex
                == obsRecordB._svObs.assemblySegmentIndex));
}

/*****************************************************************************/

std::ostream& operator<<(std::ostream& ostrm, const ObsRecord& record)
{
    ostrm << "Candidate " << record._svObs.candidateIndex
          << " AsmblAlign " << record._svObs.assemblyAlignIndex
          << " AsmblSeg " << record._svObs.assemblySegmentIndex;

    return ostrm;
}

/*****************************************************************************/
// EdgeInfoRecord
/*****************************************************************************/

EdgeInfoRecord::EdgeInfoRecord()
{
    _edgeInfo = EdgeInfo();
    _discardReason = KEPT;
}

/*****************************************************************************/

EdgeInfoRecord::EdgeInfoRecord(const EdgeInfo& edgeInfo)
{
    _edgeInfo = edgeInfo;
    _discardReason = KEPT;
}

/*****************************************************************************/

const EdgeInfo& EdgeInfoRecord::getEdgeInfo() const
{
    return _edgeInfo;
}

/*****************************************************************************/

bool EdgeInfoRecord::discard(DiscardReason discardReason)
{
    if (_discardReason != KEPT)
    {
        return false;
    }

    _discardReason = discardReason;
    return true;
}

/*****************************************************************************/

bool EdgeInfoRecord::discarded(DiscardReason& discardReason) const
{
    discardReason = _discardReason;
    return (_discardReason != KEPT);
}

/*****************************************************************************/

bool EdgeInfoRecord::operator<(const EdgeInfoRecord& edgeInfoRecordB) const
{
    if (_edgeInfo.locusIndex < edgeInfoRecordB._edgeInfo.locusIndex)
    {
        return true;
    }

    if (_edgeInfo.nodeIndex1 < edgeInfoRecordB._edgeInfo.nodeIndex1)
    {
        return true;
    }

    if (_edgeInfo.nodeIndex2 < edgeInfoRecordB._edgeInfo.nodeIndex2)
    {
        return true;
    }

    return false;
}

/*****************************************************************************/
// Based only on the EdgeInfo

bool EdgeInfoRecord::operator==(const EdgeInfoRecord& edgeInfoRecordB) const
{
    return ((_edgeInfo.locusIndex == edgeInfoRecordB._edgeInfo.locusIndex)
            && (_edgeInfo.nodeIndex1 == edgeInfoRecordB._edgeInfo.nodeIndex1)
            && (_edgeInfo.nodeIndex2 == edgeInfoRecordB._edgeInfo.nodeIndex2));
}

/*****************************************************************************/

std::ostream& operator<<(std::ostream& ostrm, const EdgeInfoRecord& record)
{
    ostrm << "Locus " << record._edgeInfo.locusIndex
          << " Node " << record._edgeInfo.nodeIndex1
          << " Node " << record._edgeInfo.nodeIndex2;

    return ostrm;
}

/*****************************************************************************/



/*****************************************************************************/
// SVLog
/*****************************************************************************/

SVLog::SVLog(unsigned int index)
    : _ind(index), _outcome(UNKNOWN)
{
    ;
}

/*****************************************************************************/

unsigned int SVLog::ind() const
{
    return _ind;
}

/*****************************************************************************/

void SVLog::reportOutcome(const Outcome outcomeVal)
{
    assert(outcomeVal != SVLog::UNKNOWN);
    assert(outcomeVal != SVLog::SPAWNED);

    assert(_outcome == UNKNOWN);
    _outcome = outcomeVal;
}

/*****************************************************************************/

SVLog::Outcome SVLog::outcome() const
{
    return _outcome;
}

/*****************************************************************************/

std::ostream& operator<<(std::ostream& ostrm, const SVLog::Outcome outcome)
{
    switch (outcome)
    {
    case SVLog::UNKNOWN:
        ostrm << "UNKNOWN";
        break;
    case SVLog::WRITTEN:
        ostrm << "Written";
        break;
    case SVLog::SPAWNED:
        ostrm << "Spawned";
        break;
    case SVLog::IMPRECISE_NON_SPANNING:
        ostrm << "Imprecise_NonSpanning";
        break;
    case SVLog::LOW_SPANNING_COUNT:
        ostrm << "LowSpanningCount";
        break;
    case SVLog::LOW_SOMATIC_SCORE:
        ostrm << "LowSomaticScore";
        break;
    }

    return ostrm;
}

/*****************************************************************************/
// AssembledSVLog
/*****************************************************************************/

AssembledSVLog::AssembledSVLog(unsigned int index)
    : SVLog(index)
{
    ;
}

/*****************************************************************************/

std::ostream& operator<<(std::ostream& ostrm, const AssembledSVLog& log)
{
    ostrm << "    assembledSV " << log.ind() << " " << log.outcome();
    return ostrm;
}

/*****************************************************************************/
// CandSVLog
/*****************************************************************************/

CandSVLog::CandSVLog(unsigned int index)
    : SVLog(index), _expectedNumAssembled(0), _numWrittenAssembled(0)
{
    ;
}

/*****************************************************************************/

void CandSVLog::reportNumAssembled(unsigned int numAssembled)
{
    assert(_expectedNumAssembled == 0);
    _expectedNumAssembled = numAssembled;
}

/*****************************************************************************/

void CandSVLog::addAssembledSV()
{
    if (_outcome == SVLog::UNKNOWN)
    {
        _outcome = SVLog::SPAWNED;
    }

    assert(_outcome == SVLog::SPAWNED);
    assert(_assembledSVLogVec.size() < _expectedNumAssembled);


    AssembledSVLog assembledSVLog(_assembledSVLogVec.size());
    _assembledSVLogVec.push_back(assembledSVLog);
}

/*****************************************************************************/

void CandSVLog::reportOutcome(const SVLog::Outcome outcomeVal)
{
    assert(outcomeVal != SVLog::UNKNOWN);
    assert(outcomeVal != SVLog::SPAWNED);

    if (_assembledSVLogVec.empty())
    {
        // No spawned assembled SVs, so outcome is for this candidate SV
        assert(_outcome == SVLog::UNKNOWN);
        _outcome = outcomeVal;
    }
    else
    {
        // Outcome is forwarded to most recently added assembled SV.
        assert(_outcome == SVLog::SPAWNED);

        _assembledSVLogVec.back().reportOutcome(outcomeVal);

        if (outcomeVal == SVLog::WRITTEN)
        {
            ++_numWrittenAssembled;
        }
    }
}

/*****************************************************************************/

unsigned int CandSVLog::numAssembledSVs() const
{
    return _assembledSVLogVec.size();
}

/*****************************************************************************/

unsigned int CandSVLog::numWrittenAssembledSVs() const
{
    return _numWrittenAssembled;
}

/*****************************************************************************/

std::ostream& operator<<(std::ostream& ostrm, const CandSVLog& log)
{
    ostrm << "  candSV " << log.ind()
          << " : " << log.numAssembledSVs() << " assembled SVs"
          << " : " << log.outcome() << std::endl;

    BOOST_FOREACH(const AssembledSVLog& assembledLog, log._assembledSVLogVec)
    {
        ostrm << assembledLog << std::endl;
    }

    return ostrm;
}

/*****************************************************************************/
// EdgeLog
/*****************************************************************************/

EdgeLog::EdgeLog(unsigned int index, const EdgeInfo& edgeInfo)
    : EdgeInfoRecord(edgeInfo), SVLog(index),
      _expectedNumCands(0), _numWrittenCands(0)
{
    ;
}

/*****************************************************************************/

void EdgeLog::addVariantKey(VariantKey variantKey)
{
    _variantKeyVec.push_back(variantKey);
}

/*****************************************************************************/

const EdgeLog::VariantKeyVec& EdgeLog::variantKeyVec() const
{
    return _variantKeyVec;
}

/*****************************************************************************/

void EdgeLog::reportNumCands(unsigned int numCands)
{
    assert(_expectedNumCands == 0);
    _expectedNumCands = numCands;
}

/*****************************************************************************/

void EdgeLog::addCandSV()
{
    assert(_candSVLogVec.size() < _expectedNumCands);

    CandSVLog candSVLog(_candSVLogVec.size());
    _candSVLogVec.push_back(candSVLog);
}

/*****************************************************************************/
// Just pass request to most recently added CandSVLog

void EdgeLog::reportNumAssembled(unsigned int numAssembled)
{
    _candSVLogVec.back().reportNumAssembled(numAssembled);
}

/*****************************************************************************/
// Just pass request to most recently added CandSVLog

void EdgeLog::addAssembledSV()
{
    _candSVLogVec.back().addAssembledSV();
}

/*****************************************************************************/

void EdgeLog::reportOutcome(const Outcome outcomeVal)
{
    assert(outcomeVal != SVLog::UNKNOWN);
    assert(outcomeVal != SVLog::SPAWNED);

    assert(!_candSVLogVec.empty());

    _candSVLogVec.back().reportOutcome(outcomeVal);
}

/*****************************************************************************/

unsigned int EdgeLog::numCandSVs() const
{
    return _candSVLogVec.size();
}

/*****************************************************************************/

unsigned int EdgeLog::numWrittenCandSVs() const
{
    unsigned int numWrittenCands(0);

    BOOST_FOREACH(const CandSVLog& candSVLog, _candSVLogVec)
    {
        if (candSVLog.numAssembledSVs() > 0)
        {
            if (candSVLog.numWrittenAssembledSVs() > 0)
            {
                ++numWrittenCands; // assembled SVs from candidate written
            }
        }
        else
        {
            if (candSVLog.outcome() == SVLog::WRITTEN)
            {
                ++numWrittenCands; // candidate written
            }
        }
    }

    return numWrittenCands;
}

/*****************************************************************************/

std::ostream& operator<<(std::ostream& ostrm, const EdgeLog& log)
{
    ostrm << "Edge " << log.ind() << " : "
          << *((EdgeInfo*) &log) << " [";

    BOOST_FOREACH(const EdgeLog::VariantKey variantKey,  log._variantKeyVec)
    {
        ostrm << " " << variantKey;
    }

    ostrm << "] numTrueMatched " << log._variantKeyVec.size()
          << " : " << log._candSVLogVec.size() << " cand SVs" << std::endl;

    BOOST_FOREACH(const CandSVLog& candSVLog, log._candSVLogVec)
    {
        ostrm << candSVLog;
    }

    return ostrm;
}

/*****************************************************************************/
// TruthTracker
/*****************************************************************************/

TruthTracker::TruthTracker(const std::string& vcfFilePath,
                           const bam_header_info& bamHeaderInfo)
    : _cset(_emptyLocusSet), _hasTruth(false), _numObs(0),
      _lastAddedEdgeLogMapEleIter(0), _numEdgesSeen(0)
{
    if (vcfFilePath.empty()) return;

    loadTruth(vcfFilePath, bamHeaderInfo);
}

/*****************************************************************************/


TruthTracker::
TruthTracker(
    const std::string& vcfFilePath,
    const SVLocusSet& cset)
    : _cset(cset), _hasTruth(false), _numObs(0),
      _lastAddedEdgeLogMapEleIter(0), _numEdgesSeen(0)
{
    if (vcfFilePath.empty()) return;

    loadTruth(vcfFilePath, cset.header);
#ifdef EASY_ITER_OVER_NODE_EDGES
    evalLocusSet(cset);
#endif
}

/*****************************************************************************/

bool TruthTracker::loadTruth(const std::string& vcfFilePath,
                             const bam_header_info& bamHeaderInfo)
{
    _chromNameTidMap = bamHeaderInfo.chrom_to_index;

    VcfFile vcfFile(vcfFilePath, _chromNameTidMap);
    vcfFile.getVariantVec(_trueVariantVec);

    _hasTruth = true;

    return true;
}

/*****************************************************************************/

static
bool bothBrkptIntersect(const Variant& variant, const SVObservation& svObs)
{
    return ((variant.brkptA().isIntersect(svObs.bp1.interval)
             && variant.brkptB().isIntersect(svObs.bp2.interval))
            ||
            (variant.brkptA().isIntersect(svObs.bp2.interval)
             && variant.brkptB().isIntersect(svObs.bp1.interval)));
}

/*****************************************************************************/

void TruthTracker::addObservation(const SVObservation& svObs)
{
    if (!_hasTruth)
    {
        return;
    }

    ++_numObs;

    EdgeLog::VariantKey variantKey(0);
    const ObsRecord obsRecord(svObs);

    // DEBUG
    // std::cerr << "addObservation : " << svObs << std::endl;

    BOOST_FOREACH(const Variant& trueVariant, _trueVariantVec)
    {
        if (bothBrkptIntersect(trueVariant, svObs))
        {
            _truthObsVecMap[variantKey].push_back(obsRecord);

            // DEBUG
            std::cerr << "Variant " << variantKey << " [" << trueVariant
                      << "] matched by " << obsRecord
                      << " bp1 " << svObs.bp1.interval
                      << " bp2 " << svObs.bp2.interval
                      << std::endl;


#if 0

            if (_lastAddedObsLogMapEleIter == _obsLogMap.end()) // not set
            {
                ObsLogMapIterBoolPr iterBoolPr
                    = _obsLogMap.
                    insert(ObsLogMapEle(edgeInfoRecord,
                                         EdgeLog(_edgeLogMap.size(), edge)));
                assert(iterBoolPr.second); // edge should be unique -> added
                _lastAddedObsLogMapEleIter = iterBoolPr.first;
                _obsLogMapIndexVec.push_back(_lastAddedObsLogMapEleIter);
            }

            _lastAddedObsLogMapEleIter->second.addVariantKey(variantKey);
#endif

        }

        ++variantKey;
    }
}

/*****************************************************************************/

#ifdef EASY_ITER_OVER_NODE_EDGES

void dumpLocusSetStats(const SVLocusSet& svLocusSet)
{
    std::cerr << "Locus Set Structure & Stats" << std::endl;

    std::cerr << "Num Loci : " << svLocusSet.size() << std::endl;

    typedef std::vector<unsigned int> FreqTable;
    typedef FreqTable::value_type FreqTableEle;

    typedef SVLocusEdgesType::value_type NodeEdgeEle;

    FreqTable numLocusNodesFreqs;

    unsigned int totalNumEdges(0);

    BOOST_FOREACH(const SVLocus& locus, svLocusSet)
    {
        const unsigned int numLocusNodes(locus.size());

        if (numLocusNodes >= numLocusNodesFreqs.size())
        {
            numLocusNodesFreqs.resize(numLocusNodes + 1);
        }

        ++numLocusNodesFreqs[numLocusNodes];

        std::cerr << "Locus[" << locus.getIndex() << "] "
                  << locus.size() << " Nodes" << std::endl;

        unsigned int nodeInd(0);

        BOOST_FOREACH(const SVLocusNode& node, locus)
        {
            std::cerr << "  Node[ " << nodeInd << "] "
                      << node.size() << " Edges [to other Nodes";
            totalNumEdges += node.size();

            // FIXME : Currently Node is no longer iterable and would need
            //         to use EdgeManager to restore functionality
            //         but this may be a temporary situation.
            //         For now commenting out everything that relies directly
            //         or indirectly on iterable Node with an
            //         EASY_ITER_OVER_NODE_EDGES ifdef.
            // const SVLocusEdgeManager nodeManager(node.getEdgeManager());
            // nodeManager.getMap()
            BOOST_FOREACH(const NodeEdgeEle& nodeEdge, node)
            {
                std::cerr << " " << nodeEdge.first;
            }

            std::cerr << "]" << std::endl;

            ++nodeInd;
        }
    }

    std::cerr << std::endl;

    unsigned int numLocusNodes(0);
    std::cerr << "Nodes per Locus:Freq";

    BOOST_FOREACH(const unsigned int freq, numLocusNodesFreqs)
    {
        if (freq > 0)
        {
            std::cerr << " " << numLocusNodes << ":" << freq;
        }

        ++numLocusNodes;
    }

    std::cerr << std::endl;

    std::cerr << "totalNumEdges " << totalNumEdges << std::endl
              << std::endl;
}

#endif

/*****************************************************************************/

typedef std::vector<SVLocusSet::NodeAddressType> NodeAddrVec;
typedef NodeAddrVec::const_iterator NodeAddrVecCIter;

/*****************************************************************************/
#ifdef EASY_ITER_OVER_NODE_EDGES
static
bool getMatchedNodes(const SVLocusSet& svLocusSet,
                     const GenomeInterval& interval,
                     NodeAddrVec& matchedNodeAddrVec)
{
    BOOST_FOREACH(const SVLocus& locus, svLocusSet)
    {
        NodeIndexType nodeIndex(0);

        BOOST_FOREACH(const SVLocusNode& node, locus)
        {
            if (interval.isIntersect(node.getInterval()))
            {
                SVLocusSet::NodeAddressType nodeAddr(locus.getIndex(),
                                                     nodeIndex);
                matchedNodeAddrVec.push_back(nodeAddr);
            }

            ++nodeIndex;
        }
    }

    return !matchedNodeAddrVec.empty();
}
#endif


/*****************************************************************************/

#ifdef EASY_ITER_OVER_NODE_EDGES

bool hasABLink(const SVLocus& svLocus,
               const NodeIndexType nodeIndA, const NodeIndexType nodeIndB)
{
    typedef SVLocusEdgesType::value_type NodeEdgeEle;
    const SVLocusNode& nodeA(svLocus.getNode(nodeIndA));

    BOOST_FOREACH(const NodeEdgeEle& nodeEdge, nodeA)
    {
        if (nodeEdge.first == nodeIndB)
        {
            return true;
        }
    }

    return false;
}

#endif

/*****************************************************************************/

#ifdef EASY_ITER_OVER_NODE_EDGES

bool pairBrkptNodes(const SVLocusSet& svLocusSet,
                    const NodeAddrVec& matchedNodeAddrVecA,
                    const NodeAddrVec& matchedNodeAddrVecB,
                    unsigned int& numMultiNodes,
                    unsigned int& numDiffLoci,
                    unsigned int& numTwoNodesSameLocus,
                    unsigned int& numTwoNodesSameLocusWithBiEdge,
                    unsigned int& numTwoNodesSameLocusWithABEdgeOnly,
                    unsigned int& numTwoNodesSameLocusWithBAEdgeOnly,
                    unsigned int& numSameNode)
{
    assert((!matchedNodeAddrVecA.empty())
           && (!matchedNodeAddrVecB.empty()));

    if ((matchedNodeAddrVecA.size() > 1) || (matchedNodeAddrVecB.size() > 1))
    {
        ++numMultiNodes;
        std::cerr << " MultiNodes";
        return false;
    }

    const SVLocusSet::NodeAddressType& nodeAddrA(matchedNodeAddrVecA[0]);
    const SVLocusSet::NodeAddressType& nodeAddrB(matchedNodeAddrVecB[0]);

    if (nodeAddrB.first != nodeAddrA.first)
    {
        ++numDiffLoci;
        std::cerr << " DiffLoci";
    }
    else if (nodeAddrB.second != nodeAddrA.second)
    {
        ++numTwoNodesSameLocus;

        const SVLocus& svLocus(svLocusSet.getLocus(nodeAddrA.first));
        const bool foundBFromA(hasABLink(svLocus, nodeAddrA.second,
                                         nodeAddrB.second));
        const bool foundAFromB(hasABLink(svLocus, nodeAddrB.second,
                                         nodeAddrA.second));

        if (foundBFromA && foundAFromB)
        {
            ++numTwoNodesSameLocusWithBiEdge;
            std::cerr << " TwoNodesSameLocusWithBiEdge";
            return true;
        }
        else if (foundBFromA)
        {
            ++numTwoNodesSameLocusWithABEdgeOnly;
            std::cerr << " TwoNodesSameLocusWithABEdgeOnly";
        }
        else if (foundAFromB)
        {
            ++numTwoNodesSameLocusWithBAEdgeOnly;
            std::cerr << " TwoNodesSameLocusWithBAEdgeOnly";
        }
    }
    else
    {
        ++numSameNode;
        std::cerr << " SameNode";
    }

    return false;
}

#endif

/*****************************************************************************/

#ifdef EASY_ITER_OVER_NODE_EDGES

bool TruthTracker::evalLocusSet(const SVLocusSet& svLocusSet)
{
    if (!_hasTruth)
    {
        return false;
    }

    // DEBUG
    dumpLocusSetStats(svLocusSet);

    std::cerr << "Truth variant match against Locus Set" << std::endl;

    unsigned int numNeither(0); // Neither brkpt matched by Node
    unsigned int numOneBrkpt(0); // Only one brkpt matched by Node
    unsigned int numMultiNodes(0); // Multiple Nodes for one/both brkpts
    unsigned int numDiffLoci(0); // Brkpts matched by Nodes in different Loci
    unsigned int numTwoNodesSameLocus(0); // " " Nodes in same Locus

    unsigned int numTwoNodesSameLocusWithBiEdge(0);
    unsigned int numTwoNodesSameLocusWithABEdgeOnly(0);
    unsigned int numTwoNodesSameLocusWithBAEdgeOnly(0);

    unsigned int numSameNode(0); // Brkpts matched by same Node

    EdgeLog ::VariantKey variantKey(0);

    BOOST_FOREACH(const Variant& trueVariant, _trueVariantVec)
    {
        NodeAddrVec matchedNodeAddrVecA;
        const bool intersectA(getMatchedNodes(svLocusSet,
                                              trueVariant.brkptA(),
                                              matchedNodeAddrVecA));

        NodeAddrVec matchedNodeAddrVecB;
        const bool intersectB(getMatchedNodes(svLocusSet,
                                              trueVariant.brkptB(),
                                              matchedNodeAddrVecB));

        std::cerr << variantKey << " : " << trueVariant << " : "
                  << " A intersects " << matchedNodeAddrVecA.size() << "[";
        BOOST_FOREACH(const SVLocusSet::NodeAddressType& nodeAddr,
                      matchedNodeAddrVecA)
        {
            std::cerr << " " << nodeAddr.first << ":" << nodeAddr.second;
        }

        std::cerr << "], B intersects " << matchedNodeAddrVecB.size() << "[";
        BOOST_FOREACH(const SVLocusSet::NodeAddressType& nodeAddr,
                      matchedNodeAddrVecB)
        {
            std::cerr << " " << nodeAddr.first << ":" << nodeAddr.second;
        }

        std::cerr << "]";

        if (intersectA && intersectB)
        {
            pairBrkptNodes(svLocusSet,
                           matchedNodeAddrVecA, matchedNodeAddrVecB,
                           numMultiNodes, numDiffLoci,
                           numTwoNodesSameLocus,
                           numTwoNodesSameLocusWithBiEdge,
                           numTwoNodesSameLocusWithABEdgeOnly,
                           numTwoNodesSameLocusWithBAEdgeOnly,
                           numSameNode);
        }
        else if (intersectA || intersectB)
        {
            ++numOneBrkpt;
        }
        else
        {
            ++numNeither;
        }

        std::cerr << std::endl;
        ++variantKey;
    }

    std::cerr << std::endl;
    std::cerr << "numNeither " << numNeither << std::endl;
    std::cerr << "numOneBrkpt " << numOneBrkpt << std::endl;
    std::cerr << "numMultiNodes " << numMultiNodes << std::endl;
    std::cerr << "numDiffLoci " << numDiffLoci << std::endl;
    std::cerr << "numTwoNodesSameLocus " << numTwoNodesSameLocus << std::endl;
    std::cerr << "numTwoNodesSameLocusWithBiEdge "
              << numTwoNodesSameLocusWithBiEdge << std::endl;
    std::cerr << "numTwoNodesSameLocusWithABEdgeOnly "
              << numTwoNodesSameLocusWithABEdgeOnly << std::endl;
    std::cerr << "numTwoNodesSameLocusWithBAEdgeOnly "
              << numTwoNodesSameLocusWithBAEdgeOnly << std::endl;
    std::cerr << "numSameNode " << numSameNode << std::endl;
    std::cerr << std::endl;

    return true;
}

#endif

/*****************************************************************************/
static
bool bothBrkptIntersect(const Variant& variant,
                        const SVLocusNode& node1, const SVLocusNode& node2)
{
    return ((variant.brkptA().isIntersect(node1.getInterval())
             && variant.brkptB().isIntersect(node2.getInterval()))
            ||
            (variant.brkptA().isIntersect(node2.getInterval())
             && variant.brkptB().isIntersect(node1.getInterval())));
}

/*****************************************************************************/

bool TruthTracker::addEdge(const EdgeInfo& edge)
{
    if (!_hasTruth)
    {
        return false;
    }

    ++_numEdgesSeen;

    const SVLocus& locus(_cset.getLocus(edge.locusIndex));
    const SVLocusNode& node1(locus.getNode(edge.nodeIndex1));
    const SVLocusNode& node2(locus.getNode(edge.nodeIndex2));

    EdgeLog::VariantKey variantKey(0);
    const EdgeInfoRecord edgeInfoRecord(edge);

    // This iter will become valid only if the Edge intersects truth.
    _lastAddedEdgeLogMapEleIter = _edgeLogMap.end(); // not set

    BOOST_FOREACH(const Variant& trueVariant, _trueVariantVec)
    {
        if (bothBrkptIntersect(trueVariant, node1, node2))
        {

            _truthEdgeVecMap[variantKey].push_back(edgeInfoRecord);

            if (_lastAddedEdgeLogMapEleIter == _edgeLogMap.end()) // not set
            {
                EdgeLogMapIterBoolPr iterBoolPr
                    = _edgeLogMap.
                      insert(EdgeLogMapEle(edgeInfoRecord,
                                           EdgeLog(_edgeLogMap.size(), edge)));
                assert(iterBoolPr.second); // edge should be unique -> added
                _lastAddedEdgeLogMapEleIter = iterBoolPr.first;
                _edgeLogMapIndexVec.push_back(_lastAddedEdgeLogMapEleIter);
            }

            _lastAddedEdgeLogMapEleIter->second.addVariantKey(variantKey);
        }

        ++variantKey;
    }

    return true;
}

/*****************************************************************************/

bool TruthTracker::discardEdge(const EdgeInfo& edge,
                               EdgeInfoRecord::DiscardReason reason)
{
    if (!_hasTruth)
    {
        return false;
    }

    const EdgeInfoRecord discardedEdgeInfoRecord(edge);

    EdgeLogMapIter
    edgeLogIter(_edgeLogMap.find(discardedEdgeInfoRecord));

    if (edgeLogIter != _edgeLogMap.end())
    {
        // First find all the true variants matched by this Edge.
        const EdgeLog::VariantKeyVec&
        variantKeyVec(edgeLogIter->second.variantKeyVec());

        BOOST_FOREACH(const EdgeLog::VariantKey variantKey, variantKeyVec)
        {
            BOOST_FOREACH(EdgeInfoRecord edgeInfoRecord,
                          _truthEdgeVecMap[variantKey])
            {
                if (edgeInfoRecord == discardedEdgeInfoRecord)
                {
                    edgeInfoRecord.discard(reason);
                }
            }
        }
    }

    return true;
}

/*****************************************************************************/
// Just pass request to most recently added EdgeLog

void TruthTracker::reportNumCands(unsigned int numCands)
{
    if (!_hasTruth)
    {
        return;
    }

    if (_lastAddedEdgeLogMapEleIter != _edgeLogMap.end())
    {
        _lastAddedEdgeLogMapEleIter->second.reportNumCands(numCands);
    }
}

/*****************************************************************************/
// Just pass request to most recently added EdgeLog

void TruthTracker::addCandSV()
{
    if (!_hasTruth)
    {
        return;
    }

    if (_lastAddedEdgeLogMapEleIter != _edgeLogMap.end())
    {
        _lastAddedEdgeLogMapEleIter->second.addCandSV();
    }
}

/*****************************************************************************/
// Just pass request to most recently added EdgeLog -> most recent CandSVLog

void TruthTracker::reportNumAssembled(unsigned int numAssembled)
{
    if (!_hasTruth)
    {
        return;
    }

    if (_lastAddedEdgeLogMapEleIter != _edgeLogMap.end())
    {
        _lastAddedEdgeLogMapEleIter->second.reportNumAssembled(numAssembled);
    }
}

/*****************************************************************************/
// Just pass request to most recently added EdgeLog -> most recent CandSVLog

void TruthTracker::addAssembledSV()
{
    if (!_hasTruth)
    {
        return;
    }

    if (_lastAddedEdgeLogMapEleIter != _edgeLogMap.end())
    {
        _lastAddedEdgeLogMapEleIter->second.addAssembledSV();
    }
}

/*****************************************************************************/
// Just pass request to most recently added EdgeLog

void TruthTracker::reportOutcome(const SVLog::Outcome outcomeVal)
{
    if (!_hasTruth)
    {
        return;
    }

    if (_lastAddedEdgeLogMapEleIter != _edgeLogMap.end())
    {
        _lastAddedEdgeLogMapEleIter->second.reportOutcome(outcomeVal);
    }
}

/*****************************************************************************/

bool TruthTracker::dumpObs()
{
    if (!_hasTruth)
    {
        return false;
    }

    std::cerr << "Num Obs added : " << _numObs << std::endl;


    return true;
}

/*****************************************************************************/

bool TruthTracker::dumpMatchedEdgeInfo()
{
    if (!_hasTruth)
    {
        return false;
    }

    BOOST_FOREACH(const EdgeLogMapIter& edgeLogMapEleIter, _edgeLogMapIndexVec)
    {
        std::cerr << edgeLogMapEleIter->second;
    }

    return true;
}

/*****************************************************************************/

bool TruthTracker::dumpStats()
{
    if (!_hasTruth)
    {
        return false;
    }

    std::cerr << "Num edges seen : " << _numEdgesSeen << std::endl;

    unsigned int numTPs(0);

    BOOST_FOREACH(const TruthEdgeVecMapEle& mapEle, _truthEdgeVecMap)
    {
        bool detected(false);

        BOOST_FOREACH(const EdgeInfoRecord& record, mapEle.second)
        {
            EdgeInfoRecord::DiscardReason discardReason(EdgeInfoRecord::KEPT);

            if (!record.discarded(discardReason))
            {
                detected = true;
                break;
            }
        }

        if (detected)
        {
            ++numTPs;
        }
    }

    std::cerr << "Num true SVs : " << _trueVariantVec.size() << std::endl;
    std::cerr << "Num orig detections : " << _truthEdgeVecMap.size()
              << std::endl;
    std::cerr << "Num TPs : " << numTPs << std::endl;

    return true;
}

/*****************************************************************************/



