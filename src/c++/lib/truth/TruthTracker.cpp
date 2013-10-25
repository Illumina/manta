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

#include "format/VcfFile.hh"

#include "truth/TruthTracker.hh"

/*****************************************************************************/
// EdgeInfoRecord
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
    if (_edgeInfo.locusIndex < edgeInfoRecordB._edgeInfo.locusIndex) {
        return true;
    }

    if (_edgeInfo.nodeIndex1 < edgeInfoRecordB._edgeInfo.nodeIndex1) {
        return true;
    }

    if (_edgeInfo.nodeIndex2 < edgeInfoRecordB._edgeInfo.nodeIndex2) {
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
// TruthTracker
/*****************************************************************************/

TruthTracker::TruthTracker()
    : _hasTruth(false), _numEdgesSeen(0)
{
    ;
}

/*****************************************************************************/

bool TruthTracker::loadTruth(const std::string& vcfFilePath,
                             const bam_header_info& bamHeaderInfo)
{
    make_chrom_tid_map(bamHeaderInfo, _chromNameTidMap);

    VcfFile vcfFile(vcfFilePath, _chromNameTidMap);
    vcfFile.getVariantVec(_trueVariantVec);

    _hasTruth = true;
    _vcfFilePath = vcfFilePath;

    return true;
}

/*****************************************************************************/

void dumpLocusSetStats(const SVLocusSet& svLocusSet)
{
    std::cerr << "Locus Set Structure & Stats" << std::endl;

    std::cerr << "Num Loci : " << svLocusSet.size() << std::endl;

    typedef std::vector<unsigned int> FreqTable;
    typedef FreqTable::value_type FreqTableEle;

    typedef SVLocusNode::edges_type::value_type NodeEdgeEle;

    FreqTable numLocusNodesFreqs;

    unsigned int totalNumEdges(0);

    BOOST_FOREACH(const SVLocus& locus, svLocusSet)
    {
        const unsigned int numLocusNodes(locus.size());

        if (numLocusNodes >= numLocusNodesFreqs.size()) {
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

/*****************************************************************************/

typedef std::vector<SVLocusSet::NodeAddressType> NodeAddrVec;
typedef NodeAddrVec::const_iterator NodeAddrVecCIter;

/*****************************************************************************/

bool getMatchedNodes(const SVLocusSet& svLocusSet,
                     const GenomeInterval& interval,
                     NodeAddrVec& matchedNodeAddrVec)
{
    BOOST_FOREACH(const SVLocus& locus, svLocusSet)
    {
        NodeIndexType nodeIndex(0);

        BOOST_FOREACH(const SVLocusNode& node, locus)
        {
            if (interval.isIntersect(node.getInterval())) {
                SVLocusSet::NodeAddressType nodeAddr(locus.getIndex(),
                                                     nodeIndex);
                matchedNodeAddrVec.push_back(nodeAddr);
            }

            ++nodeIndex;
        }
    }

    return !matchedNodeAddrVec.empty();
}

/*****************************************************************************/

bool hasABLink(const SVLocus& svLocus,
               const NodeIndexType nodeIndA, const NodeIndexType nodeIndB)
{
    typedef SVLocusNode::edges_type::value_type NodeEdgeEle;
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

/*****************************************************************************/

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
    assert((matchedNodeAddrVecA.size() > 0)
           && (matchedNodeAddrVecB.size() > 0));

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

/*****************************************************************************/

bool TruthTracker::evalLocusSet(const SVLocusSet& svLocusSet)
{
    if (!_hasTruth) {
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

    VariantKey variantKey(0);

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

/*****************************************************************************/

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

bool TruthTracker::addEdge(const EdgeInfo& edge, const SVLocusSet& cset)
{
    if (!_hasTruth) {
        return false;
    }

    ++_numEdgesSeen;

    const SVLocus& locus(cset.getLocus(edge.locusIndex));
    const SVLocusNode& node1(locus.getNode(edge.nodeIndex1));
    const SVLocusNode& node2(locus.getNode(edge.nodeIndex2));

    // DEBUG
    std::cerr << "Edge : Locus " << edge.locusIndex
              << " Node " << edge.nodeIndex1
              << " Node " << edge.nodeIndex2
              << " [";
    unsigned int numTrueMatched(0);

    VariantKey variantKey(0);

    BOOST_FOREACH(const Variant& trueVariant, _trueVariantVec)
    {
        if (bothBrkptIntersect(trueVariant, node1, node2))
        {
            const EdgeInfoRecord edgeInfoRecord(edge);

            _truthEdgeVecMap[variantKey].push_back(edgeInfoRecord);
            _edgeTruthVecMap[edgeInfoRecord].push_back(variantKey);

            ++numTrueMatched;
            
            // DEBUG
            std::cerr << " " << variantKey;
        }

        ++variantKey;
    }

    // DEBUG
    std::cerr << "] numTrueMatched " << numTrueMatched << std::endl;

    return true;
}

/*****************************************************************************/

bool TruthTracker::discardEdge(const EdgeInfo& edge,
                               EdgeInfoRecord::DiscardReason reason)
{
    const EdgeInfoRecord discardedEdgeInfoRecord(edge);

    EdgeTruthVecMapIter
        edgeTruthVecIter(_edgeTruthVecMap.find(discardedEdgeInfoRecord));

    if (edgeTruthVecIter != _edgeTruthVecMap.end())
    {
        // First find all the true variants matched by this Edge.
        const VariantKeyVec& variantKeyVec(edgeTruthVecIter->second);

        BOOST_FOREACH(const VariantKey variantKey, variantKeyVec)
        {
            BOOST_FOREACH(EdgeInfoRecord edgeInfoRecord, 
                          _truthEdgeVecMap[variantKey])
            {
                if (edgeInfoRecord == discardedEdgeInfoRecord) {
                    edgeInfoRecord.discard(reason);
                }
            }
        }
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

            if (!record.discarded(discardReason)) {
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



