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
// <https://github.com/downloads/sequencing/licenses/>.
//

///
/// \author Chris Saunders
///

#include "blt_util/log.hh"
#include "common/Exceptions.hh"
#include "manta/SVLocusSet.hh"

#include "boost/archive/binary_iarchive.hpp"
#include "boost/archive/binary_oarchive.hpp"
#include "boost/foreach.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>



std::ostream&
operator<<(std::ostream& os, const SVLocusSet::NodeAddressType& a)
{
    os << a.first << ":" << a.second;
    return os;
}



void
SVLocusSet::
merge(const SVLocus& inputLocus)
{
    //
    // test each node in input locus for intersection and insert/join to existing node
    //

    using namespace illumina::common;

    assert(! _isFinalized);

#ifdef DEBUG_SVL
    checkState(true);
    log_os << "SVLocusSet::merge inputLocus: " << inputLocus;
#endif

    const LocusIndexType startLocusIndex(insertLocus(inputLocus));
    const SVLocus& startLocus(_loci[startLocusIndex]);
    LocusIndexType headLocusIndex(startLocusIndex);
    bool isInputLocusMoved(false);

    // because we have a non-general interval overlap test, we must order search
    // nodes by begin_pos on each chromosome
    //
    typedef std::map<GenomeInterval,NodeIndexType> nodeMap_t;
    nodeMap_t nodeMap;
    {
        const NodeIndexType nodeCount(startLocus.size());
        for (NodeIndexType nodeIndex(0); nodeIndex<nodeCount; ++nodeIndex)
        {
            nodeMap.insert(std::make_pair(startLocus.getNode(nodeIndex).interval,nodeIndex));
        }
    }

    BOOST_FOREACH(const nodeMap_t::value_type& nodeVal, nodeMap)
    {
        const NodeIndexType nodeIndex(nodeVal.second);

#ifdef DEBUG_SVL
        log_os << "SVLocusSet::merge inputNode: " << NodeAddressType(std::make_pair(startLocusIndex,nodeIndex)) << " " << startLocus.getNode(nodeIndex);
#endif

        std::set<NodeAddressType> intersectNodes;
        getNodeMergeableIntersect(startLocusIndex, nodeIndex, isInputLocusMoved, intersectNodes);

#ifdef DEBUG_SVL
        log_os << "SVLocusSet::merge insersect_size: " << intersectNodes.size() << "\n";
        BOOST_FOREACH(const NodeAddressType& val, intersectNodes)
        {
            log_os << "i-index: " << val << " node: " <<  getNode(val) << "\n";
        }
#endif

        if (isInputLocusMoved)
        {
#if 0
            if (intersectNodes.empty())
            {
                const NodeAddressType val(std::make_pair(startLocusIndex,nodeIndex));
                std::ostringstream oss;
                oss << "ERROR: no intersecting nodes found during merge\n"
                    << "\tsearch node: " << val << " " << getNode(val)
                    << "\thli: " << headLocusIndex << "\n";
                BOOST_THROW_EXCEPTION(PreConditionException(oss.str()));
            }
#endif
            if (2>intersectNodes.size()) continue;
        }
        else
        {
            if (intersectNodes.empty()) continue;
        }

        bool isMultiLocus(false);
        BOOST_FOREACH(const NodeAddressType& addy, intersectNodes)
        {
            if(addy.first == headLocusIndex) continue;
            isMultiLocus=true;
            break;
        }

        if(isMultiLocus)
        {
            // if there are any intersections, copy the loci of all intersecting nodes into
            // a single locus, by convention we use the lowest locusIndex of the intersecting set
            moveIntersectToLowIndex(intersectNodes,startLocusIndex,headLocusIndex);
            if(! isInputLocusMoved) isInputLocusMoved=(headLocusIndex != startLocusIndex);

            getNodeMergeableIntersect(startLocusIndex, nodeIndex, isInputLocusMoved, intersectNodes);
            assert(! intersectNodes.empty());
        }

#ifdef DEBUG_SVL
        log_os << "intersect2_size: " << intersectNodes.size() << "\n";
        BOOST_FOREACH(const NodeAddressType& val, intersectNodes)
        {
            log_os << "i2-index: " << val << " node: " <<  getNode(val) << "\n";
        }
#endif

        // merge overlapping nodes in order from highest nodeid to lowest, so that the
        // merge process does not invalid nodeids of higher value
        //
        // to do this, we first need to find a node corresponding the input node, and sort
        // the remaining nodes by nodeIndex:
        //
        NodeAddressType inputSuperAddy;
        std::vector<NodeAddressType> nodeIndices;
        {
            bool isInputSuperFound(false);
            const known_pos_range2& inputRange(getLocus(startLocusIndex).getNode(nodeIndex).interval.range);

            BOOST_FOREACH(const NodeAddressType& val, intersectNodes)
            {
                assert(val.first==headLocusIndex);

                // one node must be a superset of the input node, find this and store separately:
                if ((! isInputSuperFound) && getNode(val).interval.range.is_superset_of(inputRange))
                {
                    inputSuperAddy=val;
                    isInputSuperFound=true;
                    continue;
                }
                nodeIndices.push_back(val);
            }
            assert(isInputSuperFound);
            std::sort(nodeIndices.rbegin(),nodeIndices.rend());
        }

        // merge this inputNode with each intersecting inputNode,
        // and eliminate the intersecting node:
        //
        NodeAddressType mergeTargetAddy(inputSuperAddy);
        BOOST_FOREACH(NodeAddressType nodeAddy, nodeIndices)
        {
            if (nodeAddy<mergeTargetAddy) std::swap(nodeAddy,mergeTargetAddy);
#ifdef DEBUG_SVL
            log_os << "MergeAndRemove: " << nodeAddy << "\n";
#endif
            mergeNodePtr(nodeAddy,mergeTargetAddy);
            removeNode(nodeAddy);
#ifdef DEBUG_SVL
            log_os << "Finished: " << nodeAddy << "\n";
            checkState();
#endif
        }
    }

    if (startLocusIndex != headLocusIndex)
    {
#ifdef DEBUG_SVL
        log_os << "clearLocusIndex: " << startLocusIndex << "\n";
#endif

        clearLocus(startLocusIndex);
    }

#ifdef DEBUG_SVL
    checkState(true,true);
#endif
}



void
SVLocusSet::
merge(const SVLocusSet& inputSet)
{
    // TODO: check for compatible bam headers between inputSet and this

    assert(getMinMergeEdgeCount() == inputSet.getMinMergeEdgeCount());

    BOOST_FOREACH(const SVLocus& locus, inputSet._loci)
    {
        try
        {
            merge(locus);
        }
        catch (...)
        {
            log_os << "ERROR: SVLocusSet merge failed.\n"
                   << "\tSVLocusSet source: " << inputSet.getSource() << "\n"
                   << "\tSVLocus index: " << locus.getIndex() << "\n";
            throw;
        }
    }
}



void
SVLocusSet::
getNodeIntersectCore(
    const LocusIndexType locusIndex,
    const NodeIndexType nodeIndex,
    const LocusSetIndexerType& searchNodes,
    const LocusIndexType filterLocusIndex,
    std::set<NodeAddressType>& intersectNodes) const
{
    typedef LocusSetIndexerType::const_iterator in_citer;

    intersectNodes.clear();

#ifdef DEBUG_SVL
    log_os << "SVLocusSet::getNodeIntersectCore inputNode: " << locusIndex << ":" << nodeIndex << " " << getNode(std::make_pair(locusIndex,nodeIndex));
    checkState();
#endif

    // get all existing nodes which intersect with this one:
    const NodeAddressType inputAddy(std::make_pair(locusIndex,nodeIndex));
    in_citer it(searchNodes.lower_bound(inputAddy));
    const GenomeInterval& inputInterval(getNode(inputAddy).interval);
    const pos_t maxRegionSize(_maxRegionSize[inputInterval.tid]);

    const in_citer it_begin(searchNodes.begin()), it_end(searchNodes.end());

    // first look forward and extend to find all nodes which this inputNode intersects:
    for (in_citer it_fwd(it); it_fwd != it_end; ++it_fwd)
    {
        if (it_fwd->first == filterLocusIndex) continue;
#ifdef DEBUG_SVL
        log_os << "FWD test: " << (*it_fwd) << " " << getNode(*it_fwd);
#endif
        if (! inputInterval.isIntersect(getNode(*it_fwd).interval)) break;
        intersectNodes.insert(*it_fwd);
#ifdef DEBUG_SVL
        log_os << "FWD insert: " << (*it_fwd) << "\n";
#endif
    }

    // now find all intersecting nodes in reverse direction:
    for (in_citer it_rev(it); it_rev != it_begin; )
    {
        --it_rev;
        if (it_rev->first == filterLocusIndex) continue;
#ifdef DEBUG_SVL
        log_os << "REV test: " << (*it_rev) << " " << getNode(*it_rev);
#endif
        const GenomeInterval& searchInterval(getNode(*it_rev).interval);
        if (! inputInterval.isIntersect(searchInterval))
        {
            if(! isOverlapAllowed()) break;

            if(inputInterval.tid != searchInterval.tid) break;
            if((searchInterval.range.begin_pos()+maxRegionSize)<inputInterval.range.begin_pos()) break;
            continue;
        }

        intersectNodes.insert(*it_rev);
#ifdef DEBUG_SVL
        log_os << "REV insert: " << (*it_rev) << "\n";
#endif
    }
}



void
SVLocusSet::
getNodeMergeableIntersect(
    const LocusIndexType locusIndex,
    const NodeIndexType nodeIndex,
    const bool isInputLocusMoved,
    std::set<NodeAddressType>& mergeIntersectNodes) const
{
    // Two ways sets of mergeable nodes can occur:
    // (1) There is a set of nodes which overlap with both input Node and one of the remote nodes that the input points to. When totaled together,
    // the edge count of this set + the inputNode edge exceeds minMergeEdgeCount.
    // (2) The input node either contains an edge which is greater than minMergeEdgeCOunt or will contain such an edge due to (1),
    // in this case the input node can be merged with an locally overlapping node which also contains an edge which is greater than minMergeEdgeCount
    //

    const NodeAddressType inputAddy(std::make_pair(locusIndex,nodeIndex));
    const SVLocusNode& inputNode(getNode(inputAddy));

#ifdef DEBUG_SVL
    log_os << "SVLocusSet::getNodeMergableIntersect inputNode: " << inputAddy << " " << inputNode;
    checkState();
#endif

    //
    // build a new index, which contains an enumeration of remote nodes for each intersecting node,
    // and a map pointing back to the intersecting locals for each remote:
    //
    typedef std::multimap<NodeAddressType, NodeIndexType> rlmap_t;
    typedef rlmap_t::const_iterator rliter_t;
    typedef std::pair<rliter_t,rliter_t> rlmap_range_t;

    rlmap_t remoteToLocal;
    LocusSetIndexerType remoteIntersect(*this);

    // these nodes intersect the input and already qualify as non-noise:
    std::set<NodeAddressType> signalIntersectNodes;
    {
        // get a standard intersection of the input node:
        std::set<NodeAddressType> intersectNodes;
        getNodeIntersect(locusIndex,nodeIndex,intersectNodes);
        BOOST_FOREACH(const NodeAddressType& addy, intersectNodes)
        {
            const SVLocus& intersectLocus(getLocus(addy.first));
            const SVLocusNode& intersectNode(getNode(addy));

            // for each remote node of the input, get all existing nodes which intersect with it:
            BOOST_FOREACH(const SVLocusNode::edges_type::value_type& intersectEdge, intersectNode)
            {
                // 1. build remote <-> local indexing structures:
                NodeAddressType remoteAddy(std::make_pair(addy.first,intersectEdge.first));
                remoteToLocal.insert(std::make_pair(remoteAddy,addy.second));
                remoteIntersect.insert(remoteAddy);
            }

            // 2. build the signal node set:
            if(! intersectLocus.isNoiseNode(_minMergeEdgeCount,addy.second))
            {
                signalIntersectNodes.insert(addy);
            }
        }

#ifdef DEBUG_SVL
        log_os << "SVLocusSet::getNodeMergableIntersect remoteIntersect.size(): " << remoteIntersect.size() << "\n";
        BOOST_FOREACH(const NodeAddressType& addy, remoteIntersect)
        {
            log_os << "\tremoteIsect node: " << addy << " " << getNode(addy);
        }
#endif
    }

    //
    // next build mergeIntersect by running through all edges of the input node
    //
    mergeIntersectNodes.clear();

    // for each remote node of the input, get all existing nodes which intersect with it:
    BOOST_FOREACH(const SVLocusNode::edges_type::value_type& inputEdge, inputNode)
    {
#ifdef DEBUG_SVL
        log_os << "SVLocusSet::getNodeMergableIntersect processing edge: " << inputEdge.first << "\n";
        checkState();
#endif

        // we need to total the edgecounts of all intersecting edges before determining if these can be added to mergeIntersect,
        // so we create an intermediate store here:
        std::set<NodeAddressType> countStore;
        getNodeIntersectCore(locusIndex,inputEdge.first,remoteIntersect,locusIndex,countStore);

        // total counts for this edge:
        unsigned mergedRemoteEdgeCount(0);
        unsigned mergedLocalEdgeCount(0);
        BOOST_FOREACH(const NodeAddressType remoteIsectAddy, countStore)
        {

#ifdef DEBUG_SVL
            log_os << "SVLocusSet::getNodeMergableIntersect countStore addy: " << remoteIsectAddy << "\n";
#endif

            const rlmap_range_t remoteIsectRange(remoteToLocal.equal_range(remoteIsectAddy));
            assert(remoteIsectRange.first != remoteToLocal.end());
            for(rliter_t riter(remoteIsectRange.first); riter != remoteIsectRange.second; ++riter)
            {
                const NodeIndexType localNodeIndex(riter->second);
                const SVLocus& remoteIsectLocus(getLocus(remoteIsectAddy.first));

                // total edge counts on the remote->local edge:
                mergedRemoteEdgeCount += remoteIsectLocus.getEdge(remoteIsectAddy.second,localNodeIndex).count;

                // total edge counts on the local->remote edge:
                mergedLocalEdgeCount += remoteIsectLocus.getEdge(localNodeIndex,remoteIsectAddy.second).count;
            }
        }

#ifdef DEBUG_SVL
        log_os << "SVLocusSet::getNodeMergableIntersect pre-input merge counts local: " << mergedLocalEdgeCount << " remote: " << mergedRemoteEdgeCount << "\n";
        checkState();
#endif

        if(! isInputLocusMoved)
        {
            { // total edge counts on the input remote->local edge
                const SVLocus& inputLocus(getLocus(inputAddy.first));
                mergedRemoteEdgeCount += inputLocus.getEdge(inputEdge.first,nodeIndex).count;
            }

            // total edge counts on the input local->remote edge
            mergedLocalEdgeCount += inputEdge.second.count;
        }


#ifdef DEBUG_SVL
        log_os << "SVLocusSet::getNodeMergableIntersect final merge counts local: " << mergedLocalEdgeCount << " remote: " << mergedRemoteEdgeCount << "\n";
        checkState();
#endif

        if((mergedLocalEdgeCount < _minMergeEdgeCount) &&
           (mergedRemoteEdgeCount < _minMergeEdgeCount)) continue;

        // Add type (1) mergeable nodes:
        std::set<NodeAddressType> thisEdgeMergeIntersectNodes;
        BOOST_FOREACH(const NodeAddressType remoteAddy, countStore)
        {
            const rlmap_range_t remoteIsectRange(remoteToLocal.equal_range(remoteAddy));
            assert(remoteIsectRange.first != remoteToLocal.end());
            for(rliter_t riter(remoteIsectRange.first); riter != remoteIsectRange.second; ++riter)
            {
                const NodeAddressType localIntersectAddy(std::make_pair(remoteAddy.first,riter->second));
                thisEdgeMergeIntersectNodes.insert(localIntersectAddy);
                mergeIntersectNodes.insert(localIntersectAddy);
            }
        }

        /// for each type (1) node, add any new intersections to the signal node set:
        ///
        /// this is not very efficient for now -- each (1) edge added in potentially expands the current node to intersect new signal nodes
        /// -- this loop looks for thos new signal nodes
        ///
        std::set<NodeAddressType> intersectNodes;
        BOOST_FOREACH(const NodeAddressType mergeAddy, thisEdgeMergeIntersectNodes)
        {
            // get a standard intersection of the input node:
            getNodeIntersectCore(mergeAddy.first,mergeAddy.second, _inodes,locusIndex,intersectNodes);
            BOOST_FOREACH(const NodeAddressType& intersectAddy, intersectNodes)
            {
                const SVLocus& intersectLocus(getLocus(intersectAddy.first));
                if(intersectLocus.isNoiseNode(_minMergeEdgeCount,intersectAddy.second)) continue;

#ifdef DEBUG_SVL
                if(signalIntersectNodes.count(intersectAddy) == 0)
                {
                    log_os << "SVLocusSet::getNodeMergableIntersect signal boost merge/new: " << mergeAddy << " " << intersectAddy << "\n";
                }
#endif

                signalIntersectNodes.insert(intersectAddy);
            }
        }

        // Add type (2) mergeable nodes:
        BOOST_FOREACH(const NodeAddressType signalAddy, signalIntersectNodes)
        {
            mergeIntersectNodes.insert(signalAddy);
        }
    }
}



void
SVLocusSet::
getRegionIntersect(
    const GenomeInterval interval,
    std::set<NodeAddressType>& intersectNodes)
{
    const LocusIndexType startLocusIndex(insertLocus(SVLocus()));
    const NodeIndexType nodeIndex = getLocus(startLocusIndex).addNode(interval);

    getNodeIntersect(startLocusIndex,nodeIndex,intersectNodes);

    clearLocus(startLocusIndex);
}



void
SVLocusSet::
moveIntersectToLowIndex(
    const std::set<NodeAddressType>& intersectNodes,
    const LocusIndexType startLocusIndex,
    LocusIndexType& locusIndex)
{
    const unsigned startHeadLocusIndex(locusIndex);

    // assign all intersect clusters to the lowest index number
    const bool isClearSource(startLocusIndex!=startHeadLocusIndex);

    // get lowest index number that is not startLocusIndex:
    bool isFirst(true);
    BOOST_FOREACH(const LocusSetIndexerType::value_type& val, intersectNodes)
    {
        if ((!isFirst) && (val.first >= locusIndex)) continue;
        locusIndex = val.first;
        isFirst=false;
    }

    combineLoci(startHeadLocusIndex,locusIndex,isClearSource);
    BOOST_FOREACH(const LocusSetIndexerType::value_type& val, intersectNodes)
    {
        combineLoci(val.first,locusIndex);
    }

#ifdef DEBUG_SVL
    log_os << "Reassigned all intersecting nodes to index: " << locusIndex << " shli: " << startHeadLocusIndex << " sli:" << startLocusIndex << "\n";
    checkState();
#endif
}



void
SVLocusSet::
combineLoci(
    const LocusIndexType fromIndex,
    const LocusIndexType toIndex,
    const bool isClearSource)
{
    assert(toIndex<_loci.size());

#ifdef DEBUG_SVL
    log_os << "combineLoci: from: " << fromIndex << " toIndex: " << toIndex << " isClear:" << isClearSource << "\n";
#endif

    if (fromIndex == toIndex) return;
    if (fromIndex>=_loci.size()) return;

    SVLocus& fromLocus(_loci[fromIndex]);
    if (fromLocus.empty()) return;

    SVLocus& toLocus(_loci[toIndex]);
    toLocus.copyLocus(fromLocus);
    if (isClearSource) clearLocus(fromIndex);
}



LocusIndexType
SVLocusSet::
insertLocus(
    const SVLocus& inputLocus)
{
    LocusIndexType locusIndex(0);
    if (_emptyLoci.empty())
    {
        static const unsigned maxIndex(std::numeric_limits<LocusIndexType>::max());
        locusIndex=_loci.size();
        assert(locusIndex<maxIndex);
        _loci.resize(locusIndex+1);
    }
    else
    {
        locusIndex=(*_emptyLoci.begin());
        assert(_loci[locusIndex].empty());
        _emptyLoci.erase(locusIndex);
    }

    SVLocus& locus(_loci[locusIndex]);
    observe_notifier(locus);
    locus.updateIndex(locusIndex);
    locus.copyLocus(inputLocus);
    return locusIndex;
}



void
SVLocusSet::
mergeNodePtr(NodeAddressType fromPtr,
             NodeAddressType toPtr)
{
#ifdef DEBUG_SVL
    log_os << "MergeNode: from: " << fromPtr << " to: " << toPtr << " fromLocusSize: " << getLocus(fromPtr.first).size() << "\n";
#endif
    LocusSetIndexerType::iterator iter(_inodes.find(toPtr));
    assert(iter != _inodes.end());
    assert(fromPtr.first == toPtr.first);
    getLocus(fromPtr.first).mergeNode(fromPtr.second,toPtr.second);
}



void
SVLocusSet::
clean()
{
    BOOST_FOREACH(SVLocus& locus, _loci)
    {
        if(locus.empty()) continue;
        _totalCleaned += locus.clean(getMinMergeEdgeCount());
        if(locus.empty()) _emptyLoci.insert(locus.getIndex());
    }
}



void
SVLocusSet::
cleanRegion(const GenomeInterval interval)
{
#ifdef DEBUG_SVL
    log_os << "cleanRegion interval: " << interval << "\n";
#endif

    std::set<NodeAddressType> intersectNodes;
    getRegionIntersect(interval,intersectNodes);

    BOOST_FOREACH(const NodeAddressType& val, intersectNodes)
    {
        SVLocus& locus(getLocus(val.first));
        if(locus.empty()) continue;
        _totalCleaned += locus.cleanNode(getMinMergeEdgeCount(), val.second);
        if(locus.empty()) _emptyLoci.insert(locus.getIndex());

#ifdef DEBUG_SVL
        log_os << "cleanRegion intersect: " << val << " empty_after_clean: " << locus.empty() << "\n";
#endif
    }
}



void
SVLocusSet::
dump(std::ostream& os) const
{
    os << "LOCUSSET_START\n";
    BOOST_FOREACH(const SVLocus& locus, _loci)
    {
        os << locus;
    }
    os << "LOCUSSET_END\n";
}



void
SVLocusSet::
dumpRegion(std::ostream& os,
           const GenomeInterval interval)
{
    std::set<NodeAddressType> intersectNodes;
    getRegionIntersect(interval,intersectNodes);

    LocusSetIndexerType sortedNodes(*this);
    BOOST_FOREACH(const NodeAddressType& val, intersectNodes)
    {
        sortedNodes.insert(val);
    }

    BOOST_FOREACH(const NodeAddressType& val, sortedNodes)
    {
        os << "SVNode LocusIndex:NodeIndex : " << val << "\n";
        os << getNode(val);
    }
}



void
SVLocusSet::
dumpStats(std::ostream& os) const
{
    static const char sep('\t');

    os << "locusIndex"
       << sep << "nodeCount"
       << sep << "nodeObsCount"
       << sep << "maxNodeObsCount"
       << sep << "regionSize"
       << sep << "maxRegionSize"
       << sep << "edgeCount"
       << sep << "maxEdgeCount"
       << sep << "edgeObsCount"
       << sep << "maxEdgeObsCount"
       << '\n';

    LocusIndexType locusIndex(0);
    BOOST_FOREACH(const SVLocus& locus, _loci)
    {
        unsigned locusNodeObsCount(0), maxNodeObsCount(0);
        unsigned locusRegionSize(0), maxRegionSize(0);
        unsigned locusEdgeCount(0), maxEdgeCount(0), locusEdgeObsCount(0), maxEdgeObsCount(0);
        BOOST_FOREACH(const SVLocusNode& node, locus)
        {
            // nodes:
            const unsigned nodeObsCount(node.count);
            maxNodeObsCount = std::max(maxNodeObsCount,nodeObsCount);
            locusNodeObsCount += nodeObsCount;

            // regions:
            const unsigned regionSize(node.interval.range.size());
            maxRegionSize = std::max(maxRegionSize,regionSize);
            locusRegionSize += regionSize;

            // edges:
            maxEdgeCount = std::max(maxEdgeCount,node.size());
            locusEdgeCount += node.size();
            BOOST_FOREACH(const SVLocusNode::edges_type::value_type& edge, node)
            {
                const unsigned edgeObsCount(edge.second.count);
                maxEdgeObsCount = std::max(maxEdgeObsCount,edgeObsCount);
                locusEdgeObsCount += edgeObsCount;
            }
        }
        os << locusIndex
           << sep << locus.size()
           << sep << locusNodeObsCount
           << sep << maxNodeObsCount
           << sep << locusRegionSize
           << sep << maxRegionSize
           << sep << locusEdgeCount
           << sep << maxEdgeCount
           << sep << locusEdgeObsCount
           << sep << maxEdgeObsCount
           << "\n";
        locusIndex++;
    }
}



void
SVLocusSet::
save(const char* filename) const
{
    using namespace boost::archive;

    assert(NULL != filename);
    std::ofstream ofs(filename, std::ios::binary);
    binary_oarchive oa(ofs);

    oa << header;
    oa << _minMergeEdgeCount;
    oa << _isFinalized;
    oa << _totalCleaned;
    BOOST_FOREACH(const SVLocus& locus, _loci)
    {
        if (locus.empty()) continue;
        oa << locus;
    }
}



void
SVLocusSet::
load(const char* filename)
{
    using namespace boost::archive;

    clear();

    assert(NULL != filename);
    std::ifstream ifs(filename, std::ios::binary);
    binary_iarchive ia(ifs);

    _source=filename;

    ia >> header;
    ia >> _minMergeEdgeCount;
    ia >> _isFinalized;
    ia >> _totalCleaned;
    SVLocus locus;
    while (ifs.peek() != EOF)
    {
        locus.clear();
        ia >> locus;
        if (locus.empty()) continue;
        const LocusIndexType locusIndex(size());
        _loci.push_back(locus);
        SVLocus& locusCopy(_loci.back());
        observe_notifier(locusCopy);
        locusCopy.updateIndex(locusIndex);
    }

    reconstructIndex();
    checkState(true,true);
}



void
SVLocusSet::
reconstructIndex()
{
    clearIndex();

    LocusIndexType locusIndex(0);
    BOOST_FOREACH(SVLocus& locus, _loci)
    {
        const unsigned nodeCount(locus.size());
        for (NodeIndexType nodeIndex(0); nodeIndex<nodeCount; ++nodeIndex)
        {
            const NodeAddressType addy(std::make_pair(locusIndex,nodeIndex));
            _inodes.insert(addy);
            updateMaxRegionSize(getNode(addy).interval);
        }
        if (locus.empty()) _emptyLoci.insert(locusIndex);
        locusIndex++;
    }
}



void
SVLocusSet::
dumpIndex(std::ostream& os) const
{
    os << "SVLocusSet Index START\n";
    BOOST_FOREACH(const NodeAddressType& in, _inodes)
    {
        os << "SVNodeIndex: " << in << "\n";
    }
    os << "SVLocusSet Index END\n";
}



void
SVLocusSet::
checkState(
        const bool isCheckOverlap,
        const bool isCheckLocusConnected) const
{
    using namespace illumina::common;

    unsigned locusIndex(0);
    unsigned totalNodeCount(0);
    BOOST_FOREACH(const SVLocus& locus, _loci)
    {
        locus.checkState(isCheckLocusConnected);

        const unsigned nodeCount(locus.size());
        totalNodeCount += nodeCount;

        for (NodeIndexType nodeIndex(0); nodeIndex<nodeCount; ++nodeIndex)
        {
            LocusSetIndexerType::const_iterator citer(_inodes.find(std::make_pair(locusIndex,nodeIndex)));
            if (citer == _inodes.end())
            {
                std::ostringstream oss;
                oss << "ERROR: locus node is missing from node index\n"
                    << "\tNode index: " << locusIndex << " node: " << getNode(std::make_pair(locusIndex,nodeIndex));
                BOOST_THROW_EXCEPTION(PreConditionException(oss.str()));
            }
            if ((citer->first != locusIndex) || (citer->second != nodeIndex))
            {
                std::ostringstream oss;
                oss << "ERROR: locus node is has conflicting index number in node index\n"
                    << "\tinode index_value: " << citer->first << ":" << citer->second << "\n"
                    << "\tNode index: " << locusIndex << ":" << locusIndex << " node: " << getNode(std::make_pair(locusIndex,nodeIndex));
                BOOST_THROW_EXCEPTION(PreConditionException(oss.str()));
            }
        }
        locusIndex++;
    }

    if (totalNodeCount != _inodes.size())
    {
        using namespace illumina::common;
        std::ostringstream oss;
        oss << "ERROR: SVLocusSet conflicting internal node counts. totalNodeCount: " << totalNodeCount << " inodeSize: " << _inodes.size() << "n";
        BOOST_THROW_EXCEPTION(PreConditionException(oss.str()));
    }

    if (! isCheckOverlap) return;

    if(isOverlapAllowed()) return;

    bool isFirst(true);
    GenomeInterval lastInterval;
    NodeAddressType lastAddy;
    BOOST_FOREACH(const NodeAddressType& addy, _inodes)
    {
        const GenomeInterval& interval(getNode(addy).interval);

        // don't allow zero-length or negative intervals:
        assert(interval.range.begin_pos() < interval.range.end_pos());

        // don't allow overlapping intervals:
        if (isFirst)
        {
            isFirst=false;
        }
        else if (interval.tid == lastInterval.tid)
        {
            if (lastInterval.range.end_pos() > interval.range.begin_pos())
            {
                std::ostringstream oss;
                oss << "ERROR: Overlapping nodes in graph\n"
                    << "\tlast_index: " << lastAddy << " interval: " << lastInterval << "\n"
                    << "\tthis_index: " << addy << " interval: " << interval << "\n"
                    << "\tlast_node: " << lastAddy << " "<< getNode(lastAddy) << "\n"
                    << "\tthis_node: " << addy << " "<< getNode(addy) << "\n";
                BOOST_THROW_EXCEPTION(PreConditionException(oss.str()));
            }
        }
        lastAddy = addy;
        lastInterval = interval;
    }
}

