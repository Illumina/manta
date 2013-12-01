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
/// \author Chris Saunders
///

#include "blt_util/log.hh"
#include "common/Exceptions.hh"
#include "svgraph/SVLocusSet.hh"

#include "boost/archive/binary_iarchive.hpp"
#include "boost/archive/binary_oarchive.hpp"
#include "boost/foreach.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>


void
SampleReadCounts::
write(
    std::ostream& os,
    const char* label) const
{
    static const char sep('\t');
    os << label << "_totalAnomalousConsidered:" << sep << anom << '\n';
    os << label << "_totalNonAnomalousConsidered:" << sep << nonAnom << '\n';
}



std::ostream&
operator<<(std::ostream& os, const SVLocusSet::NodeAddressType& a)
{
    os << a.first << ":" << a.second;
    return os;
}



void
SVLocusSet::
locusHurl(const LocusIndexType index, const char* label) const
{
    using namespace illumina::common;

    std::ostringstream oss;
    oss << "ERROR: Attempting to access locus: " << index << " in locusSet with size: " << size() << " accessLabel: " << label << "\n";
    BOOST_THROW_EXCEPTION(LogicException(oss.str()));
}



/// is the node set from multiple loci?
static
bool
isMultiLocus(
    const LocusIndexType locusIndex,
    const std::set<SVLocusSet::NodeAddressType>& nodes)
{
    BOOST_FOREACH(const SVLocusSet::NodeAddressType& addy, nodes)
    {
        if (addy.first != locusIndex) return true;
    }
    return false;
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

    // meaningless input indicates an error in client code:
    assert(! inputLocus.empty());

    if (inputLocus.empty()) return;

#ifdef DEBUG_SVL
    static const std::string logtag("SVLocusSet::merge");
    log_os << logtag << " inputLocus: " << inputLocus;
    checkState(true);
#endif

    inputLocus.checkState(true);

    const LocusIndexType startLocusIndex(insertLocus(inputLocus));
    const SVLocus& startLocus(_loci[startLocusIndex]);
    LocusIndexType headLocusIndex(startLocusIndex);

    // indicates if the input locus has been 'moved' into another locus in the graph:
    bool isInputLocusMoved(false);

    // indicates that the locus will not be inserted into the graph.
    // if true, skip merge and clear out the startLocus
    bool isAbortMerge(false);

    // because we have a non-general interval overlap test, we must order search
    // nodes by begin_pos on each chromosome
    //
    typedef std::map<GenomeInterval,NodeIndexType> nodeMap_t;
    nodeMap_t nodeMap;
    {
        const NodeIndexType nodeCount(startLocus.size());
        for (NodeIndexType nodeIndex(0); nodeIndex<nodeCount; ++nodeIndex)
        {
            nodeMap.insert(std::make_pair(startLocus.getNode(nodeIndex).getInterval(),nodeIndex));
        }
    }

    // reuse this intersectNodes object throughout the merge:
    std::set<NodeAddressType> intersectNodes;

    // test if the graph has grown too complex in these regions. If so, abort the insertion of this locus:
    BOOST_FOREACH(const nodeMap_t::value_type& nodeVal, nodeMap)
    {
        static const bool isTestUsability(true);

        // get a standard intersection of the input node:
        const bool isUsable(getNodeIntersect(startLocusIndex, nodeVal.second, intersectNodes, isTestUsability));

        if (! isUsable)
        {
#ifdef DEBUG_SVL
            log_os << logtag << "Aborting merge\n";
#endif
            isAbortMerge=true;
            break;
        }
    }

    BOOST_FOREACH(const nodeMap_t::value_type& nodeVal, nodeMap)
    {
        if (isAbortMerge) break;

        const NodeIndexType nodeIndex(nodeVal.second);

#ifdef DEBUG_SVL
        log_os << logtag << " inputNode: " << NodeAddressType(std::make_pair(startLocusIndex,nodeIndex)) << " " << startLocus.getNode(nodeIndex);
#endif

        getNodeMergeableIntersect(startLocusIndex, nodeIndex, isInputLocusMoved, intersectNodes);

#ifdef DEBUG_SVL
        log_os << logtag << " insersect_size: " << intersectNodes.size() << "\n";
        BOOST_FOREACH(const NodeAddressType& val, intersectNodes)
        {
            log_os << logtag << " intersect address: " << val << " node: " <<  getNode(val) << "\n";
        }
#endif

        if (isInputLocusMoved)
        {
            if (2>intersectNodes.size()) continue;
        }
        else
        {
            if (intersectNodes.empty()) continue;
        }

        while ( isMultiLocus(headLocusIndex, intersectNodes) )
        {
            // if there are any intersections, copy the loci of all intersecting nodes into
            // a single locus, by convention we use the lowest locusIndex of the intersecting set
            moveIntersectToLowIndex(intersectNodes,startLocusIndex,headLocusIndex);
            if (! isInputLocusMoved) isInputLocusMoved=(headLocusIndex != startLocusIndex);

            getNodeMergeableIntersect(startLocusIndex, nodeIndex, isInputLocusMoved, intersectNodes);
            assert(! intersectNodes.empty());
        }

#ifdef DEBUG_SVL
        log_os << logtag << " intersect2_size: " << intersectNodes.size() << "\n";
        BOOST_FOREACH(const NodeAddressType& val, intersectNodes)
        {
            log_os << logtag << " intersect2 address: " << val << " node: " <<  getNode(val) << "\n";
        }
#endif

        // merge overlapping nodes in order from highest nodeid to lowest, so that the
        // merge process does not invalidate nodeids of higher value
        //
        // we first need to find a node corresponding to the input node (but possibly merged to a larger region:
        //
        NodeAddressType inputSuperAddy;
        {
            bool isInputSuperFound(false);
            const known_pos_range2& inputRange(getLocus(startLocusIndex).getNode(nodeIndex).getInterval().range);

            BOOST_FOREACH(const NodeAddressType& val, intersectNodes)
            {
                assert(val.first == headLocusIndex);

                // one node must be a superset of the input node, find this and store separately:
                if (getNode(val).getInterval().range.is_superset_of(inputRange))
                {
                    inputSuperAddy=val;
                    isInputSuperFound=true;
                    break;
                }
            }
            assert(isInputSuperFound);
        }

        // merge this inputNode with each intersecting Node,
        // and eliminate the intersecting node:
        //
        NodeAddressType mergeTargetAddy(inputSuperAddy);
        BOOST_REVERSE_FOREACH(NodeAddressType nodeAddy, intersectNodes)
        {
            if (nodeAddy == inputSuperAddy) continue;
            if (nodeAddy < mergeTargetAddy) std::swap(nodeAddy,mergeTargetAddy);
#ifdef DEBUG_SVL
            log_os << logtag << " MergeAndRemove: " << nodeAddy << "\n";
#endif
            mergeNodePtr(nodeAddy,mergeTargetAddy);
            removeNode(nodeAddy);
#ifdef DEBUG_SVL
            log_os << logtag << " Finished: " << nodeAddy << "\n";
            checkState();
#endif
        }
    }

    if (isAbortMerge || isInputLocusMoved)
    {
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

    _totalCleaned += inputSet._totalCleaned;
    _normalReads.merge(inputSet._normalReads);
    _tumorReads.merge(inputSet._tumorReads);
    _highestSearchCount = std::max(_highestSearchCount, inputSet._highestSearchCount);
    _isMaxSearchCount = (_isMaxSearchCount || inputSet._isMaxSearchCount);
    _highestSearchDensity = std::max(_highestSearchDensity, inputSet._highestSearchDensity);
    _isMaxSearchDensity = (_isMaxSearchDensity || inputSet._isMaxSearchDensity);

}



bool
SVLocusSet::
getNodeIntersectCore(
    const LocusIndexType inputLocusIndex,
    const NodeIndexType inputNodeIndex,
    const LocusSetIndexerType& searchNodes,
    const LocusIndexType filterLocusIndex,
    std::set<NodeAddressType>& intersectNodes,
    const bool isTestUsability) const
{
    typedef LocusSetIndexerType::const_iterator in_citer;

#ifdef DEBUG_SVL
    static const std::string logtag("SVLocusSet::getNodeIntersectCore");
    log_os << logtag << " inputNode: " << inputLocusIndex << ":" << inputNodeIndex << " " << getNode(std::make_pair(inputLocusIndex,inputNodeIndex));
    checkState();
#endif

    assert(_isIndexed);

    intersectNodes.clear();

    // get all existing nodes which intersect with this one:
    const NodeAddressType inputAddy(std::make_pair(inputLocusIndex,inputNodeIndex));
    const in_citer it(searchNodes.lower_bound(inputAddy));
    const GenomeInterval& inputInterval(getNode(inputAddy).getInterval());
    const pos_t maxRegionSize(_maxRegionSize[inputInterval.tid]);

    const in_citer it_begin(searchNodes.begin()), it_end(searchNodes.end());

    // diagnostics to determine if graph is growing too dense in one region:
    bool isUsable(true);
    unsigned searchCount(0);

    // first look forward and extend to find all nodes which this inputNode intersects:
    for (in_citer it_fwd(it); it_fwd != it_end; ++it_fwd)
    {
        if (isTestUsability)
        {
            searchCount++;
            if (searchCount > _opt.maxSearchCount)
            {
                isUsable = false;
                _isMaxSearchCount=true;
                break;
            }
        }

        if (it_fwd->first == filterLocusIndex) continue;
#ifdef DEBUG_SVL
        log_os << logtag << "\tFWD test: " << (*it_fwd) << " " << getNode(*it_fwd);
#endif
        if (! inputInterval.isIntersect(getNode(*it_fwd).getInterval())) break;
        intersectNodes.insert(*it_fwd);
#ifdef DEBUG_SVL
        log_os << logtag << "\tFWD insert: " << (*it_fwd) << "\n";
#endif
    }

    // now find all intersecting nodes in reverse direction:
    for (in_citer it_rev(it); it_rev != it_begin; )
    {
        --it_rev;

        if (isTestUsability)
        {
            if (! isUsable) break;
            searchCount++;
            if (searchCount > _opt.maxSearchCount)
            {
                isUsable = false;
                _isMaxSearchCount=true;
                break;
            }
        }

        if (it_rev->first == filterLocusIndex) continue;
#ifdef DEBUG_SVL
        log_os << logtag << "\tREV test: " << (*it_rev) << " " << getNode(*it_rev);
#endif
        const GenomeInterval& searchInterval(getNode(*it_rev).getInterval());
        if (! inputInterval.isIntersect(searchInterval))
        {
            if (! isOverlapAllowed()) break;

            if (inputInterval.tid != searchInterval.tid) break;
            if ((searchInterval.range.begin_pos()+maxRegionSize)<inputInterval.range.begin_pos()) break;
            continue;
        }

        intersectNodes.insert(*it_rev);
#ifdef DEBUG_SVL
        log_os << logtag << "\tREV insert: " << (*it_rev) << "\n";
#endif
    }

    if (! isTestUsability) return true;

    _highestSearchCount = std::max(_highestSearchCount, searchCount);

    pos_t searchSize(inputInterval.range.end_pos() - std::max(0, inputInterval.range.begin_pos()-maxRegionSize));

    assert(searchSize>=0);
    if (0 != searchSize)
    {
        static const pos_t minSearchSize(40);
        searchSize = std::max(searchSize, minSearchSize);
        const float searchDensity(static_cast<float>(searchCount)/static_cast<float>(searchSize));
        _highestSearchDensity = std::max(_highestSearchDensity, searchDensity);

        if (searchDensity > _opt.maxSearchDensity)
        {
            isUsable = false;
            _isMaxSearchDensity=true;
        }
    }

    return isUsable;
}



void
SVLocusSet::
getIntersectingEdgeNodes(
    const LocusIndexType inputLocusIndex,
    const NodeIndexType inputRemoteNodeIndex,
    const EdgeMapType& remoteToLocal,
    const LocusSetIndexerType& remoteIntersect,
    std::vector<EdgeInfoType>& edges) const
{
    typedef EdgeMapType::const_iterator rliter_t;
    typedef std::pair<rliter_t,rliter_t> rlmap_range_t;

    edges.clear();

    // find all remote nodes that are part of edges which intersect the input:
    std::set<NodeAddressType> edgeIntersectRemoteTemp;
    getNodeIntersectCore(inputLocusIndex,inputRemoteNodeIndex,remoteIntersect,inputLocusIndex,edgeIntersectRemoteTemp);

    BOOST_FOREACH(const NodeAddressType remoteIsectAddy, edgeIntersectRemoteTemp)
    {
        // find what local nodes the remote notes trace back to:
        const rlmap_range_t remoteIsectRange(remoteToLocal.equal_range(remoteIsectAddy));
        assert(remoteIsectRange.first != remoteToLocal.end());
        for (rliter_t riter(remoteIsectRange.first); riter != remoteIsectRange.second; ++riter)
        {
            const NodeAddressType localIntersectAddy(std::make_pair(remoteIsectAddy.first,riter->second));
            edges.push_back(std::make_pair(localIntersectAddy,remoteIsectAddy.second));
        }
    }
}



void
SVLocusSet::
findSignalNodes(
    const LocusIndexType inputLocusIndex,
    const NodeAddressType findSignalAddy,
    std::set<NodeAddressType>& signalIntersectNodes,
    const std::set<NodeAddressType>& inputIntersectRemotes,
    bool& isIntersectRemotes) const
{
#ifdef DEBUG_SVL
    static const std::string logtag("SVLocusSet::findSignalNodes");
    log_os << logtag << "\tsignal_boost: findSignalAddy: " << findSignalAddy << "\n";
#endif
    // get a standard intersection of the input node:
    std::set<NodeAddressType> intersectNodes;
    getNodeIntersectCore(findSignalAddy.first, findSignalAddy.second, _inodes, inputLocusIndex, intersectNodes);
    BOOST_FOREACH(const NodeAddressType intersectAddy, intersectNodes)
    {
#ifdef DEBUG_SVL
        log_os << logtag << "\tsignal_boost: intersectAddy: " << intersectAddy << "\n";
#endif
        if (isNoiseNode(intersectAddy))
        {
            // check for the rare remote intersect condition:
            if (! isIntersectRemotes)
            {
                if (inputIntersectRemotes.count(intersectAddy))
                {
                    isIntersectRemotes=true;
                }
            }
            continue;
        }

#ifdef DEBUG_SVL
        if (signalIntersectNodes.count(intersectAddy) == 0)
        {
            log_os << logtag << " signal boost merge/new: " << findSignalAddy << " " << intersectAddy << "\n";
        }
#endif

        signalIntersectNodes.insert(intersectAddy);
    }
}



void
SVLocusSet::
getNodeMergeableIntersect(
    const LocusIndexType inputLocusIndex,
    const NodeIndexType inputNodeIndex,
    const bool isInputLocusMoved,
    std::set<NodeAddressType>& mergeIntersectNodes) const
{
    //
    // TODO: There's room for significant optimization of these methods. The improvements are not trivial,
    //   but they would allow us to filter fewer nodes from being merged when node intersection counts become large.
    //

    //
    // There are two ways sets of mergeable nodes can occur:
    //
    // (1) There is a set of nodes which overlap with both input node and one of the remote nodes that the input points to (ie they have a shared edge).
    // When totaled together, the edge count of this set + the inputNode edge exceeds minMergeEdgeCount.
    //
    // (2) The input node either contains an edge which is greater than minMergeEdgeCount or will contain such an edge due to (1),
    // in this case the input node can be merged with a locally overlapping node which also contains an edge which is greater than
    // minMergeEdgeCount. Note that in case (2) remote node intersection is not required.
    //

    const NodeAddressType inputAddy(std::make_pair(inputLocusIndex,inputNodeIndex));
    const SVLocusNode& inputNode(getNode(inputAddy));

#ifdef DEBUG_SVL
    static const std::string logtag("SVLocusSet::getNodeMergableIntersect");
    log_os << logtag << " inputNode: " << inputAddy << " " << inputNode;
    checkState();
#endif

    // reuse this intersectNodes as a temporary throughout the methods below
    std::set<NodeAddressType> intersectNodes;

    //
    // build a new index, which contains an enumeration of remote nodes for each intersecting node,
    // and a map pointing back to the intersecting locals for each remote:
    //
    EdgeMapType remoteToLocal;
    LocusSetIndexerType remoteIntersect(*this);

    // these nodes intersect the input and already qualify as non-noise:
    std::set<NodeAddressType> signalIntersectNodes;
    {
        // get a standard intersection of the input node:
        getNodeIntersect(inputLocusIndex, inputNodeIndex, intersectNodes);

        BOOST_FOREACH(const NodeAddressType& intersectAddy, intersectNodes)
        {
            const SVLocusNode& intersectNode(getNode(intersectAddy));

            // get the remotes of each node which intersect with the query node,
            // place these in remoteIntersect
            const SVLocusEdgeManager edgeMap(intersectNode.getEdgeManager());
            BOOST_FOREACH(const SVLocusEdgesType::value_type& intersectEdge, edgeMap.getMap())
            {
                // 1. build remote <-> local indexing structures:
                NodeAddressType remoteAddy(std::make_pair(intersectAddy.first,intersectEdge.first));
                remoteToLocal.insert(std::make_pair(remoteAddy,intersectAddy.second));
                remoteIntersect.insert(remoteAddy);
            }

            // 2. build the signal node set:
            if (isNoiseNode(intersectAddy)) continue;

            signalIntersectNodes.insert(intersectAddy);
        }

#ifdef DEBUG_SVL
        log_os << logtag << " remoteIntersect.size(): " << remoteIntersect.size() << "\n";
        BOOST_FOREACH(const NodeAddressType& addy, remoteIntersect)
        {
            log_os << logtag << "\tremoteIsect node: " << addy << " " << getNode(addy);
        }

        log_os << logtag << " signalIntersect.size(): " << signalIntersectNodes.size() << "\n";
        BOOST_FOREACH(const NodeAddressType& addy, signalIntersectNodes)
        {
            log_os << logtag << "\tsignalIsect node: " << addy << " " << getNode(addy);
        }
#endif
    }

    //
    // next build mergeIntersect by running through all edges of the input node
    //
    mergeIntersectNodes.clear();

    // for each edge from the input node, get all intersecting edges
    const SVLocusEdgeManager edgeMap(inputNode.getEdgeManager());
    BOOST_FOREACH(const SVLocusEdgesType::value_type& inputEdge, edgeMap.getMap())
    {
#ifdef DEBUG_SVL
        log_os << logtag << " processing edge: " << inputAddy << "->" << inputLocusIndex << ":" << inputEdge.first << "\n";
        checkState();
#endif

        // get the set of nodes that intersect the input node *and* have remotes which intersect the input remotes:
        std::vector<EdgeInfoType> inputIntersectEdges;
        getIntersectingEdgeNodes(inputLocusIndex, inputEdge.first, remoteToLocal, remoteIntersect, inputIntersectEdges);

        // total counts for this edge:
        unsigned mergedRemoteEdgeCount(0);
        unsigned mergedLocalEdgeCount(0);
        BOOST_FOREACH(const EdgeInfoType edgeInfo, inputIntersectEdges)
        {
            const SVLocus& edgeLocus(getLocus(edgeInfo.first.first));
            const NodeIndexType localNodeIndex(edgeInfo.first.second);
            const NodeIndexType remoteNodeIndex(edgeInfo.second);

            // total edge counts on the remote->local edge:
            mergedRemoteEdgeCount += edgeLocus.getEdge(remoteNodeIndex,localNodeIndex).getCount();

            // total edge counts on the local->remote edge:
            mergedLocalEdgeCount += edgeLocus.getEdge(localNodeIndex,remoteNodeIndex).getCount();
        }

#ifdef DEBUG_SVL
        log_os << logtag << " pre-input merge counts"
               << " local: " << mergedLocalEdgeCount
               << " remote: " << mergedRemoteEdgeCount
               << "\n";
        checkState();
#endif

        if (! isInputLocusMoved)
        {
            // if the input hasn't been moved into the primary locus graph yet, then we need to include the inputLocus
            // in order to get an accurate edge intersection count:
            {
                // total edge counts on the input remote->local edge
                const SVLocus& inputLocus(getLocus(inputAddy.first));
                mergedRemoteEdgeCount += inputLocus.getEdge(inputEdge.first,inputNodeIndex).getCount();
            }

            // total edge counts on the input local->remote edge
            mergedLocalEdgeCount += inputEdge.second.getCount();
        }


#ifdef DEBUG_SVL
        log_os << logtag << " final merge counts"
               << " local: " << mergedLocalEdgeCount
               << " remote: " << mergedRemoteEdgeCount
               << "\n";
        checkState();
#endif

        if ((mergedLocalEdgeCount < getMinMergeEdgeCount()) &&
            (mergedRemoteEdgeCount < getMinMergeEdgeCount())) continue;

        //
        // Add type1 mergeable nodes:
        //
        BOOST_FOREACH(const EdgeInfoType edgeInfo, inputIntersectEdges)
        {
            mergeIntersectNodes.insert(edgeInfo.first);
        }

        /// for each type1 node, add any new intersections to the signal node set:
        ///
        /// this is not very efficient for now -- each type1 edge added in potentially expands the current node to intersect new signal nodes
        /// -- this loop looks for those new signal nodes
        ///

        {
            /// this is used to search for the (rare) case where the intersection set locals overlap with the intersection set remotes
            std::set<NodeAddressType> inputIntersectRemotes;
            BOOST_FOREACH(const EdgeInfoType edgeInfo, inputIntersectEdges)
            {
                inputIntersectRemotes.insert(std::make_pair(edgeInfo.first.first,edgeInfo.second));
            }

            bool isIntersectRemotes(false);

            // check both the original node and intersected nodes for intersection to any of the group's remotes, and for new type2 signal intersect:
            findSignalNodes(inputLocusIndex, inputAddy, signalIntersectNodes, inputIntersectRemotes, isIntersectRemotes);
            BOOST_FOREACH(const EdgeInfoType edgeInfo, inputIntersectEdges)
            {
                findSignalNodes(inputLocusIndex, edgeInfo.first, signalIntersectNodes, inputIntersectRemotes, isIntersectRemotes);
            }

            if (isIntersectRemotes)
            {
                BOOST_FOREACH(const NodeAddressType intersectAddy, inputIntersectRemotes)
                {
#ifdef DEBUG_SVL
                    log_os << logtag << " adding ownRemote: " << intersectAddy << "\n";
#endif
                    mergeIntersectNodes.insert(intersectAddy);

                    // check to see if this adds even more signal nodes!
                    findSignalNodes(inputLocusIndex, intersectAddy, signalIntersectNodes, inputIntersectRemotes, isIntersectRemotes);
                }
            }
        }
        //
        // Add type2 mergeable nodes:
        //
        BOOST_FOREACH(const NodeAddressType signalAddy, signalIntersectNodes)
        {
            mergeIntersectNodes.insert(signalAddy);
        }
    }

#ifdef DEBUG_SVL
    log_os << logtag << " END. IntersectNodeSize: " << mergeIntersectNodes.size() << " Nodes:\n";
    BOOST_FOREACH(const NodeAddressType addy, mergeIntersectNodes)
    {
        log_os << logtag << "\tInode: " << addy << "\n";
    }
#endif
}



void
SVLocusSet::
getRegionIntersect(
    const GenomeInterval interval,
    std::set<NodeAddressType>& intersectNodes)
{
    const LocusIndexType startLocusIndex(insertLocus(SVLocus()));
    const NodeIndexType nodeIndex(getLocus(startLocusIndex).addNode(interval, this));

    getNodeIntersect(startLocusIndex, nodeIndex, intersectNodes);

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
    static const std::string logtag("SVLocusSet::moveIntersectToLowIndex");
    log_os << logtag << " Reassigned all intersecting nodes to locusIndex: " << locusIndex << " startHeadLocusIndex: " << startHeadLocusIndex << " startLocusIndex:" << startLocusIndex << "\n";
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
    static const std::string logtag("SVLocusSet::combineLoci");
    log_os << logtag << " from: " << fromIndex << " toIndex: " << toIndex << " isClear:" << isClearSource << "\n";
#endif

    if (fromIndex == toIndex) return;
    if (fromIndex>=_loci.size()) return;

    SVLocus& fromLocus(_loci[fromIndex]);
    if (fromLocus.empty()) return;

    SVLocus& toLocus(_loci[toIndex]);
    toLocus.copyLocus(fromLocus, this);
    if (isClearSource) clearLocus(fromIndex);
}



LocusIndexType
SVLocusSet::
insertLocus(
    const SVLocus& inputLocus)
{
    assert(_isIndexed);

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
    locus.updateIndex(locusIndex);
    locus.copyLocus(inputLocus, this);
    return locusIndex;
}



void
SVLocusSet::
mergeNodePtr(NodeAddressType fromPtr,
             NodeAddressType toPtr)
{
#ifdef DEBUG_SVL
    static const std::string logtag("SVLocusSet::mergeNodePtr");
    log_os << logtag << " from: " << fromPtr << " to: " << toPtr << " fromLocusSize: " << getLocus(fromPtr.first).size() << "\n";
#endif
    assert(_isIndexed);

    LocusSetIndexerType::iterator iter(_inodes.find(toPtr));
    assert(iter != _inodes.end());
    assert(fromPtr.first == toPtr.first);
    getLocus(fromPtr.first).mergeNode(fromPtr.second, toPtr.second, this);
}



void
SVLocusSet::
clean()
{
    BOOST_FOREACH(SVLocus& locus, _loci)
    {
        if (locus.empty()) continue;
        _totalCleaned += locus.clean(getMinMergeEdgeCount(), this);
        if (locus.empty()) _emptyLoci.insert(locus.getIndex());
    }
#ifdef DEBUG_SVL
    checkForOverlapNodes(true);
#endif
}



void
SVLocusSet::
cleanRegion(const GenomeInterval interval)
{
#ifdef DEBUG_SVL
    static const std::string logtag("SVLocusSet::cleanRegion");
    log_os << logtag << " interval: " << interval << "\n";
#endif

    std::set<NodeAddressType> intersectNodes;
    getRegionIntersect(interval, intersectNodes);

    // process nodes in reverse to properly handle instances when a locus has
    // multiple intersect nodes. This way we won't try to iterate into an
    // address which has been shifted by node deletion:
    BOOST_REVERSE_FOREACH(const NodeAddressType& val, intersectNodes)
    {
        SVLocus& locus(getLocus(val.first));
        if (locus.empty()) continue;
        _totalCleaned += locus.cleanNode(getMinMergeEdgeCount(), val.second, this);
        if (locus.empty()) _emptyLoci.insert(locus.getIndex());

#ifdef DEBUG_SVL
        log_os << logtag << " intersect: " << val << " is_empty_after_clean: " << locus.empty() << "\n";
#endif
    }
#ifdef DEBUG_SVL
    checkForOverlapNodes(true);
#endif
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

    os << "disjointSubgraphs:" << sep << nonEmptySize() << "\n";
    os << "nodes:" << sep << totalNodeCount() << "\n";
    os << "directedEdges:" << sep << totalEdgeCount() << "\n";
    os << "totalGraphEvidence:" << sep << totalObservationCount() << "\n";
    os << "totalCleaned:" << sep << _totalCleaned << "\n";
    os << "highestSearchCount:" << sep << _highestSearchCount << "\n";
    os << "isMaxSearchCount:" << sep << _isMaxSearchCount << "\n";
    os << "highestSearchDensity:" << sep << _highestSearchDensity << "\n";
    os << "isMaxSearchDensity:" << sep << _isMaxSearchDensity << "\n";

    _normalReads.write(os,"normal");
    _tumorReads.write(os,"tumor");
}


void
SVLocusSet::
dumpLocusStats(std::ostream& os) const
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
            const unsigned nodeObsCount(node.outCount());
            maxNodeObsCount = std::max(maxNodeObsCount,nodeObsCount);
            locusNodeObsCount += nodeObsCount;

            // regions:
            const unsigned regionSize(node.getInterval().range.size());
            maxRegionSize = std::max(maxRegionSize,regionSize);
            locusRegionSize += regionSize;

            // edges:
            maxEdgeCount = std::max(maxEdgeCount,node.size());
            locusEdgeCount += node.size();
            const SVLocusEdgeManager edgeMap(node.getEdgeManager());
            BOOST_FOREACH(const SVLocusEdgesType::value_type& edge, edgeMap.getMap())
            {
                const unsigned edgeObsCount(edge.second.getCount());
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
    oa << _opt;
    oa << _isFinalized;
    oa << _totalCleaned;
    oa << _normalReads;
    oa << _tumorReads;
    oa << _highestSearchCount;
    oa << _highestSearchDensity;
    oa << _isMaxSearchCount;
    oa << _isMaxSearchDensity;

    BOOST_FOREACH(const SVLocus& locus, _loci)
    {
        if (locus.empty()) continue;
        oa << locus;
    }
}



void
SVLocusSet::
load(
    const char* filename,
    const bool isSkipIndex)
{
    using namespace boost::archive;

#ifdef DEBUG_SVL
    log_os << "SVLocusSet::load BEGIN\n";
#endif

    clear();

    assert(NULL != filename);
    std::ifstream ifs(filename, std::ios::binary);
    binary_iarchive ia(ifs);

    _source=filename;

    ia >> header;
    ia >> _opt;
    ia >> _isFinalized;
    ia >> _totalCleaned;
    ia >> _normalReads;
    ia >> _tumorReads;
    ia >> _highestSearchCount;
    ia >> _highestSearchDensity;
    ia >> _isMaxSearchCount;
    ia >> _isMaxSearchDensity;

    SVLocus locus;
    while (ifs.peek() != EOF)
    {
        locus.clear(this);
        ia >> locus;
        if (locus.empty()) continue;
        const LocusIndexType locusIndex(size());
        _loci.push_back(locus);
        SVLocus& locusCopy(_loci.back());
        locusCopy.updateIndex(locusIndex);
    }

    if (! isSkipIndex)
    {
        reconstructIndex();
        checkState(true,true);
    }
    else
    {
        _isIndexed = false;
    }

#ifdef DEBUG_SVL
    log_os << "SVLocusSet::load END\n";
#endif
}



void
SVLocusSet::
reconstructIndex()
{
#ifdef DEBUG_SVL
    log_os << "reconstructIndex BEGIN\n";
#endif
    clearIndex();

#ifdef DEBUG_SVL
    log_os << "reconstructIndex cleared\n";
#endif

    LocusIndexType locusIndex(0);
    BOOST_FOREACH(SVLocus& locus, _loci)
    {
        const unsigned nodeCount(locus.size());
        for (NodeIndexType nodeIndex(0); nodeIndex<nodeCount; ++nodeIndex)
        {
            const NodeAddressType addy(std::make_pair(locusIndex,nodeIndex));
            _inodes.insert(addy);
            updateMaxRegionSize(getNode(addy).getInterval());
        }
        if (locus.empty()) _emptyLoci.insert(locusIndex);
        locusIndex++;
    }

    _isIndexed=true;

#ifdef DEBUG_SVL
    log_os << "reconstructIndex END\n";
#endif
}



void
SVLocusSet::
dumpIndex(std::ostream& os) const
{
    assert(_isIndexed);

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

    assert(_isIndexed);

    unsigned locusIndex(0);
    unsigned checkStateTotalNodeCount(0);
    BOOST_FOREACH(const SVLocus& locus, _loci)
    {
        locus.checkState(isCheckLocusConnected);

        const unsigned nodeCount(locus.size());
        checkStateTotalNodeCount += nodeCount;

        if (nodeCount == 0)
        {
            if (_emptyLoci.count(locusIndex) == 0)
            {
                std::ostringstream oss;
                oss << "ERROR: empty locus is not updated in the empty index. Locus index: " << locusIndex << "\n";
                BOOST_THROW_EXCEPTION(LogicException(oss.str()));
            }
        }

        for (NodeIndexType nodeIndex(0); nodeIndex<nodeCount; ++nodeIndex)
        {
            LocusSetIndexerType::const_iterator citer(_inodes.find(std::make_pair(locusIndex,nodeIndex)));
            if (citer == _inodes.end())
            {
                std::ostringstream oss;
                oss << "ERROR: locus node is missing from node index\n"
                    << "\tNode index: " << locusIndex << " node: " << getNode(std::make_pair(locusIndex,nodeIndex));
                BOOST_THROW_EXCEPTION(LogicException(oss.str()));
            }
            if ((citer->first != locusIndex) || (citer->second != nodeIndex))
            {
                std::ostringstream oss;
                oss << "ERROR: locus node has conflicting index number in node index\n"
                    << "\tinode index_value: " << citer->first << ":" << citer->second << "\n"
                    << "\tNode index: " << locusIndex << ":" << locusIndex << " node: " << getNode(std::make_pair(locusIndex,nodeIndex));
                BOOST_THROW_EXCEPTION(LogicException(oss.str()));
            }
        }
        locusIndex++;
    }

    if (checkStateTotalNodeCount != _inodes.size())
    {
        using namespace illumina::common;
        std::ostringstream oss;
        oss << "ERROR: SVLocusSet conflicting internal node counts. TotalNodeCount: " << checkStateTotalNodeCount << " inodeSize: " << _inodes.size() << "n";
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }

    if (! isCheckOverlap) return;

    // if isOverlapAllowed() then we should expect noise nodes to overlap, but we can still check signal nodes:
    const bool isFilterNoise(isOverlapAllowed());
    checkForOverlapNodes(isFilterNoise);
}



void
SVLocusSet::
checkForOverlapNodes(
    const bool isFilterNoise) const
{
    using namespace illumina::common;

    bool isFirst(true);
    GenomeInterval lastInterval;
    NodeAddressType lastAddy;
    BOOST_FOREACH(const NodeAddressType& addy, _inodes)
    {
        if (isFilterNoise)
        {
            if (isNoiseNode(addy)) continue;
        }

        const GenomeInterval& interval(getNode(addy).getInterval());

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
                    << "\tthis_node: " << addy << " "<< getNode(addy) << "\n"
                    << "\n"
                    << header << "\n";
                BOOST_THROW_EXCEPTION(LogicException(oss.str()));
            }
        }
        lastAddy = addy;
        lastInterval = interval;
    }
}

