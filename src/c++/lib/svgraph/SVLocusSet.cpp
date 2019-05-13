//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2019 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

/// \file
/// \author Chris Saunders
///

#include "svgraph/SVLocusSet.hpp"
#include "blt_util/SizeDistribution.hpp"
#include "blt_util/log.hpp"
#include "common/Exceptions.hpp"

#include "blt_util/thirdparty_push.h"

#include "boost/archive/binary_iarchive.hpp"
#include "boost/archive/binary_oarchive.hpp"
#include "boost/foreach.hpp"

#include "blt_util/thirdparty_pop.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

std::ostream& operator<<(std::ostream& os, const SVLocusSet::NodeAddressType& a)
{
  os << a.first << ":" << a.second;
  return os;
}

#ifdef DEBUG_SVL
void SVLocusSet::locusHurl(const LocusIndexType index, const char* label) const
{
  using namespace illumina::common;

  std::ostringstream oss;
  oss << "Attempting to access locus: " << index << " in locusSet with size: " << size()
      << " accessLabel: " << label;
  BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
}
#endif

/// Test if the set of node addresses come from only one, expected locus index.
///
/// \return True if all nodes addressed in \p nodeAddresses are in the expected locus.
static bool isExpectedLocusOnly(
    const LocusIndexType expectedLocusIndex, const std::set<SVLocusSet::NodeAddressType>& nodeAddresses)
{
  for (const SVLocusSet::NodeAddressType& nodeAddress : nodeAddresses) {
    if (nodeAddress.first != expectedLocusIndex) return false;
  }
  return true;
}

SVLocusSet::SVLocusSet(
    const SVLocusSetOptions&        opt,
    const bam_header_info&          bamHeaderInfo,
    const std::vector<std::string>& alignmentFilenames)
  : _bamHeaderInfo(bamHeaderInfo),
    _opt(opt),
    _inodes(*this),
    _source("UNKNOWN"),
    _isFinalized(false),
    _totalCleaned(0),
    _highestSearchCount(0),
    _highestSearchDensity(0),
    _isMaxSearchCount(false),
    _isMaxSearchDensity(false),
    _isIndexed(true)
{
  /// Initialize read counts:
  const unsigned sampleCount(alignmentFilenames.size());
  _counts.setSampleCount(sampleCount);
  for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex) {
    _counts.getSampleCounts(sampleIndex).sampleSource = alignmentFilenames[sampleIndex];
  }
}

void SVLocusSet::merge(const SVLocus& inputLocus)
{
  using namespace illumina::common;

  //
  // 1. Sanity check the inputLocus and this SVLocusSet to ensure input conditions for merging are met.
  //
  assert(!_isFinalized);

  // An empty inputLocus is meaningless input, and thus indicates an error in client code.
  assert(!inputLocus.empty());

  // In case we ever decide that an empty inputLocus is not a failed assertion, there is nothing to merge so
  // return.
  if (inputLocus.empty()) return;

#ifdef DEBUG_SVL
  static const std::string logtag("SVLocusSet::merge");
  log_os << logtag << " inputLocus: " << inputLocus;
  checkState(true);
#endif

  inputLocus.checkState(true);

  //
  // 2. Add the inputLocus as a new locus in SVLocusSet, called startLocus. In this step no merging takes
  // place.
  //
  const LocusIndexType startLocusIndex(insertLocus(inputLocus));
  const SVLocus&       startLocus(_loci[startLocusIndex]);
  LocusIndexType       headLocusIndex(startLocusIndex);

  // True if the startLocus has been copied from its original position at startLocusIndex into another locus
  // in the graph
  bool isStartLocusDuplicatedInAnotherGraphLocus(false);

  // True if this locus will not be merged into the graph because such a merge would exceed complexity limits.
  bool isAbortMerge(false);

  // Setup data structures for node intersection search (step 3) below.
  //
  // Search nodes must be ordered by begin position on each chromosome (this is an artifact of using a
  // non-general interval overlap test).
  typedef std::map<GenomeInterval, NodeIndexType> nodeMap_t;
  nodeMap_t                                       startLocusNodeMap;
  {
    const NodeIndexType nodeCount(startLocus.size());
    for (NodeIndexType nodeIndex(0); nodeIndex < nodeCount; ++nodeIndex) {
      startLocusNodeMap.insert(std::make_pair(startLocus.getNode(nodeIndex).getInterval(), nodeIndex));
    }
  }

  // reuse this intersectingNodeAddresses object throughout the merge:
  std::set<NodeAddressType> intersectingNodeAddresses;

  //
  // 3. Test if the graph's node intersections with any of startLocus' nodes exceeds complexity limits. If so,
  // abort the merge for startLocus.
  //
  // Note that the complexity measurement below is a function of individual startLocus nodes, but the
  // procedure has been setup to filter out the entire startLocus in case of even a single high-complexity
  // node. This could lead to non-complex/true variant graph edges being removed becuase they transitively
  // associate with a complex node.
  //
  // An improved procedure would identify only the high complexity nodes, and edit startLocus to remove them
  // (as well as any orphan nodes resulting from the complex node removal).
  //
  // Given the improved procedure would be difficult to implement, a stopgap measure is used below to
  // approximate this: the overwhelming majority of merge input are simple 2 node loci. If the complexity
  // check continues to check and filter these simple cases while skipping filtration of more complex loci,
  // there is no risk of throwing out a good edge by transitive association with a complex node.
  //
  const bool isTestUsability(inputLocus.size() <= 2);
  for (const nodeMap_t::value_type& startLocusNodeVal : startLocusNodeMap) {
    const bool isUsable(getIntersectingNodeAddresses(
        startLocusIndex, startLocusNodeVal.second, intersectingNodeAddresses, isTestUsability));

    if (!isUsable) {
#ifdef DEBUG_SVL
      log_os << logtag << "Aborting merge\n";
#endif
      isAbortMerge = true;
      break;
    }
  }

  //
  // 4. Test each node in the startLocus for mergeable intersections to other nodes in this graph. If such
  // cases are found, merge these nodes.
  //
  for (const nodeMap_t::value_type& startLocusNodeVal : startLocusNodeMap) {
    if (isAbortMerge) break;

    const NodeIndexType startLocusNodeIndex(startLocusNodeVal.second);

#ifdef DEBUG_SVL
    log_os << logtag
           << " startLocusNode: " << NodeAddressType(std::make_pair(startLocusIndex, startLocusNodeIndex))
           << " " << startLocus.getNode(startLocusNodeIndex);
#endif

    // get all nodes which can be merged with start node.
    getMergeableIntersectingNodeAddresses(
        startLocusIndex,
        startLocusNodeIndex,
        isStartLocusDuplicatedInAnotherGraphLocus,
        intersectingNodeAddresses);

#ifdef DEBUG_SVL
    log_os << logtag << " intersect_size: " << intersectingNodeAddresses.size() << "\n";
    for (const NodeAddressType& intersectingNodeAddress : intersectingNodeAddresses) {
      log_os << logtag << " intersect address: " << intersectingNodeAddress
             << " node: " << getNode(intersectingNodeAddress) << "\n";
    }
#endif

    // Skip to the next start locus node if no mergeable nodes are found.
    if (intersectingNodeAddresses.empty()) continue;

    // If the start locus has been copied to another graph locus, then the copy of the start locus will be
    // among the set of mergeable nodes returned. In this case we should also skip to the next node in the
    // start locus if only one intersecting node was found, because the start node is just finding a copy of
    // itself.
    if (isStartLocusDuplicatedInAnotherGraphLocus) {
      if (2 > intersectingNodeAddresses.size()) continue;
    }

    // If the intersecting nodes (and the start node they intersect) span multiple loci, consolidate all
    // intersecting nodes into a single locus in preparation for merging.
    //
    while (!isExpectedLocusOnly(headLocusIndex, intersectingNodeAddresses)) {
      // Since the intersecting nodes are not already in the same locus, the contents of
      // all loci represented in the intersecting node set are moved to a single locus.
      // By convention the target of this copy is the locus with the lowest index among the
      // loci in the intersecting node set.
      moveIntersectingNodesToLowestLocusIndex(intersectingNodeAddresses, startLocusIndex, headLocusIndex);

      // Record whether the start locus was duplicated in the above locus consolidation step
      if (!isStartLocusDuplicatedInAnotherGraphLocus)
        isStartLocusDuplicatedInAnotherGraphLocus = (headLocusIndex != startLocusIndex);

      // Rerun the search for nodes intersecting the start node:
      // 1. This is a lazy/safe way to update the intersecting nodes with their new addresses resulting from
      // the above locus consolidation step.
      // 2. It appears to be possible for this search to reveal new mergeable nodes, for which additional
      // rounds of locus consolidation are required.
      //
      // \TODO Document an example where this could occur
      getMergeableIntersectingNodeAddresses(
          startLocusIndex,
          startLocusNodeIndex,
          isStartLocusDuplicatedInAnotherGraphLocus,
          intersectingNodeAddresses);
      assert(!intersectingNodeAddresses.empty());

#ifdef DEBUG_SVL
      log_os << logtag << " multilocus detected, nodes moved and re-intersected. Updated intersect_size: "
             << intersectingNodeAddresses.size() << "\n";
#endif
    }

#ifdef DEBUG_SVL
    log_os << logtag << " intersect2_size: " << intersectingNodeAddresses.size() << "\n";
    for (const NodeAddressType& intersectingNodeAddress : intersectingNodeAddresses) {
      log_os << logtag << " intersect2 address: " << intersectingNodeAddress
             << " node: " << getNode(intersectingNodeAddress) << "\n";
    }
#endif

    // Find the node corresponding to the startLocusNode already duplicated into another locus in the graph.
    //
    // The startLocusNode may already have been merged into a larger node already, so we look for a node which
    // is a superset of startLocusNode.
    //
    NodeAddressType startLocusNodeSupersetAddress;
    {
      bool                    isStartLocusNodeSupersetFound(false);
      const known_pos_range2& startLocusNodeRange(
          getLocus(startLocusIndex).getNode(startLocusNodeIndex).getInterval().range);

      for (const NodeAddressType& intersectingNodeAddress : intersectingNodeAddresses) {
        assert(intersectingNodeAddress.first == headLocusIndex);

        // one node must be a superset of the start node, find this and store separately:
        if (getNode(intersectingNodeAddress).getInterval().range.is_superset_of(startLocusNodeRange)) {
          startLocusNodeSupersetAddress = intersectingNodeAddress;
          isStartLocusNodeSupersetFound = true;
          break;
        }
      }
      assert(isStartLocusNodeSupersetFound);
    }

    // Merge each intersecting node with the node found at startLocusNodeSupersetAddress. Target each merge to
    // the lower of the two node index values, and eliminate the node at the higher index value.
    //
    // Merge must be ordered from highest to lowest node index, so that the merge process does not invalidate
    // node indexes at higher value as it goes.
    //
    NodeAddressType mergeTargetNodeAddress(startLocusNodeSupersetAddress);
    BOOST_REVERSE_FOREACH(NodeAddressType intersectingNodeAddress, intersectingNodeAddresses)
    {
      if (intersectingNodeAddress == startLocusNodeSupersetAddress) continue;
      if (intersectingNodeAddress < mergeTargetNodeAddress)
        std::swap(intersectingNodeAddress, mergeTargetNodeAddress);
#ifdef DEBUG_SVL
      log_os << logtag << " MergeAndRemove: " << intersectingNodeAddress << "\n";
#endif
      mergeNodePtr(intersectingNodeAddress, mergeTargetNodeAddress);
      eraseNode(intersectingNodeAddress);
#ifdef DEBUG_SVL
      log_os << logtag << " Finished: " << intersectingNodeAddress << "\n";
      checkState();
#endif
    }
  }

  if (isAbortMerge || isStartLocusDuplicatedInAnotherGraphLocus) {
    clearLocus(startLocusIndex);
  }

#ifdef DEBUG_SVL
  checkState(true, true);
#endif
}

void SVLocusSet::merge(const SVLocusSet& inputSet)
{
  // TODO: check for compatible bam headers between inputSet and this

  assert(getMinMergeEdgeCount() == inputSet.getMinMergeEdgeCount());

  for (const SVLocus& locus : inputSet._loci) {
    if (locus.empty()) continue;

    try {
      merge(locus);
    } catch (...) {
      log_os << "ERROR: SVLocusSet merge failed.\n"
             << "\tSVLocusSet source: " << inputSet.getSource() << "\n"
             << "\tSVLocus index: " << locus.getIndex() << "\n";
      throw;
    }
  }

  _totalCleaned += inputSet._totalCleaned;
  _counts.merge(inputSet._counts);
  _highestSearchCount   = std::max(_highestSearchCount, inputSet._highestSearchCount);
  _isMaxSearchCount     = (_isMaxSearchCount || inputSet._isMaxSearchCount);
  _highestSearchDensity = std::max(_highestSearchDensity, inputSet._highestSearchDensity);
  _isMaxSearchDensity   = (_isMaxSearchDensity || inputSet._isMaxSearchDensity);
  _buildTime.merge(inputSet._buildTime);
  _mergeTime.merge(inputSet._mergeTime);  // this one is more of a formality...
}

bool SVLocusSet::getIntersectingNodeAddressesCore(
    const LocusIndexType       queryLocusIndex,
    const NodeIndexType        queryNodeIndex,
    const LocusSetIndexerType& searchNodes,
    const LocusIndexType       filterLocusIndex,
    std::set<NodeAddressType>& intersectingNodeAddresses,
    const bool                 isTestUsability) const
{
#ifdef DEBUG_SVL
  static const std::string logtag("SVLocusSet::getNodeIntersectCore");
  log_os << logtag << " inputNode: " << queryLocusIndex << ":" << queryNodeIndex << " "
         << getNode(std::make_pair(queryLocusIndex, queryNodeIndex));
  checkState();
#endif

  assert(_isIndexed);

  intersectingNodeAddresses.clear();

  // Get all nodes in searchNodes which intersect with the query node:
  const NodeAddressType queryNodeAddress(std::make_pair(queryLocusIndex, queryNodeIndex));
  const auto            searchNodeIterStart(searchNodes.data().lower_bound(queryNodeAddress));
  const GenomeInterval& queryInterval(getNode(queryNodeAddress).getInterval());
  const pos_t           maxRegionSize(_maxRegionSize[queryInterval.tid]);

  // diagnostics to determine if graph is growing too dense in one region:
  bool     isUsable(true);
  unsigned searchCount(0);

  // Look for all intersecting nodes with begin position >= queryInterval's begin position
  const auto searchNodeIterEnd(searchNodes.data().end());
  for (auto searchNodeIter(searchNodeIterStart); searchNodeIter != searchNodeIterEnd; ++searchNodeIter) {
    if (isTestUsability) {
      searchCount++;
      if (searchCount > _opt.maxSearchCount) {
        isUsable          = false;
        _isMaxSearchCount = true;
        break;
      }
    }

    const auto searchNodeLocusIndex(searchNodeIter->first);
    if (searchNodeLocusIndex == filterLocusIndex) continue;

#ifdef DEBUG_SVL
    log_os << logtag << "\tFWD test: " << (*searchNodeIter) << " " << getNode(*searchNodeIter);
#endif

    if (!queryInterval.isIntersect(getNode(*searchNodeIter).getInterval())) break;
    intersectingNodeAddresses.insert(*searchNodeIter);

#ifdef DEBUG_SVL
    log_os << logtag << "\tFWD insert: " << (*searchNodeIter) << "\n";
#endif
  }

  // Look for all intersecting nodes with begin position < queryInterval's begin position
  const auto searchNodeIterBegin(searchNodes.data().begin());
  for (auto searchNodeIter(searchNodeIterStart); searchNodeIter != searchNodeIterBegin;) {
    --searchNodeIter;

    if (isTestUsability) {
      if (!isUsable) break;
      searchCount++;
      if (searchCount > _opt.maxSearchCount) {
        isUsable          = false;
        _isMaxSearchCount = true;
        break;
      }
    }

    if (searchNodeIter->first == filterLocusIndex) continue;
#ifdef DEBUG_SVL
    log_os << logtag << "\tREV test: " << (*searchNodeIter) << " " << getNode(*searchNodeIter);
#endif
    const GenomeInterval& searchInterval(getNode(*searchNodeIter).getInterval());
    if (!queryInterval.isIntersect(searchInterval)) {
      if (!isOverlapAllowed()) break;

      if (queryInterval.tid != searchInterval.tid) break;
      if ((searchInterval.range.begin_pos() + maxRegionSize) < queryInterval.range.begin_pos()) break;
      continue;
    }

    intersectingNodeAddresses.insert(*searchNodeIter);
#ifdef DEBUG_SVL
    log_os << logtag << "\tREV insert: " << (*searchNodeIter) << "\n";
#endif
  }

  if (!isTestUsability) return true;

  _highestSearchCount = std::max(_highestSearchCount, searchCount);

  pos_t searchSize(
      queryInterval.range.end_pos() - std::max(0, queryInterval.range.begin_pos() - maxRegionSize));

  assert(searchSize >= 0);
  if (0 != searchSize) {
    static const pos_t minSearchSize(40);
    searchSize = std::max(searchSize, minSearchSize);
    const float searchDensity(static_cast<float>(searchCount) / static_cast<float>(searchSize));
    _highestSearchDensity = std::max(_highestSearchDensity, searchDensity);

    if (searchDensity > _opt.maxSearchDensity) {
      isUsable            = false;
      _isMaxSearchDensity = true;
    }
  }

  return isUsable;
}

void SVLocusSet::getIntersectingEdges(
    const LocusIndexType       queryLocusIndex,
    const NodeIndexType        queryRemoteNodeIndex,
    const EdgeMapType&         remoteIntersectNodeToLocalNodeMap,
    const LocusSetIndexerType& remoteIntersectNodes,
    std::vector<EdgeInfoType>& edges) const
{
  typedef EdgeMapType::const_iterator   rliter_t;
  typedef std::pair<rliter_t, rliter_t> rlmap_range_t;

  edges.clear();

  // Find all nodes, from the remoteIntersectNodes set, which intersect the query node:
  //
  // for this application, the query locus is an input set isolated from the rest of the graph, so nodes
  // intersected in the query locus are filtered out
  //
  std::set<NodeAddressType> edgeIntersectRemoteTemp;
  getIntersectingNodeAddressesCore(
      queryLocusIndex, queryRemoteNodeIndex, remoteIntersectNodes, queryLocusIndex, edgeIntersectRemoteTemp);

  for (const NodeAddressType& remoteIsectAddy : edgeIntersectRemoteTemp) {
    // find what local nodes the remote nodes trace back to:
    const rlmap_range_t remoteIsectRange(remoteIntersectNodeToLocalNodeMap.equal_range(remoteIsectAddy));
    assert(remoteIsectRange.first != remoteIntersectNodeToLocalNodeMap.end());
    for (rliter_t riter(remoteIsectRange.first); riter != remoteIsectRange.second; ++riter) {
      const NodeAddressType localIntersectAddy(std::make_pair(remoteIsectAddy.first, riter->second));
      edges.push_back(std::make_pair(localIntersectAddy, remoteIsectAddy.second));
    }
  }
}

void SVLocusSet::getIntersectingSignalNodeAddresses(
    const LocusIndexType             filterLocusIndex,
    const NodeAddressType            targetNodeAddress,
    std::set<NodeAddressType>&       intersectingSignalNodeAddresses,
    const std::set<NodeAddressType>& intersectingNoiseNodeTestTargets,
    bool&                            isIntersectingNoiseNodeOverlapTestTargets) const
{
#ifdef DEBUG_SVL
  static const std::string logtag("SVLocusSet::findSignalNodes");
  log_os << logtag << " targetNodeAddress: " << targetNodeAddress << "\n";
#endif
  // Get all nodes which intersect the node at targetNodeAddress.
  std::set<NodeAddressType> intersectingNodeAddresses;
  getIntersectingNodeAddressesCore(
      targetNodeAddress.first,
      targetNodeAddress.second,
      _inodes,
      filterLocusIndex,
      intersectingNodeAddresses);
  for (const NodeAddressType& intersectingNodeAddress : intersectingNodeAddresses) {
#ifdef DEBUG_SVL
    log_os << logtag << " intersectingNodeAddress: " << intersectingNodeAddress << "\n";
#endif
    if (!isNoiseNode(intersectingNodeAddress)) {
#ifdef DEBUG_SVL
      if (intersectingSignalNodeAddresses.count(intersectingNodeAddress) == 0) {
        log_os << logtag << " merge/new: " << targetNodeAddress << " " << intersectingNodeAddress << "\n";
      }
#endif
      intersectingSignalNodeAddresses.insert(intersectingNodeAddress);
    } else {
      if (!isIntersectingNoiseNodeOverlapTestTargets) {
        if (intersectingNoiseNodeTestTargets.count(intersectingNodeAddress)) {
          isIntersectingNoiseNodeOverlapTestTargets = true;
        }
      }
    }
  }
}

void SVLocusSet::getMergeableIntersectingNodeAddresses(
    const LocusIndexType       queryLocusIndex,
    const NodeIndexType        queryNodeIndex,
    const bool                 isQueryLocusDuplicatedInAnotherLocus,
    std::set<NodeAddressType>& mergeableIntersectingNodeAddresses) const
{
  //
  // TODO: There's room for significant optimization of these methods. The improvements are not trivial,
  //   but they would allow us to filter fewer nodes from being merged when node intersection counts become
  //   large.
  //

  //
  // There are two ways sets of mergeable nodes can occur:
  //
  // 1. There is a set of nodes which overlap with both query node and one
  // of the remote nodes that the query has an edge connecting to (ie they have a shared edge).
  // When totaled, the evidence count of these edges plus the query node edge
  // is at least minMergeEdgeCount. See schematic below.
  //
  // 2. The query node either contains an edge which is at least minMergeEdgeCount
  // or will contain such an edge due to (1), in this case the query node can be merged
  // with a locally overlapping node which also contains an edge which is greater than
  // minMergeEdgeCount. Note that in case (2) remote node intersection is not required.
  // See schematic below.
  //
  // # Schematic illustrations of intersection cases
  //
  // ## Key
  // - Assume a single genome coordinate in column space.
  // - Node is: |--NodeX--|
  // - Noise edge is:  ------>
  // - Signal edge is: ======>
  //
  //
  // ## Case 1:
  //
  // QueryLocus:      |-Node1-|--------------->|-Node2-|
  //
  // GraphLocus:           |-Node1-|-------->|-Node2-|
  //                           |-Node3-|--->|-Node4-|
  //
  // MergedLocus:     |-----Node1------|===>|--Node2---|
  //
  //
  // ## Case 2:
  //
  // QueryLocus:  |-Node1-|=======>|-Node2-|
  //
  // GraphLocus:                         |-Node1-|=======>|-Node2-|
  //
  // MergedLocus: |-Node3-|=======>|----Node1----|=======>|-Node2-|
  //

  const NodeAddressType queryNodeAddress(std::make_pair(queryLocusIndex, queryNodeIndex));
  const SVLocusNode&    queryNode(getNode(queryNodeAddress));

#ifdef DEBUG_SVL
  static const std::string logtag("SVLocusSet::getMergeableIntersectingNodeAddresses");
  log_os << logtag << " queryNode: " << queryNodeAddress << " " << queryNode;
  checkState();
#endif

  // Reuse this as a temporary throughout the methods below, it isn't actually needed at this scope.
  std::set<NodeAddressType> intersectingNodeAddresses;

  //
  // Build a new searchable node data structure which contains, for the set of nodes X intersecting the query
  // node,
  // all of the remote nodes Y connected by edges to X (searchableIntersectingNodeConnections)
  //
  // Also build a map for each node y \in Y pointing back to node x \in X (connectedNodeToIntersectingNodeMap)
  //
  // Schematic Example. Given the following queryNode and Graph Loci:
  //
  // queryNode:                 |--Node1--|
  //
  // GraphLocus1:                |-Node1-|<-------->|-Node2-|
  //                                |-Node3-|<--->|-Node4-|<---->| Node5 |
  //
  // GraphLocus2:  |-Node1-|<------->|-Node2-|
  //
  // searchableIntersectingNodeConnections should contain:
  // (Locus1,Node2),(Locus1,Node4),(Locus2,Node1)
  //
  // connectedNodeToIntersectingNodeMap should contain:
  // (Locus1,Node2) => (Locus1,Node1)
  // (Locus1,Node4) => (Locus1,Node3)
  // (Locus2,Node1) => (Locus2,Node2)
  //
  //
  LocusSetIndexerType searchableIntersectingNodeConnections(*this);
  EdgeMapType         connectedNodeToIntersectingNodeMap;

  // Nodes which intersect the query node and have already met the signal evidence threshold:
  std::set<NodeAddressType> intersectingSignalNodeAddresses;
  {
    // get a standard intersection of the query node:
    getIntersectingNodeAddresses(queryLocusIndex, queryNodeIndex, intersectingNodeAddresses);

    //
    // 1. Add data to searchableIntersectingNodeConnections and connectedNodeToIntersectingNodeMap
    //
    for (const NodeAddressType& intersectingNodeAddress : intersectingNodeAddresses) {
      const SVLocusNode& intersectingNode(getNode(intersectingNodeAddress));

      // Iterate through intersectingNode's edges to get its connecting node addresses
      const SVLocusEdgeManager intersectingNodeEdgeMap(intersectingNode.getEdgeManager());
      for (const SVLocusEdgesType::value_type& intersectingNodeEdge : intersectingNodeEdgeMap.getMap()) {
        NodeAddressType connectingNodeAddress(
            std::make_pair(intersectingNodeAddress.first, intersectingNodeEdge.first));
        searchableIntersectingNodeConnections.data().insert(connectingNodeAddress);
        connectedNodeToIntersectingNodeMap.insert(
            std::make_pair(connectingNodeAddress, intersectingNodeAddress.second));
      }
    }

#ifdef DEBUG_SVL
    log_os << logtag << " searchableIntersectingNodeConnections.size(): "
           << searchableIntersectingNodeConnections.data().size() << "\n";
    for (const NodeAddressType& nodeAddress : searchableIntersectingNodeConnections.data()) {
      log_os << logtag << "\tintersectingNodeConnection: " << nodeAddress << " " << getNode(nodeAddress);
    }
#endif

    //
    // 2. Add data to intersectingSignalNodeAddresses:
    //
    // Note that the signal node search is not transitive b/c we have required all signal nodes
    // in the graph to have merged already.
    //
    for (const NodeAddressType& intersectingNodeAddress : intersectingNodeAddresses) {
      if (!isNoiseNode(intersectingNodeAddress)) {
        intersectingSignalNodeAddresses.insert(intersectingNodeAddress);
      }
    }

#ifdef DEBUG_SVL
    log_os << logtag << " intersectingSignalNodeAddresses.size(): " << intersectingSignalNodeAddresses.size()
           << "\n";
    for (const NodeAddressType& intersectingSignalNodeAddress : intersectingSignalNodeAddresses) {
      log_os << logtag << "\tintersectingSignalNode: " << intersectingSignalNodeAddress << " "
             << getNode(intersectingSignalNodeAddress);
    }
#endif
  }

  //
  // Begin building this function's primary output, mergeableIntersectingNodeAddresses,
  // by enumerating all outgoing edges from the query node
  //
  mergeableIntersectingNodeAddresses.clear();

  // loop through each edge connected to the query node
  const SVLocusEdgeManager queryNodeOutgoingEdgeMap(queryNode.getEdgeManager());
  for (const SVLocusEdgesType::value_type& queryNodeOutgoingEdge : queryNodeOutgoingEdgeMap.getMap()) {
#ifdef DEBUG_SVL
    log_os << logtag << " processing edge: " << queryNodeAddress << "->" << queryLocusIndex << ":"
           << queryNodeOutgoingEdge.first << "\n";
    checkState();
#endif

    //
    // For each of the query node's outgoing edges, get all intersecting edges.
    //
    // 'intersecting edge' means that the from->to nodes connected by the intersecting edge overlap
    // with the from->to nodes (respectively) connected by the query node's outgoing edge
    //
    // Schematic Example. Given that Node1 is the queryNode and the following QueryLocus and Graph Loci:
    //
    // QueryLocus:                |--Node1--|<--------->|-Node2-|
    //
    // GraphLocus1:                |-Node1-|<-------->|-Node2-|
    //                |-Node4-|<---->|-Node3-|<--->|-Node4-|<------->| Node5 |
    //
    // GraphLocus2:  |-Node1-|<------->|-Node2-|
    //
    // The intersecting edges of the query node's only outgoing edge are:
    // (Locus1,Node1) -> (Locus1,Node2)
    // (Locus1,Node3) -> (Locus1,Node4)
    //
    std::vector<EdgeInfoType> intersectingEdges;
    getIntersectingEdges(
        queryLocusIndex,
        queryNodeOutgoingEdge.first,
        connectedNodeToIntersectingNodeMap,
        searchableIntersectingNodeConnections,
        intersectingEdges);

    unsigned intersectingEdgeCount(intersectingEdges.size());

    // If the query locus has not been copied into another graph locus, then the query node's outgoing edge
    // is not reflected in the intersecting edge count, so it needs to be added here to find the total count
    // of intersecting edges.
    if (!isQueryLocusDuplicatedInAnotherLocus) {
      intersectingEdgeCount++;
    }

    // isRegionCheck initiates a more detailed evidence signal threshold check process
    //
    // - The default process checks the total evidence summed over the entire intersecting edge set.
    // This neglects to account for the possibility that transitive overlap of many nodes could have
    // a high total evidence count but low evidence density.
    //
    // - The RegionCheck process sums up evidence in each genomic sub-interval. It more accurately
    // reflects peak evidence but is somewhat slower to compute.
    //
    // Example:
    //
    // Assume each node below has an evidence count of 1.
    //
    // |----Node1----|
    //           |-----Node2-----|
    //                        |----Node3----|
    //
    // Default evidence count:
    // 33333333333333333333333333333333333333
    //
    // RegionCheck evidence count:
    // 11111111112222211111111222211111111111
    //
    //

    // The peak evidence count from the RegionCheck process will always equal the evidence count from the
    // default process when 2 or fewer nodes exist, so there's no reason to turn the RegionCheck process on in
    // this case.
    const bool isRegionCheck(intersectingEdgeCount > 2);

    if (isRegionCheck) {
      _mergeRegions.clear();
    }

    // Store total evidence count for the outgoing and incoming edge as part of the default (ie.
    // non-RegionCheck) process to determine if the intersection set contains sufficient evidence to initiate
    // a merge
    //
    // Note: Outgoing/incoming orientation are defined relative to the target node.
    //
    unsigned mergedOutgoingEdgeEvidenceCount(0);
    unsigned mergedIncomingEdgeEvidenceCount(0);

    //
    // enumerate node evidence using either the default or RegionCheck process:
    //
    auto addEdgeEvidenceCount = [&](const SVLocus&      edgeLocus,
                                    const NodeIndexType edgeFromNodeIndex,
                                    const NodeIndexType edgeToNodeIndex) {
      // total edge counts on the to->from edge:
      const unsigned incomingEdgeEvidenceCount =
          edgeLocus.getEdge(edgeToNodeIndex, edgeFromNodeIndex).getCount();

      // total edge counts on the from->to edge:
      const unsigned outgoingEdgeEvidenceCount =
          edgeLocus.getEdge(edgeFromNodeIndex, edgeToNodeIndex).getCount();

      if (isRegionCheck) {
        // In the default process the local outgoing evidence and the remote incoming evidence are the same.
        // In the RegionCheck process the sum of all evidence for these two cases is still the same, but the
        // density of evidence in the overlapping local nodes and overlapping remote nodes may be very
        // different. For this reason the evidence density of both ends of both edges are enumerated below.
        //
        const known_pos_range2& localNodeRange(edgeLocus.getNode(edgeFromNodeIndex).getInterval().range);
        const known_pos_range2& remoteNodeRange(edgeLocus.getNode(edgeToNodeIndex).getInterval().range);

        _mergeRegions.localNodeOutgoingEdgeEvidence.add(localNodeRange, outgoingEdgeEvidenceCount);
        _mergeRegions.localNodeIncomingEdgeEvidence.add(localNodeRange, incomingEdgeEvidenceCount);
        _mergeRegions.remoteNodeOutgoingEdgeEvidence.add(remoteNodeRange, incomingEdgeEvidenceCount);
        _mergeRegions.remoteNodeIncomingEdgeEvidence.add(remoteNodeRange, outgoingEdgeEvidenceCount);
      } else {
        mergedOutgoingEdgeEvidenceCount += outgoingEdgeEvidenceCount;
        mergedIncomingEdgeEvidenceCount += incomingEdgeEvidenceCount;
      }
    };

    for (const EdgeInfoType& intersectingEdgeInfo : intersectingEdges) {
      const SVLocus&      intersectingEdgeLocus(getLocus(intersectingEdgeInfo.first.first));
      const NodeIndexType intersectingEdgeFromNodeIndex(intersectingEdgeInfo.first.second);
      const NodeIndexType intersectingEdgeToNodeIndex(intersectingEdgeInfo.second);
      addEdgeEvidenceCount(intersectingEdgeLocus, intersectingEdgeFromNodeIndex, intersectingEdgeToNodeIndex);
    }

    // If the query locus has not been duplicated in another graph locus, then the query node edge must
    // also be included to get accurate edge intersection count:
    if (!isQueryLocusDuplicatedInAnotherLocus) {
      addEdgeEvidenceCount(getLocus(queryLocusIndex), queryNodeIndex, queryNodeOutgoingEdge.first);
    }

    if (isRegionCheck) {
      // The density of evidence on the local and remote ends of an edge may be different, even though
      // the total evidence carried by the edge is the same. Here we take the lower peak evidence density from
      //  either of the local or remote nodes.
      //
      // Example:
      // In the following schematic there are 3 intersecting edges. The local node group is the set of nodes
      // intersecting Node1. The remote node is the set of nodes connected to local nodes by an edge. All 3
      // edges have an evidence count of 1 on the local->remote edge and 0 on the remote->local edge.
      //
      //       |-Node1-|----->|-Node2-|
      //        |-Node3-|---------->|-Node4-|
      //       |-Node5-|---------------->|-Node6-|
      //
      // In this case the default process would total the edge evidence count to 3. In the RegionCheck process
      // the local node outgoing edge evidence has a peak density of 3, but the remote node incoming edge
      // evidence has a peak density of 2. The min value of the two ends is used below, thus the final
      // outgoing edge evidence count is 2.
      //
      mergedOutgoingEdgeEvidenceCount = (std::min(
          _mergeRegions.localNodeOutgoingEdgeEvidence.maxVal(),
          _mergeRegions.remoteNodeIncomingEdgeEvidence.maxVal()));
      mergedIncomingEdgeEvidenceCount = (std::min(
          _mergeRegions.localNodeIncomingEdgeEvidence.maxVal(),
          _mergeRegions.remoteNodeOutgoingEdgeEvidence.maxVal()));
    }

#ifdef DEBUG_SVL
    log_os << logtag << " isRegionCheck: " << isRegionCheck << "\n";
    log_os << logtag << " final merge counts"
           << " local: " << mergedOutgoingEdgeEvidenceCount << " remote: " << mergedIncomingEdgeEvidenceCount
           << "\n";
    checkState();
#endif

    // if the total outgoing and incoming evidence count for this edge is below the signal threshold, then
    // move on to the query node's next edge.
    if ((mergedOutgoingEdgeEvidenceCount < getMinMergeEdgeCount()) &&
        (mergedIncomingEdgeEvidenceCount < getMinMergeEdgeCount()))
      continue;

    //
    // Add type1 mergeable nodes:
    //
    for (const EdgeInfoType& intersectingEdgeInfo : intersectingEdges) {
      mergeableIntersectingNodeAddresses.insert(intersectingEdgeInfo.first);
    }

    // for each type1 node, add any new intersections to the signal node set:
    //
    // this is not very efficient for now -- each type1 edge added in potentially
    // expands the current node to intersect new signal nodes
    // -- this loop looks for those new signal nodes
    //

    {
      // Get all of the remote node addresses from the intersecting edge set.
      //
      // This is used to search for the (rare) case where the intersecting edge set
      // locals overlap with the remotes.
      std::set<NodeAddressType> intersectingEdgeRemoteNodeAddresses;
      for (const EdgeInfoType& intersectingEdgeInfo : intersectingEdges) {
        intersectingEdgeRemoteNodeAddresses.insert(
            std::make_pair(intersectingEdgeInfo.first.first, intersectingEdgeInfo.second));
      }

      bool isIntersectRemotes(false);

      // Check both the query node and the nodes intersecting the query node for signal node intersections.
      //
      // Also check if any of these nodes intersect the remote node set (a rare event), in which case the
      // remote nodes will be added to the mergeable node set for this query node. This will mean that the
      // local and remote nodes will later be merged together into a larger node connected by a self-edge.
      //
      // Example of remote node intersection given Node1 is the query node:
      //
      // Intersecting Edges:
      //          |-Node1-|--->|-Node2-|
      //         |-Node3-|---------->|-Node4-|
      //               |--Node5--|------>|-Node6-|
      //
      // In this case Nodes 3 and 5 intersect the query node (Node1) and the set of three nodes has already
      // been found to be mergeable based on the total evidence of the three intersecting edges in the
      // left->right direction. The fact that one of the mergeable source (local) nodes (Node5) intersects one
      // of the sink (remote) nodes (Node2), is the condition that results in isIntersectRemotes being set to
      // true.
      //
      getIntersectingSignalNodeAddresses(
          queryLocusIndex,
          queryNodeAddress,
          intersectingSignalNodeAddresses,
          intersectingEdgeRemoteNodeAddresses,
          isIntersectRemotes);
      for (const EdgeInfoType& intersectingEdgeInfo : intersectingEdges) {
        getIntersectingSignalNodeAddresses(
            queryLocusIndex,
            intersectingEdgeInfo.first,
            intersectingSignalNodeAddresses,
            intersectingEdgeRemoteNodeAddresses,
            isIntersectRemotes);
      }

      if (isIntersectRemotes) {
        for (const NodeAddressType& intersectingEdgeRemoteNodeAddress : intersectingEdgeRemoteNodeAddresses) {
#ifdef DEBUG_SVL
          log_os << logtag << " adding ownRemote: " << intersectingEdgeRemoteNodeAddress << "\n";
#endif
          mergeableIntersectingNodeAddresses.insert(intersectingEdgeRemoteNodeAddress);

          // check to see if the remote nodes transitively intersect even more signal nodes!
          getIntersectingSignalNodeAddresses(
              queryLocusIndex,
              intersectingEdgeRemoteNodeAddress,
              intersectingSignalNodeAddresses,
              intersectingEdgeRemoteNodeAddresses,
              isIntersectRemotes);
        }
      }
    }

    //
    // Add type2 mergeable nodes:
    //
    for (const NodeAddressType& signalNodeAddress : intersectingSignalNodeAddresses) {
      mergeableIntersectingNodeAddresses.insert(signalNodeAddress);
    }
  }

#ifdef DEBUG_SVL
  log_os << logtag << " END. IntersectNodeSize: " << mergeableIntersectingNodeAddresses.size() << " Nodes:\n";
  for (const NodeAddressType addy : mergeableIntersectingNodeAddresses) {
    log_os << logtag << "\tInode: " << addy << "\n";
  }
#endif
}

void SVLocusSet::getRegionIntersect(const GenomeInterval interval, std::set<NodeAddressType>& intersectNodes)
{
  const LocusIndexType startLocusIndex(insertLocus(SVLocus()));
  const NodeIndexType  nodeIndex(getLocus(startLocusIndex).addNode(interval, this));

  getIntersectingNodeAddresses(startLocusIndex, nodeIndex, intersectNodes);

  clearLocus(startLocusIndex);
}

void SVLocusSet::moveIntersectingNodesToLowestLocusIndex(
    const std::set<NodeAddressType>& intersectingNodeAddresses,
    const LocusIndexType             startLocusIndex,
    LocusIndexType&                  headLocusIndex)
{
  // Capture the headLocusIndex value on input.
  const unsigned inputHeadLocusIndex(headLocusIndex);

  // Now reassign headLocusIndex to the lowest locus value found in the intersecting node set:
  bool isFirst(true);
  for (const NodeAddressType& intersectingNodeAddress : intersectingNodeAddresses) {
    if ((!isFirst) && (intersectingNodeAddress.first >= headLocusIndex)) continue;
    headLocusIndex = intersectingNodeAddress.first;
    isFirst        = false;
  }

  // Move all locus content to locus at headLocusIndex.
  //
  // Don't clear the source locus if it is startLocusIndex
  const bool isClearSourceLocus(startLocusIndex != inputHeadLocusIndex);
  combineLoci(inputHeadLocusIndex, headLocusIndex, isClearSourceLocus);
  for (const NodeAddressType& intersectingNodeAddress : intersectingNodeAddresses) {
    combineLoci(intersectingNodeAddress.first, headLocusIndex);
  }

#ifdef DEBUG_SVL
  static const std::string logtag("SVLocusSet::moveIntersectingNodesToLowestLocusIndex");
  log_os << logtag << " Reassigned all intersecting nodes to headLocusIndex: " << headLocusIndex
         << " inputHeadLocusIndex: " << inputHeadLocusIndex << " startLocusIndex:" << startLocusIndex << "\n";
  checkState();
#endif
}

void SVLocusSet::combineLoci(
    const LocusIndexType fromIndex, const LocusIndexType toIndex, const bool isClearSource)
{
  assert(toIndex < _loci.size());

#ifdef DEBUG_SVL
  static const std::string logtag("SVLocusSet::combineLoci");
  log_os << logtag << " from: " << fromIndex << " toIndex: " << toIndex << " isClear:" << isClearSource
         << "\n";
#endif

  if (fromIndex == toIndex) return;
  if (fromIndex >= _loci.size()) return;

  SVLocus& fromLocus(_loci[fromIndex]);
  if (fromLocus.empty()) return;

  SVLocus& toLocus(_loci[toIndex]);
  toLocus.copyLocus(fromLocus, this);
  if (isClearSource) clearLocus(fromIndex);
}

LocusIndexType SVLocusSet::insertLocus(const SVLocus& inputLocus)
{
  assert(_isIndexed);

  // The locus index (which is also the position in which we will place this locus in _loci), is found as
  // follows.
  //
  // If the locus vector contains any empty locus entries, we take the index from the lowest empty position.
  //
  // If the locus vector doesn't contain any empty entries, we append a new entry to the end of the vector,
  // and take the index reflecting the appended position.
  //
  LocusIndexType locusIndex(0);
  if (_emptyLoci.empty()) {
    static const unsigned maxIndex(std::numeric_limits<LocusIndexType>::max());
    locusIndex = _loci.size();
    assert(locusIndex < maxIndex);
    _loci.resize(locusIndex + 1);
  } else {
    locusIndex = (*_emptyLoci.begin());
    assert(_loci[locusIndex].empty());
    _emptyLoci.erase(locusIndex);
  }

  SVLocus& locus(_loci[locusIndex]);
  locus.updateIndex(locusIndex);
  locus.copyLocus(inputLocus, this);
  return locusIndex;
}

void SVLocusSet::mergeNodePtr(NodeAddressType fromPtr, NodeAddressType toPtr)
{
#ifdef DEBUG_SVL
  static const std::string logtag("SVLocusSet::mergeNodePtr");
  log_os << logtag << " from: " << fromPtr << " to: " << toPtr
         << " fromLocusSize: " << getLocus(fromPtr.first).size() << "\n";
#endif
  assert(_isIndexed);

  LocusSetIndexerType::iterator iter(_inodes.data().find(toPtr));
  assert(iter != _inodes.data().end());
  assert(fromPtr.first == toPtr.first);
  getLocus(fromPtr.first).mergeNode(fromPtr.second, toPtr.second, this);
}

void SVLocusSet::clean()
{
  for (SVLocus& locus : _loci) {
    if (locus.empty()) continue;
    _totalCleaned += locus.clean(getMinMergeEdgeCount(), this);

    // if true, this locus is newly empty after cleaning:
    if (locus.empty()) _emptyLoci.insert(locus.getIndex());
  }
#ifdef DEBUG_SVL
  checkForOverlapNodes(true);
#endif
}

void SVLocusSet::cleanRegion(const GenomeInterval interval)
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

void SVLocusSet::dump(std::ostream& os) const
{
  os << "LOCUSSET_START\n";
  for (const SVLocus& locus : _loci) {
    os << locus;
  }
  os << "LOCUSSET_END\n";
}

void SVLocusSet::dumpRegion(std::ostream& os, const GenomeInterval interval)
{
  std::set<NodeAddressType> intersectNodes;
  getRegionIntersect(interval, intersectNodes);

  LocusSetIndexerType sortedNodes(*this);
  for (const NodeAddressType& val : intersectNodes) {
    sortedNodes.data().insert(val);
  }

  for (const NodeAddressType& val : sortedNodes.data()) {
    os << "SVNode LocusIndex:NodeIndex : " << val << "\n";
    os << getNode(val);
  }
}

void SVLocusSet::dumpStats(std::ostream& os) const
{
  static const char sep('\t');

  os << "GraphBuildTime" << sep;
  _buildTime.reportHr(os);
  os << "\n";
  os << "GraphMergeTime" << sep;
  _mergeTime.reportHr(os);
  os << "\n";
  os << "disjointSubgraphs" << sep << nonEmptySize() << "\n";
  os << "nodes" << sep << totalNodeCount() << "\n";
  os << "directedEdges" << sep << totalEdgeCount() << "\n";
  os << "selfEdges" << sep << selfEdgeCount() << "\n";
  os << "totalGraphEvidence" << sep << totalObservationCount() << "\n";
  os << "totalCleaned" << sep << _totalCleaned << "\n";
  os << "highestSearchCount" << sep << _highestSearchCount << "\n";
  os << "isMaxSearchCount" << sep << _isMaxSearchCount << "\n";
  os << "highestSearchDensity" << sep << _highestSearchDensity << "\n";
  os << "isMaxSearchDensity" << sep << _isMaxSearchDensity << "\n";

  /// TODO: Add real sample labels:
  {
    std::vector<std::string> labels;
    const unsigned           csize(_counts.size());
    for (unsigned i(0); i < csize; ++i) {
      std::ostringstream oss;
      oss << "Sample" << i;
      labels.push_back(oss.str());
    }
    _counts.write(os, labels);
  }
  os << "\n";

  // node region size quantiles
  {
    SizeDistribution nodeSize;
    for (const SVLocus& locus : *this) {
      for (const SVLocusNode& node : locus) {
        const unsigned regionSize(node.getInterval().range.size());
        nodeSize.addObservation(regionSize);
      }
    }

    static const float    quantLevel[] = {0.25f, 0.5f, 0.75f, 0.9f, 0.95f, 0.99f};
    static const unsigned quantLevelCount(sizeof(quantLevel) / sizeof(float));
    os << "NodeRegionSizequantile:\n";
    for (unsigned i(0); i < quantLevelCount; ++i) {
      os << quantLevel[i] << sep << nodeSize.quantile(quantLevel[i]) << "\n";
    }
  }

  {
    // node edge count distro: 0.1,2,3... X+
    static const unsigned maxEdgeCount(10);
    std::vector<unsigned> edgeCount(maxEdgeCount);
    getNodeEdgeCountDistro(edgeCount);
    os << "NodeEdgeCount:\n";
    for (unsigned i(0); i < maxEdgeCount; ++i) {
      os << i;
      if ((i + 1) == maxEdgeCount) os << '+';
      os << sep << edgeCount[i] << "\n";
    }
  }

  {
    // node obs distro: 0,1,2,3... X+
    static const unsigned maxObsCount(30);
    std::vector<unsigned> obsCount(maxObsCount);
    getNodeObsCountDistro(obsCount);
    os << "NodeObservationCount:\n";
    for (unsigned i(0); i < maxObsCount; ++i) {
      os << i;
      if ((i + 1) == maxObsCount) os << '+';
      os << sep << obsCount[i] << "\n";
    }
  }
}

void SVLocusSet::dumpLocusStats(std::ostream& os) const
{
  static const char sep('\t');

  os << "locusIndex" << sep << "nodeCount" << sep << "nodeObsCount" << sep << "maxNodeObsCount" << sep
     << "regionSize" << sep << "maxRegionSize" << sep << "edgeCount" << sep << "maxEdgeCount" << sep
     << "edgeObsCount" << sep << "maxEdgeObsCount" << '\n';

  LocusIndexType locusIndex(0);
  for (const SVLocus& locus : _loci) {
    unsigned locusNodeObsCount(0), maxNodeObsCount(0);
    unsigned locusRegionSize(0), maxRegionSize(0);
    unsigned locusEdgeCount(0), maxEdgeCount(0), locusEdgeObsCount(0), maxEdgeObsCount(0);
    for (const SVLocusNode& node : locus) {
      // nodes:
      const unsigned nodeObsCount(node.outCount());
      maxNodeObsCount = std::max(maxNodeObsCount, nodeObsCount);
      locusNodeObsCount += nodeObsCount;

      // regions:
      const unsigned regionSize(node.getInterval().range.size());
      maxRegionSize = std::max(maxRegionSize, regionSize);
      locusRegionSize += regionSize;

      // edges:
      maxEdgeCount = std::max(maxEdgeCount, node.size());
      locusEdgeCount += node.size();
      const SVLocusEdgeManager edgeMap(node.getEdgeManager());
      for (const SVLocusEdgesType::value_type& edge : edgeMap.getMap()) {
        const unsigned edgeObsCount(edge.second.getCount());
        maxEdgeObsCount = std::max(maxEdgeObsCount, edgeObsCount);
        locusEdgeObsCount += edgeObsCount;
      }
    }
    os << locusIndex << sep << locus.size() << sep << locusNodeObsCount << sep << maxNodeObsCount << sep
       << locusRegionSize << sep << maxRegionSize << sep << locusEdgeCount << sep << maxEdgeCount << sep
       << locusEdgeObsCount << sep << maxEdgeObsCount << "\n";
    locusIndex++;
  }
}

void SVLocusSet::save(const char* filename) const
{
  using namespace boost::archive;

  assert(nullptr != filename);
  std::ofstream   ofs(filename, std::ios::binary);
  binary_oarchive oa(ofs);

  oa << getBamHeader();
  oa << _opt;
  oa << _isFinalized;
  oa << _totalCleaned;
  oa << _counts;
  oa << _highestSearchCount;
  oa << _highestSearchDensity;
  oa << _isMaxSearchCount;
  oa << _isMaxSearchDensity;
  oa << _buildTime;
  oa << _mergeTime;

  for (const SVLocus& locus : _loci) {
    if (locus.empty()) continue;
    oa << locus;
  }
}

SVLocusSet::SVLocusSet(const char* filename, const bool isSkipIndex)
  : SVLocusSet(SVLocusSetOptions(), bam_header_info(), {})
{
  using namespace boost::archive;

#ifdef DEBUG_SVL
  log_os << "SVLocusSet::load BEGIN\n";
#endif

  assert(filename);

  try {
    std::ifstream   ifs(filename, std::ios::binary);
    binary_iarchive ia(ifs);

    _source = filename;

    ia >> _bamHeaderInfo;
    ia >> _opt;
    ia >> _isFinalized;
    ia >> _totalCleaned;
    ia >> _counts;
    ia >> _highestSearchCount;
    ia >> _highestSearchDensity;
    ia >> _isMaxSearchCount;
    ia >> _isMaxSearchDensity;
    ia >> _buildTime;
    ia >> _mergeTime;

    SVLocus locus;
    while (ifs.peek() != EOF) {
      locus.clear(this);
      ia >> locus;
      if (locus.empty()) continue;
      const LocusIndexType locusIndex(size());
      _loci.push_back(locus);
      SVLocus& locusCopy(_loci.back());
      locusCopy.updateIndex(locusIndex);
    }
  } catch (...) {
    log_os << "ERROR: Exception caught while attempting to deserialize Manta SV locus graph file:\n"
           << "'" << filename << "'"
           << "\n";
    throw;
  }

  if (!isSkipIndex) {
    reconstructIndex();
    checkState(true, true);
  } else {
    _isIndexed = false;
  }

#ifdef DEBUG_SVL
  log_os << "SVLocusSet::load END\n";
#endif
}

void SVLocusSet::reconstructIndex()
{
#ifdef DEBUG_SVL
  log_os << "reconstructIndex BEGIN\n";
#endif
  clearIndex();

#ifdef DEBUG_SVL
  log_os << "reconstructIndex cleared\n";
#endif

  LocusIndexType locusIndex(0);
  for (SVLocus& locus : _loci) {
    const unsigned nodeCount(locus.size());
    for (NodeIndexType nodeIndex(0); nodeIndex < nodeCount; ++nodeIndex) {
      const NodeAddressType addy(std::make_pair(locusIndex, nodeIndex));
      _inodes.data().insert(addy);
      updateMaxRegionSize(getNode(addy).getInterval());
    }
    if (locus.empty()) _emptyLoci.insert(locusIndex);
    locusIndex++;
  }

  _isIndexed = true;

#ifdef DEBUG_SVL
  log_os << "reconstructIndex END\n";
#endif
}

#ifdef DEBUG_SVL
void SVLocusSet::dumpIndex(std::ostream& os) const
{
  assert(_isIndexed);

  os << "SVLocusSet Index START\n";
  for (const NodeAddressType& in : _inodes.data()) {
    os << "SVNodeIndex: " << in << "\n";
  }
  os << "SVLocusSet Index END\n";
}
#endif

void SVLocusSet::checkState(const bool isCheckOverlap, const bool isCheckLocusConnected) const
{
  using namespace illumina::common;

  assert(_isIndexed);

  unsigned locusIndex(0);
  unsigned checkStateTotalNodeCount(0);
  for (const SVLocus& locus : _loci) {
    locus.checkState(isCheckLocusConnected);

    const unsigned nodeCount(locus.size());
    checkStateTotalNodeCount += nodeCount;

    if (nodeCount == 0) {
      if (_emptyLoci.count(locusIndex) == 0) {
        std::ostringstream oss;
        oss << "Empty locus is not updated in the empty index. Locus index: " << locusIndex;
        BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
      }
    }

    for (NodeIndexType nodeIndex(0); nodeIndex < nodeCount; ++nodeIndex) {
      LocusSetIndexerType::const_iterator citer(_inodes.data().find(std::make_pair(locusIndex, nodeIndex)));
      if (citer == _inodes.data().end()) {
        std::ostringstream oss;
        oss << "Locus node is missing from node index\n"
            << "\tNode index: " << locusIndex << " node: " << getNode(std::make_pair(locusIndex, nodeIndex));
        BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
      }
      if ((citer->first != locusIndex) || (citer->second != nodeIndex)) {
        std::ostringstream oss;
        oss << "Locus node has conflicting index number in node index\n"
            << "\tinode index_value: " << citer->first << ":" << citer->second << "\n"
            << "\tNode index: " << locusIndex << ":" << locusIndex
            << " node: " << getNode(std::make_pair(locusIndex, nodeIndex));
        BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
      }
    }
    locusIndex++;
  }

  if (checkStateTotalNodeCount != _inodes.data().size()) {
    using namespace illumina::common;
    std::ostringstream oss;
    oss << "SVLocusSet conflicting internal node counts. TotalNodeCount: " << checkStateTotalNodeCount
        << " inodeSize: " << _inodes.data().size();
    BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
  }

  if (!isCheckOverlap) return;

  // if isOverlapAllowed() then we should expect noise nodes to overlap, but we can still check signal nodes:
  const bool isFilterNoise(isOverlapAllowed());
  checkForOverlapNodes(isFilterNoise);
}

#if 0
void
SVLocusSet::
compressSingletonNodes() const
{
    using namespace illumina::common;

    bool isFirst(true);
    GenomeInterval lastInterval;
    NodeAddressType lastAddy;
    for (const NodeAddressType& addy : _inodes.data())
    {
        if (isNoiseNode(addy)) continue;

        if (! isSingletonNode(addy)) continue;

        const GenomeInterval& interval(getNode(addy).getInterval());

        // don't allow zero-length or negative intervals:
        assert(interval.range.begin_pos() < interval.range.end_pos());

        // compress nearby singleton nodes into one:
        if (isFirst)
        {
            isFirst=false;
        }
        else if (interval.tid == lastInterval.tid)
        {
            if (lastInterval.range.end_pos() > interval.range.begin_pos())
            {
                std::ostringstream oss;
                oss << "Overlapping nodes in graph\n"
                    << "\tlast_index: " << lastAddy << " interval: " << lastInterval << "\n"
                    << "\tthis_index: " << addy << " interval: " << interval << "\n"
                    << "\tlast_node: " << lastAddy << " "<< getNode(lastAddy) << "\n"
                    << "\tthis_node: " << addy << " "<< getNode(addy) << "\n"
                    << "\n"
                    << header << "\n";
                BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
            }
        }
        lastAddy = addy;
        lastInterval = interval;
    }
}
#endif

void SVLocusSet::checkForOverlapNodes(const bool isFilterNoise) const
{
  using namespace illumina::common;

  bool            isFirst(true);
  GenomeInterval  lastInterval;
  NodeAddressType lastAddy;
  for (const NodeAddressType& addy : _inodes.data()) {
    if (isFilterNoise) {
      if (isNoiseNode(addy)) continue;
    }

    const GenomeInterval& interval(getNode(addy).getInterval());

    // don't allow zero-length or negative intervals:
    assert(interval.range.begin_pos() < interval.range.end_pos());

    // don't allow overlapping intervals:
    if (isFirst) {
      isFirst = false;
    } else if (interval.tid == lastInterval.tid) {
      if (lastInterval.range.end_pos() > interval.range.begin_pos()) {
        std::ostringstream oss;
        oss << "Overlapping nodes in graph\n"
            << "\tlast_index: " << lastAddy << " interval: " << lastInterval << "\n"
            << "\tthis_index: " << addy << " interval: " << interval << "\n"
            << "\tlast_node: " << lastAddy << " " << getNode(lastAddy) << "\n"
            << "\tthis_node: " << addy << " " << getNode(addy) << "\n"
            << "\n"
            << getBamHeader() << "\n";
        BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
      }
    }
    lastAddy     = addy;
    lastInterval = interval;
  }
}
