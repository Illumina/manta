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

#pragma once

#include "blt_util/flyweight_observer.hpp"
#include "svgraph/SVLocusNode.hpp"

typedef unsigned LocusIndexType;

/// The move message is composed of a bool, indicating if the node is being added (true) or deleted (false)
/// from the index, and the id of the node itself.
///
typedef std::pair<bool, std::pair<LocusIndexType, NodeIndexType>> SVLocusNodeMoveMessage;

struct SVLocusSet;

/// \brief a connected subgraph of the full SV locus graph
///
/// An SVLocus is a set of graph nodes which are all (transitivivly) associated by some level of sequencing
/// evidence, thus when the full SVLocus set is analyzed, each SVLocus can be treated as an independent
/// variant calling problem without introducing any approximations.
///
/// The SVLocus acts as a "container of SVLocusNode objects" and thus a significant portion of the API surface
/// is container-like (empty/size/begin/end). Note that edges are not explicitly stored in this graph format,
/// nodes contain all edge information.
///
/// The index number of this SVLocs is stored in the class and needs to be synced with the parent SVLocusSet
/// object.
///
/// The parent SVLocusSet object tracks individual SVLocusNode objects in this class (to enable a fast
/// region-based query of all nodes). Thus any operation which changes node index numbers has to be synced
/// with SVLocusSet via the flyweight_notifier base class.
///
/// The SVLocusNode index numbers also need to be synced with the node index number stored in the SVLocusNode
/// edges, all edges are recorded between nodes via node index numbers, so any change to the node index number
/// requires all in/out edges to be updated appropriately.
///
/// All SVLocusNodes will be non-overlapping when the graph is finalized, but during the graph building stage
/// this may not be the case.
///
struct SVLocus : public flyweight_notifier<SVLocusNodeMoveMessage> {
  typedef std::vector<SVLocusNode> graph_type;

  typedef graph_type::iterator       iterator;
  typedef graph_type::const_iterator const_iterator;

  friend struct SVLocusSet;

  bool empty() const { return _graph.empty(); }

  /// Total number of Locus Node
  unsigned size() const { return _graph.size(); }

  const_iterator begin() const { return _graph.begin(); }

  const_iterator end() const { return _graph.end(); }

  LocusIndexType getIndex() const { return _index; }

  /// Get a node by index
  const SVLocusNode& getNode(const NodeIndexType nodePtr) const
  {
#ifdef DEBUG_SVL
    if (nodePtr >= _graph.size()) nodeHurl(nodePtr);
#endif
    assert(nodePtr < _graph.size());
    return _graph[nodePtr];
  }

  /// Add a node
  /// \param[in] interval breakend interval
  /// \param[in] obs observer subscribes to the addNode events
  NodeIndexType addNode(const GenomeInterval interval, flyweight_observer_t* obs = nullptr)
  {
    assert(interval.tid >= 0);

    NodeIndexType nodePtr(newGraphNode());
    SVLocusNode&  node(getNode(nodePtr));
    node.setInterval(interval);
    // default _evidenceRange to the breakend interval unless a better estimate is provided
    node.setEvidenceRange(interval.range);
    notifyAdd(obs, nodePtr);
    return nodePtr;
  }

  /// Link two nodes
  /// \param[in] fromIndex index of from node
  /// \param[in] toIndex index of to node
  /// \param[in] fromCount edge count (weight) assigned to from -> to edge
  /// \param[in] toCount edge count (weight) assigned to to -> from edge
  /// by default, an edge count is only added on from -> to
  void linkNodes(
      const NodeIndexType fromIndex,
      const NodeIndexType toIndex,
      const unsigned      fromCount = 1,
      const unsigned      toCount   = 0)
  {
    SVLocusNode& fromNode(getNode(fromIndex));
    SVLocusNode& toNode(getNode(toIndex));
    assert(!fromNode.isEdge(toIndex));
    assert(!toNode.isEdge(fromIndex));

    SVLocusEdge fromEdge;
    fromEdge.setCount(fromCount);
    SVLocusEdge toEdge;
    toEdge.setCount(toCount);
    fromNode.mergeEdge(toIndex, fromEdge);
    toNode.mergeEdge(fromIndex, toEdge);
  }

  /// Set evidence range of a node
  void setNodeEvidence(const NodeIndexType nodeIndex, const known_pos_range2& evidenceRange)
  {
    getNode(nodeIndex).setEvidenceRange(evidenceRange);
  }

  /// Find the indices of all nodes connected to the \p startIndex node
  ///
  /// \param[in] startIndex Index of node to start connected node search from
  ///
  /// \param[out] connectedNodeIndexSet Container of indices for all nodes connected to \p startIndex node,
  /// including \p startIndex itself
  void findAllConnectedNodes(
      const NodeIndexType startIndex, std::set<NodeIndexType>& connectedNodeIndexSet) const;

  /// the total observations in all nodes of the locus
  unsigned totalObservationCount() const
  {
    unsigned sum(0);
    for (const SVLocusNode& node : *this) {
      sum += node.outCount();
    }
    return sum;
  }

  /// the total number of edges between all nodes of the locus
  unsigned totalEdgeCount() const
  {
    unsigned sum(0);
    for (const SVLocusNode& node : *this) {
      sum += node.size();
    }
    return sum;
  }

  /// the total number of self edges in the locus
  unsigned selfEdgeCount() const
  {
    unsigned       sum(0);
    const unsigned nodeSize(size());
    for (unsigned nodeIndex(0); nodeIndex < nodeSize; ++nodeIndex) {
      if (getNode(nodeIndex).isEdge(nodeIndex)) sum++;
    }
    return sum;
  }

  /// fill edge count histogram up to edgeCount.size()
  void getNodeEdgeCountDistro(std::vector<unsigned>& edgeCount) const
  {
    if (edgeCount.empty()) return;
    const unsigned maxEdge(edgeCount.size() - 1);
    for (const SVLocusNode& node : *this) {
      edgeCount[std::min(node.size(), maxEdge)]++;
    }
  }

  /// fill obs count histogram up to obsCount.size()
  void getNodeObsCountDistro(std::vector<unsigned>& obsCount) const
  {
    if (obsCount.empty()) return;
    const unsigned maxObs(obsCount.size() - 1);
    for (const SVLocusNode& node : *this) {
      obsCount[std::min(node.outCount(), maxObs)]++;
    }
  }

  /// return from->to edge
  const SVLocusEdge& getEdge(const NodeIndexType fromNodeIndex, const NodeIndexType toNodeIndex) const
  {
    const SVLocusNode& fromNode(getNode(fromNodeIndex));
    try {
      return fromNode.getEdge(toNodeIndex);
    } catch (...) {
      // throw a richer exception message than node can produce on its own:
      getEdgeException(fromNodeIndex, toNodeIndex);
    }

    // handle compiler warning for return val, this code should never run:
    static SVLocusEdge bogusWarning;
    return bogusWarning;
  }

  /// clear the graph
  void clear(flyweight_observer_t* obs)
  {
    for (NodeIndexType i(0); i < size(); ++i) {
      notifyDelete(obs, i);
    }
    _graph.clear();
  }

  /// find any self-overlapping nodes within the locus and merge
  void mergeSelfOverlap();

  /// Assert that internal data-structures are in a consistent state
  ///
  /// \param[in] isCheckConnected If true, check that all nodes in the locus are connected.
  void checkState(const bool isCheckConnected = false) const;

  /// total the evidence count of all in-edges to this node
  unsigned getNodeInCount(const LocusIndexType nodeIndex) const;

  void clearNodeEdges(const NodeIndexType nodePtr);

  /// a fancier version of the SVLocusNode dumper which can
  /// report in-edge information
  void dumpNode(std::ostream& os, const LocusIndexType nodeIndex) const;

  template <class Archive>
  void save(Archive& ar, const unsigned /* version */) const
  {
    ar << _graph;
  }

  template <class Archive>
  void load(Archive& ar, const unsigned /* version */)
  {
    clear(nullptr);
    ar >> _graph;
  }

  BOOST_SERIALIZATION_SPLIT_MEMBER()

private:
  SVLocusNode& getNode(const NodeIndexType nodePtr)
  {
#ifdef DEBUG_SVL
    if (nodePtr >= _graph.size()) nodeHurl(nodePtr);
#endif
    assert(nodePtr < _graph.size());
    return _graph[nodePtr];
  }

  void nodeHurl(const NodeIndexType nodePtr) const;

  void updateIndex(const LocusIndexType& index) { _index = index; }

  /// Test if all (out+in) edges of the indexed node have less than \p minMergeEdgeCount evidence.
  /// If this is true the node will be removed when the graph is denoised, we refer to this as a
  /// noise node.
  ///
  /// \param[in] minMergeEdgeCount The minimum evidence on an edge below which it is considered noise.
  /// \param[in] nodeIndex Index of the query node within this locus.
  /// \return True if the index points to a noise node
  bool isNoiseNode(const unsigned minMergeEdgeCount, const NodeIndexType nodeIndex) const;

  /// remove all unmerged noise in-edges of node and
  /// provide list of nodes which should be deleted
  ///
  /// return amount of evidence cleaned
  unsigned cleanNodeCore(
      const unsigned minMergeEdgeCount, const NodeIndexType nodeIndex, std::set<NodeIndexType>& emptyNodes);

  /// remove all unmerged noise in-edges of node and possibly
  /// delete empty nodes
  ///
  /// return amount of evidence cleaned
  unsigned cleanNode(
      const unsigned minMergeEdgeCount, const NodeIndexType nodeIndex, flyweight_observer_t* obs);

  /// remove all unmerged noise edges and nodes
  ///
  /// return amount of evidence cleaned
  unsigned clean(const unsigned minMergeEdgeCount, flyweight_observer_t* obs);

  void getEdgeException(const NodeIndexType fromIndex, const NodeIndexType toIndex) const;

  /// erase edges in both directions:
  void eraseEdgePair(const NodeIndexType index1, const NodeIndexType index2)
  {
    eraseEdge(index1, index2);
    if (index1 == index2) return;
    eraseEdge(index2, index1);
  }

  /// erase edge in one direction
  void eraseEdge(const NodeIndexType fromIndex, const NodeIndexType toIndex)
  {
    getNode(fromIndex).eraseEdge(toIndex);
  }

  /// \brief Copy the contents of \p fromLocus into this locus without merging.
  ///
  /// This method copies each node from \p fromLocus into this locus object. The method
  /// makes no attempt to merge overlapping nodes, it is meant to be run as an intermediate
  /// step in the merge process.
  ///
  /// The fromLocus nodes that are copied in are mostly unchanged, except that all node index
  /// numbers used have to be offset such that they are all higher than the node index numbers
  /// already in use by this locus object.
  ///
  /// \param obs Observer object is required to notify the parent SVLocusSet of all node changes so that its
  /// 'search for nodes by genome interval' logic works correctly.
  void copyLocus(const SVLocus& fromLocus, flyweight_observer_t* obs)
  {
    assert(&fromLocus != this);

    // simple method -- copy everything in with an offset in all index numbers:
    const unsigned offset(_graph.size());
    for (const SVLocusNode& fromNode : fromLocus) {
      const NodeIndexType nodeIndex(newGraphNode());
      getNode(nodeIndex) = SVLocusNode(fromNode, offset);
      notifyAdd(obs, nodeIndex);
    }
  }

  /// join from node into to node
  ///
  /// from node is effectively destroyed,
  //// because all incoming edges will be updated
  ///
  void mergeNode(const NodeIndexType fromIndex, const NodeIndexType toIndex, flyweight_observer_t* obs);

  /// \brief Remove node \p nodePtr from this Locus
  ///
  void eraseNode(const NodeIndexType nodePtr, flyweight_observer_t* obs);

  // remove a list of node ids
  void eraseNodes(const std::set<NodeIndexType>& nodes, flyweight_observer_t* obs);

  NodeIndexType newGraphNode()
  {
    static const unsigned maxIndex(std::numeric_limits<NodeIndexType>::max());
    unsigned              index(_graph.size());
    assert(index < maxIndex);
    _graph.resize(index + 1);
    return static_cast<NodeIndexType>(index);
  }

  void notifyAdd(flyweight_observer_t* obs, const NodeIndexType nodePtr)
  {
    if (nullptr == obs) return;
#ifdef DEBUG_SVL
    log_os << "SVLocusNotifier: Add node: " << _index << ":" << nodePtr << "\n";
#endif
    notify_flyweight_observer(obs, std::make_pair(true, std::make_pair(_index, nodePtr)));
  }

  void notifyDelete(flyweight_observer_t* obs, const NodeIndexType nodePtr)
  {
    if (nullptr == obs) return;
#ifdef DEBUG_SVL
    log_os << "SVLocusNotifier: Delete node: " << _index << ":" << nodePtr << "\n";
#endif
    notify_flyweight_observer(obs, std::make_pair(false, std::make_pair(_index, nodePtr)));
  }

  graph_type     _graph;
  LocusIndexType _index = 0;
};

std::ostream& operator<<(std::ostream& os, const SVLocus& locus);

BOOST_CLASS_IMPLEMENTATION(SVLocus, boost::serialization::object_serializable)
