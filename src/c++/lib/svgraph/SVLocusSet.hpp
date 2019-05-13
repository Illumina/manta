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

#include <algorithm>
#include <iosfwd>
#include <string>
#include <vector>

#include "blt_util/RegionSum.hpp"
#include "blt_util/time_util.hpp"
#include "htsapi/bam_header_info.hpp"
#include "manta/SVBreakend.hpp"
#include "options/SVLocusSetOptions.hpp"
#include "svgraph/SVLocus.hpp"
#include "svgraph/SVLocusSampleCounts.hpp"

#ifdef DEBUG_SVL
#include <iostream>
#include "blt_util/log.hpp"
#endif

/// \brief The parent object used to manage Manta's SV Locus graph.
///
/// This object includes all data on the graph itself, in addition to various meta-data on the graph or graph
/// creation process.
///
/// The primary role of this object, per its name, is to act as a "bag of SVLocus objects", thus there is some
/// standard container interface supporting this role (size/begin/end...), and special interface reflecting
/// the region based nature of the stored objects (getRegionIntersect). The rationale of the SVLocus object is
/// that one should exist for each disjoint subgraph of the SV Locus graph. This is true when the graph is
/// finalized, but may not hold during intermediate stages of the graph building and merging process.
///
/// Additional interface methods added to this object include:
///
/// merge() - Merging 2 SVLocusSet objects is the fundamental operation which allows for the creation of a
/// genome-wide graph. It is one of the more complex operations in Manta due to the signal vs. noise edge
/// evidence distinction.
///
/// cleanRegion()/clean()/finalize() - 'cleaning' is an operation to remove noise edges. For the merge
/// operation to work without information loss, this cannot be done until all sequences from a region have
/// been observed and merged in, after which any noise edge connecting two nodes which both lie in the
/// observed region can be removed. cleanRegion() does this for a specific region (executed as information is
/// scanned in for the purpose of memory management), while clean() removes all noise edges from the entire
/// graph. finalize() is used to declare intent to stop merging information into the graph and begin variant
/// calling, its primary role is to run clean().
///
/// save()/load() - write/read the graph
///
/// Various additional stats functions report on:
/// 1. graph summary stats (no of nodes, edges, self-edges, etc...)
/// 2. graph creation summary stats (no of anom/split/filtered read-pairs observed, etc...)
/// 3. graph processing info (no of edges/nodes cleaned so far), reading merging runtime cost,...
/// 4. object parameters
///
/// ## Region based node query
///
/// Internal methods mostly provide details to support the above interface. Of special note is the node
/// indexing system used by this object. For certain operations we (internally) need to provide range based
/// queries for nodes intersecting points or segments of the genome. While an interval tree is an obvious
/// starting point for this, it is pessimistic in time/memory, given that nodes are very small relative to
/// chrom size. Instead, the maximum node size is stored for every chromosome, then a binary search is used to
/// search nodes indexed on their beginning positions, and linear search backwards up to the maximum node size
/// guarantees that all nodes intersecting a point can be found in (effectively) log-N time.
///
/// This object inherits from flyweight_observer to help keep this node indexing data structure up to date as
/// nodes are added/deleted within the SVLocus objects managed by this graph.
///
struct SVLocusSet : public flyweight_observer<SVLocusNodeMoveMessage>, private boost::noncopyable {
  typedef std::vector<SVLocus>          locusset_type;
  typedef locusset_type::const_iterator const_iterator;

  /// For test purposes SVLocusSet can be constructed with no arguments,
  /// but in this case it will have no chromosome or sample information.
  explicit SVLocusSet(
      const SVLocusSetOptions&        opt                = SVLocusSetOptions(),
      const bam_header_info&          bamHeaderInfo      = bam_header_info(),
      const std::vector<std::string>& alignmentFilenames = {});

  /// \brief Deserialize object from binary file format
  /// \param[in] isSkipIndex If true, don't build the graph index, and only allow a limited set of operations
  ///
  explicit SVLocusSet(const char* filename, const bool isSkipIndex = false);

  /// Test if the set of SVLocus objects is empty.
  bool empty() const { return _loci.empty(); }

  /// Count the number of SVLocus objects in the set.
  unsigned size() const { return _loci.size(); }

  /// Count the number of SVLocus objects in the set which are not empty.
  unsigned nonEmptySize() const
  {
    assert(_isIndexed);
    return size() - _emptyLoci.size();
  }

  /// Return an iterator for the SVLocus set pointing at the beginning.
  const_iterator begin() const { return _loci.begin(); }

  /// Return an iterator for the SVLocus set pointing at the end.
  const_iterator end() const { return _loci.end(); }

  /// Return a reference of the SVLocus at the given index.
  const SVLocus& getLocus(const LocusIndexType index) const
  {
#ifdef DEBUG_SVL
    if (index >= _loci.size()) locusHurl(index, "const");
#endif
    assert(index < _loci.size());
    return _loci[index];
  }

  /// \brief Merge SVLocus into the SVLocusSet.
  ///
  /// The referenced SVLocus is destroyed in this process.
  void merge(const SVLocus& locus);

  /// \brief Merge referenced SVLocusSet into the calling SVLocusSet.
  ///
  /// The referenced SVLocusSet is destroyed in this process.
  void merge(const SVLocusSet& set);

  /// Indicate that the set is complete
  void finalize()
  {
    clean();
    _isFinalized = true;
  }

  /// Remove all existing edges with less than minMergeEdgeCount support:
  void clean();

  /// Remove all existing edges with less than minMergeEdgeCount
  void cleanRegion(const GenomeInterval interval);

  /// Return the number of nodes that have been removed from Locus objects by the clean and cleanRegion
  /// operations
  unsigned totalCleaned() const { return _totalCleaned; }

  /// Binary serialization
  void save(const char* filename) const;

  /// Debug output.
  void dump(std::ostream& os) const;

  /// \brief Debug output for a specific region.
  ///
  /// Note this is conceptually const but can't be marked const to the compiler right now due to
  /// implementation details (see getRegionIntersect).
  void dumpRegion(std::ostream& os, const GenomeInterval interval);

  /// Debug stats output on the whole SVLocus set.
  void dumpStats(std::ostream& os) const;

  /// Debug stats on each locus in tsv format.
  void dumpLocusStats(std::ostream& os) const;

  /// Return a debug string of the source of the SVLocusSet.
  const std::string& getSource() const { return _source; }

  /// Return the setting for the minMergeEdgeCount.
  unsigned getMinMergeEdgeCount() const { return _opt.getMinMergeEdgeCount(); }

  /// Total number of reads used as supporting evidence in the graph
  unsigned totalObservationCount() const
  {
    unsigned sum(0);
    for (const SVLocus& locus : *this) {
      sum += locus.totalObservationCount();
    }
    return sum;
  }

  /// Return the total number of nodes in the SVLocus objects.
  unsigned totalNodeCount() const
  {
    unsigned sum(0);
    for (const SVLocus& locus : *this) {
      sum += locus.size();
    }
    return sum;
  }

  /// Return the total number of directed edges in the graph
  unsigned totalEdgeCount() const
  {
    unsigned sum(0);
    for (const SVLocus& locus : *this) {
      sum += locus.totalEdgeCount();
    }
    return sum;
  }

  /// Return the total number of self-edges in the graph
  unsigned selfEdgeCount() const
  {
    unsigned sum(0);
    for (const SVLocus& locus : *this) {
      sum += locus.selfEdgeCount();
    }
    return sum;
  }

  /// Fill node edge count histogram up to edgeCount.size()
  void getNodeEdgeCountDistro(std::vector<unsigned>& edgeCount) const
  {
    for (const SVLocus& locus : *this) {
      locus.getNodeEdgeCountDistro(edgeCount);
    }
  }

  /// Fill node observation count histogram up to obsCount.size()
  void getNodeObsCountDistro(std::vector<unsigned>& obsCount) const
  {
    for (const SVLocus& locus : *this) {
      locus.getNodeObsCountDistro(obsCount);
    }
  }

  /// Check that internal data-structures are in
  /// a consistent state, throw on error
  void checkState(const bool isCheckOverlap = false, const bool isCheckLocusConnected = false) const;

  const AllSampleReadCounts& getAllSampleReadCounts() const { return _counts; }

  SampleReadInputCounts& getSampleReadInputCounts(const unsigned index)
  {
    return _counts.getSampleCounts(index).input;
  }

  SampleEvidenceCounts& getSampleEvidenceCounts(const unsigned index)
  {
    return _counts.getSampleCounts(index).evidence;
  }

  /// Set the buildTime  SVLocusSet object
  void addBuildTime(const CpuTimes& t) { _buildTime.merge(t); }

  /// Set the mergeTime SVLocusSet object.
  void addMergeTime(const CpuTimes& t) { _mergeTime.merge(t); }

  typedef std::pair<LocusIndexType, NodeIndexType> NodeAddressType;

  /// \brief Fill \p intersectNodes with all node addresses in this object which intersect with \p interval
  ///
  /// \param[in] interval Target genome region for node intersection search
  /// \param[out] intersectNodes Intersecting nodes are put in this object. This object is cleared on input.
  ///
  /// This method is conceptually const. For programming convenience it temporarily changes the graph
  /// structure (then undoes the change), so it can't be marked const to the compiler.
  void getRegionIntersect(const GenomeInterval interval, std::set<NodeAddressType>& intersectNodes);

  /// |brief Provide const access the bam header info which this object stores
  const bam_header_info& getBamHeader() const { return _bamHeaderInfo; }

private:
  typedef NodeAddressType EdgeMapKeyType;
  typedef NodeIndexType   EdgeMapValueType;

  typedef std::multimap<EdgeMapKeyType, EdgeMapValueType> EdgeMapType;

  typedef std::pair<EdgeMapKeyType, EdgeMapValueType> EdgeInfoType;

  /// \brief  NodeAddressSorter wraps a SVLocusSet in order to add support for
  /// testing if 2 NodeAddressTypes are sorted a < b.
  struct NodeAddressSorter {
    NodeAddressSorter(const SVLocusSet& set) : _set(set) {}

    bool operator()(const NodeAddressType& a, const NodeAddressType& b) const
    {
      if (getInterval(a) < getInterval(b)) return true;
      if (getInterval(a) == getInterval(b)) {
        // If the GenomeIntervals are the same, compare the contents of the NodeAddressTypes directly.
        return (a < b);
      }
      return false;
    }

  private:
    const GenomeInterval& getInterval(const NodeAddressType& n) const
    {
      return (_set.getLocus(n.first).getNode(n.second).getInterval());
    }

    const SVLocusSet& _set;
  };

  /// \brief Container to hold a set of node addresses which support range-based node intersect queries.
  ///
  /// This object holds a set of node addresses with a sorting scheme designed to support
  /// range-based intersection queries.
  ///
  /// The core data type here, data_t, is wrapped by this object because this object requires a custom
  /// copy-ctor/assign-op, whereas data_t can use defaults.
  struct LocusSetIndexerType {
    typedef std::set<NodeAddressType, NodeAddressSorter> data_t;
    typedef data_t::iterator                             iterator;
    typedef data_t::const_iterator                       const_iterator;

    LocusSetIndexerType(const SVLocusSet& set) : _data(NodeAddressSorter(set)) {}

    LocusSetIndexerType(const LocusSetIndexerType& rhs) = delete;

    LocusSetIndexerType& operator=(const LocusSetIndexerType& rhs)
    {
      if (this == &rhs) return *this;
      _data.clear();
      _data.insert(rhs._data.begin(), rhs._data.end());
      return *this;
    }

    data_t& data() { return _data; }

    const data_t& data() const { return _data; }

  private:
    data_t _data;
  };

  friend std::ostream& operator<<(std::ostream& os, const NodeAddressType& a);

  SVLocus& getLocus(const LocusIndexType index)
  {
#ifdef DEBUG_SVL
    if (index >= _loci.size()) locusHurl(index, "non-const");
#endif

    assert(index < _loci.size());
    return _loci[index];
  }

#ifdef DEBUG_SVL
  void locusHurl(const LocusIndexType index, const char* label) const;
#endif

  const SVLocusNode& getNode(const NodeAddressType n) const { return getLocus(n.first).getNode(n.second); }

  void clearLocus(const LocusIndexType index)
  {
#ifdef DEBUG_SVL
    log_os << "SVLocusSet::clearLocus index: " << index << "\n";
#endif
    assert(index < _loci.size());

    _loci[index].clear(this);
    _emptyLoci.insert(index);
    _source = "UNKNOWN";
  }

  /// \brief Get addresses of all nodes in the graph which intersect with the query node. (A more general
  /// version of getNodeIntersect)
  ///
  /// \param[in] queryLocusIndex First part of query node address.
  ///
  /// \param[in] queryNodeIndex Second part of query node address.
  ///
  /// \param[in] searchNodes The set of nodes which will be searched for intersections with the query node.
  ///
  /// \param[in] filterLocusIndex Intersections to nodes in this locus will be filtered out of the results.
  ///
  /// \param[out] intersectingNodeAddresses Set of all intersecting node addresses. Any set contents are
  /// erased on input.
  ///
  /// \param[in] isTestUsability If true, check whether a node intersection exceeds computability limits.
  ///
  /// \return True if the query node is usable. This can only be false when isTestUsability is true.
  bool getIntersectingNodeAddressesCore(
      const LocusIndexType       queryLocusIndex,
      const NodeIndexType        queryNodeIndex,
      const LocusSetIndexerType& searchNodes,
      const LocusIndexType       filterLocusIndex,
      std::set<NodeAddressType>& intersectingNodeAddresses,
      const bool                 isTestUsability = false) const;

  /// \brief Get addresses of all nodes in the graph which intersect with the query node.
  ///
  /// \param[in] queryLocusIndex First part of query node address.
  ///
  /// \param[in] queryNodeIndex Second part of query node address.
  ///
  /// \param[out] intersectingNodeAddresses Set of all intersecting node addresses. Any set contents are
  /// erased on input.
  ///
  /// \param[in] isTestUsability If true, check whether a node intersection exceeds computability limits.
  ///
  /// \return True if the query node is usable. This can only be false when isTestUsability is true.
  bool getIntersectingNodeAddresses(
      const LocusIndexType       queryLocusIndex,
      const NodeIndexType        queryNodeIndex,
      std::set<NodeAddressType>& intersectingNodeAddresses,
      const bool                 isTestUsability = false) const
  {
    return getIntersectingNodeAddressesCore(
        queryLocusIndex,
        queryNodeIndex,
        _inodes,
        queryLocusIndex,
        intersectingNodeAddresses,
        isTestUsability);
  }

  /// edges returned are in local_addy->remote_node orientation
  void getIntersectingEdges(
      const LocusIndexType       queryLocusIndex,
      const NodeIndexType        queryRemoteNodeIndex,
      const EdgeMapType&         remoteIntersectNodeToLocalNodeMap,
      const LocusSetIndexerType& remoteIntersectNodes,
      std::vector<EdgeInfoType>& edges) const;

  /// \brief Get addresses of nodes which could be merged with the query node, accounting for edge overlap and
  /// signal/noise thresholds used in the merging process.
  ///
  /// This will return only a subset of the nodes which intersect the query node.
  ///
  /// \param[in] queryLocusIndex First part of query node address.
  /// \param[in] queryNodeIndex Second part of query node address.
  /// \param[in] isQueryLocusDuplicatedInAnotherLocus True if the query locus been copied into a graph locus
  /// besides \p queryLocusIndex.
  /// \param[out] mergeableIntersectingNodeAddresses Addresses of nodes which could be merged with the query
  /// node. Any set contents are erased on input.
  void getMergeableIntersectingNodeAddresses(
      const LocusIndexType       queryLocusIndex,
      const NodeIndexType        queryNodeIndex,
      const bool                 isQueryLocusDuplicatedInAnotherLocus,
      std::set<NodeAddressType>& mergeableIntersectingNodeAddresses) const;

  /// \brief Move all intersecting nodes to be co-located in the same locus.
  ///
  /// The set of intersecting loci are the loci containing nodes in the intersecting node set. In this method
  /// the contents of all intersecting loci are copied to a target locus. The target locus is the member of
  /// the intersecting loci set with the lowest locus index.
  ///
  /// All loci are cleared after copying into the target locus, except for the locus at \p startLocusIndex.
  ///
  /// By co-locating all of the intersecting nodes into one locus, this method achieves a preliminary step in
  /// the process of merging the intersecting node set.
  ///
  /// \param[in] intersectingNodeAddresses Addresses of a set of nodes which will be moved such that they are
  /// all in the same locus.
  /// \param[in] startLocusIndex This locus is not cleared after copying its contents to the lowest locus
  /// index.
  /// \param[in,out] headLocusIndex The value on input should contain the locus holding the duplicated start
  /// locus.If the start locus hasn't been duplicated yet it will equal \p startLocusIndex. The value
  /// returned is the target locus of the move operation, which after this function completes is the neow
  /// locus holding the duplication of the start locus.
  void moveIntersectingNodesToLowestLocusIndex(
      const std::set<NodeAddressType>& intersectingNodeAddresses,
      const LocusIndexType             startLocusIndex,
      LocusIndexType&                  headLocusIndex);

  /// \brief Combine all content from \p fromIndex locus into \p toIndex locus
  ///
  /// This is typically required as a first step to merging a set of nodes which span two different loci.
  ///
  /// \param[in] isClearSource If true, clear the locus at \p fromIndex after copying.
  void combineLoci(
      const LocusIndexType fromIndex, const LocusIndexType toIndex, const bool isClearSource = true);

  /// \brief Add \p inputLocus to this SVLocusSet
  ///
  /// Copies the inputLocus into this object without attempting to do any merging. This is an intermediate
  /// (private) step in the process of merging the \p inputLocus into this graph.
  ///
  /// \return The locus index assigned to the copy of inputLocus inserted into this object
  LocusIndexType insertLocus(const SVLocus& inputLocus);

  /// \brief Erase the node at \p nodeAddress by calling SVLocus::eraseNode() on the appropriate locus.
  void eraseNode(const NodeAddressType nodeAddress)
  {
    assert(_isIndexed);

    LocusSetIndexerType::iterator iter(_inodes.data().find(nodeAddress));
    if (iter == _inodes.data().end()) return;

    SVLocus& locus(getLocus(nodeAddress.first));
    locus.eraseNode(nodeAddress.second, this);
  }

  /// Test if the addressed node has no edges with sufficient evidence to be retained when
  /// noisy graph elements are removed.
  ///
  /// \return True if the address points to a noise node.
  bool isNoiseNode(const NodeAddressType nodeAddress) const
  {
    return getLocus(nodeAddress.first).isNoiseNode(getMinMergeEdgeCount(), nodeAddress.second);
  }

  /// node with a self edge only:
  bool isSingletonNode(const NodeAddressType inputAddy) const
  {
    const SVLocusNode& inputNode(getNode(inputAddy));
    return ((inputNode.size() == 1) && inputNode.isEdge(inputAddy.second));
  }

  bool isOverlapAllowed() const { return (!_isFinalized); }

  void mergeNodePtr(NodeAddressType fromPtr, NodeAddressType toPtr);

  /// update index when nodes are moved:
  void recieve_flyweight_notification(const SVLocusNodeMoveMessage& msg)
  {
    assert(_isIndexed);

    if (msg.first) {
      // add
#ifdef DEBUG_SVL
      log_os << "SVLocusSetObserver: Adding node: " << msg.second.first << ":" << msg.second.second << "\n";
#endif
      _inodes.data().insert(msg.second);
      updateMaxRegionSize(getNode(msg.second).getInterval());
    } else {
      // delete
#ifdef DEBUG_SVL
      log_os << "SVLocusSetObserver: Deleting node: " << msg.second.first << ":" << msg.second.second << "\n";
#endif
      _inodes.data().erase(msg.second);
    }
  }

  void updateMaxRegionSize(const GenomeInterval& interval)
  {
    assert(interval.tid >= 0);
    const unsigned tid(interval.tid);
    if (tid >= _maxRegionSize.size()) {
      _maxRegionSize.resize((tid + 1), 0);
    }
    _maxRegionSize[tid] = std::max(_maxRegionSize[tid], interval.range.size());
  }

  void reconstructIndex();

  void clearIndex()
  {
    _emptyLoci.clear();
    _inodes.data().clear();
    _maxRegionSize.clear();
  }

#ifdef DEBUG_SVL
  void dumpIndex(std::ostream& os) const;
#endif

  /// \brief Throw an exception if any nodes are overlapping
  ///
  /// \param[in] isFilterNoise If true, consider only signal nodes
  void checkForOverlapNodes(const bool isFilterNoise) const;

  /// \brief Get all non-noise nodes intersecting the node at \p targetNodeAddress.
  ///
  /// This is similar to getIntersectingNodeAddresses, except that we only get the addresses of non-noise
  /// (ie. signal) nodes. Also note that in this method the output address set is appended instead of
  /// replaced.
  ///
  /// A second operation tacked onto this method is that noise nodes intersecting the target node are checked
  /// for intersection to intersectingNoiseNodeTestTargets, if intersection is found, then
  /// isIntersectingNoiseNodeOverlapTestTargets is set true.
  ///
  /// \param[in] filterLocusIndex Exclude all nodes in this locus from the intersecting node set.
  /// \param[in] targetNodeAddress Address of the target node which is the target of the intersection search.
  /// \param[out] intersectingSignalNodeAddresses Append to this set the addresses of all signal nodes found
  /// which intersect the target node.
  /// \param[in] intersectingNoiseNodeTestTargets The node addresses against which all intersecting noise
  /// nodes will be tested for overlap.
  /// \param[out] isIntersectingNoiseNodeOverlapTestTargets True if at least one overlapping noise node
  /// is found to overlap at least one test target.
  void getIntersectingSignalNodeAddresses(
      const LocusIndexType             filterLocusIndex,
      const NodeAddressType            targetNodeAddress,
      std::set<NodeAddressType>&       intersectingSignalNodeAddresses,
      const std::set<NodeAddressType>& intersectingNoiseNodeTestTargets,
      bool&                            isIntersectingNoiseNodeOverlapTestTargets) const;

  ///////////////////// data

  /// This is used to evaluate peak SV evidence density among the overlapping nodes of a set
  /// of intersecting edges.
  struct MergeRegionSumData {
    void clear()
    {
      localNodeOutgoingEdgeEvidence.clear();
      localNodeIncomingEdgeEvidence.clear();
      remoteNodeOutgoingEdgeEvidence.clear();
      remoteNodeIncomingEdgeEvidence.clear();
    }

    using rsum_t = RegionSum<unsigned>;
    rsum_t localNodeOutgoingEdgeEvidence;
    rsum_t localNodeIncomingEdgeEvidence;
    rsum_t remoteNodeOutgoingEdgeEvidence;
    rsum_t remoteNodeIncomingEdgeEvidence;
  };

  /// This object keeps its own copy of bam header info (1) for convenience (2) as an QC mechanism to ensure
  /// it is not scored against a mismatched reference (3) as an auditing mechanism so the serialized object is
  /// better documented.
  bam_header_info _bamHeaderInfo;

  SVLocusSetOptions _opt;

  /// Contains the full set of loci in this graph.
  locusset_type _loci;

  /// \brief The indexes of loci which are empty.
  ///
  /// This data is useful because empty loci are never removed from the graph, otherwise many locus index
  /// numbers (which is set for each locus to the locus's index in the _loci vector) would have to be
  /// frequenty updated. Instead of removing loci from the vector, the positions are marked as empty and
  /// prioritized for use the next time a new locus needs to be stored in the graph.
  std::set<unsigned> _emptyLoci;

  // provides an intersection search of overlapping nodes given a bound node size:
  LocusSetIndexerType _inodes;

  /// \brief The largest node breakend region in this graph for each chromosome.
  ///
  /// This is used to support the graph's region-based query scheme to find all nodes overlapping a given a
  /// certain interval from a chromosome. Instead of using a fully general indexing scheme, such as an
  /// interval tree, the scheme leverages the fact that the largest breakend region stored from each
  /// chromosome is relatively small, by supplementing a simple binary search with a linear walk of the size
  /// provided by this object, which guarantees that all intersecting nodes are found using a relatively
  /// lightweight procedure.
  std::vector<unsigned> _maxRegionSize;

  /// simple debug string describing the source of this
  std::string _source;

  /// the graph has intermediate states (during build) when overlapping regions are allowed,
  /// once complete, overlaps are not present and disallowed:
  bool _isFinalized;

  AllSampleReadCounts _counts;

  /// Total number of observations removed on edges with less than minMergeEdgeCount counts
  unsigned _totalCleaned;

  mutable unsigned _highestSearchCount;    ///< highest search count observed during graph build
  mutable float    _highestSearchDensity;  ///< highest node density observed during graph build

  mutable bool _isMaxSearchCount;    ///< has input been filtered because we hit the maximum search count
  mutable bool _isMaxSearchDensity;  ///< has input been filtered because we hit the maximum node density

  /// True if indexing is setup to support region-based queries of graph nodes. Such queries are required for
  /// merging but not used for variant calling.
  bool _isIndexed;

  CpuTimes _buildTime;
  CpuTimes _mergeTime;

  /// \brief A temporary data structure used by the RegionCheck process.
  ///
  /// The RegionCheck process searches for peak SV evidence density among a set of overlapping nodes.
  mutable MergeRegionSumData _mergeRegions;
};

std::ostream& operator<<(std::ostream& os, const SVLocusSet::NodeAddressType& a);
