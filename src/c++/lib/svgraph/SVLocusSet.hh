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

#pragma once

#include "blt_util/bam_header_info.hh"
#include "svgraph/SVLocus.hh"

#include <iosfwd>
#include <string>
#include <vector>



// A set of non-overlapping SVLocus objects
//
struct SVLocusSet : public observer<SVLocusNodeMoveMessage>
{
    typedef std::vector<SVLocus> locusset_type;
    typedef locusset_type::const_iterator const_iterator;

    SVLocusSet(
        const unsigned minMergeEdgeCount = 2) :
        _inodes(NodeAddressSorter(*this)),
        _source("UNKNOWN"),
        _minMergeEdgeCount(minMergeEdgeCount),
        _isFinalized(false),
        _totalCleaned(0),
        _totalAnom(0),
        _totalNonAnom(0)
    {}

    bool
    empty() const
    {
        return _loci.empty();
    }

    unsigned
    size() const
    {
        return _loci.size();
    }

    unsigned
    nonEmptySize() const
    {
        return size()-_emptyLoci.size();
    }

    const_iterator
    begin() const
    {
        return _loci.begin();
    }

    const_iterator
    end() const
    {
        return _loci.end();
    }

    const SVLocus&
    getLocus(const LocusIndexType index) const
    {
#ifdef DEBUG_SVL
        if (index>=_loci.size()) locusHurl(index,"const");
#endif

        assert(index<_loci.size());
        return _loci[index];
    }

    /// merge locus into this:
    ///
    /// locus is destroyed in this process
    ///
    void
    merge(const SVLocus& locus);

    /// merge locus set into this:
    ///
    /// locus set is destroyed in this process
    ///
    void
    merge(const SVLocusSet& set);

    /// indicates the total count of non-filtered
    /// anomolous reads used to construct the graph
    ///
    /// useful for tracking read stats/chimera rates
    void
    addAnomCount(const unsigned count = 1)
    {
        _totalAnom += count;
    }

    /// indicates the total count of non-filtered
    /// non-anomolous reads used to construct the graph
    ///
    /// useful for tracking read stats/chimera rates
    void
    addNonAnomCount(const unsigned count = 1)
    {
        _totalNonAnom += count;
    }

    void
    clear()
    {
        _loci.clear();
        clearIndex();
        _isFinalized=false;
        _totalCleaned=0;
        _totalAnom=0;
        _totalNonAnom=0;
    }

    /// indicate that the set is complete
    void
    finalize()
    {
        clean();
        _isFinalized=true;
    }

    /// remove all existing edges with less than minMergeEdgeCount support:
    void
    clean();

    void
    cleanRegion(const GenomeInterval interval);

    unsigned
    totalCleaned() const
    {
        return _totalCleaned;
    }

    // binary serialization
    void
    save(const char* filename) const;

    // restore from binary serialization
    void
    load(const char* filename);

    // debug output
    void
    dump(std::ostream& os) const;

    // debug output
    void
    dumpRegion(
        std::ostream& os,
        const GenomeInterval interval);

    // dump stats on the whole SVLocus set:
    void
    dumpStats(std::ostream& os) const;

    // dump stats on each locus in tsv format:
    void
    dumpLocusStats(std::ostream& os) const;
    const std::string&
    getSource() const
    {
        return _source;
    }

    unsigned
    getMinMergeEdgeCount() const
    {
        return _minMergeEdgeCount;
    }

    // total number of reads used as supporting evidence in the graph
    unsigned
    totalObservationCount() const
    {
        unsigned sum(0);
        BOOST_FOREACH(const SVLocus& locus, *this)
        {
            sum += locus.totalObservationCount();
        }
        return sum;
    }

    // total nodes in the graph
    unsigned
    totalNodeCount() const
    {
        unsigned sum(0);
        BOOST_FOREACH(const SVLocus& locus, *this)
        {
            sum += locus.size();
        }
        return sum;
    }

    // total number of directed edges in the graph
    unsigned
    totalEdgeCount() const
    {
        unsigned sum(0);
        BOOST_FOREACH(const SVLocus& locus, *this)
        {
            sum += locus.totalEdgeCount();
        }
        return sum;
    }

    /// check that internal data-structures are in
    /// a consistent state, throw on error
    void
    checkState(
        const bool isCheckOverlap = false,
        const bool isCheckLocusConnected = false) const;

private:

    typedef std::pair<LocusIndexType,NodeIndexType> NodeAddressType;

    typedef NodeAddressType EdgeMapKeyType;
    typedef NodeIndexType EdgeMapValueType;

    typedef std::multimap<EdgeMapKeyType, EdgeMapValueType> EdgeMapType;

    typedef std::pair<EdgeMapKeyType, EdgeMapValueType> EdgeInfoType;

    struct NodeAddressSorter
    {

        NodeAddressSorter(const SVLocusSet& set) :
            _set(set)
        {}

        bool
        operator()(
            const NodeAddressType& a,
            const NodeAddressType& b) const
        {
            if (getInterval(a)<getInterval(b)) return true;
            if (getInterval(a)==getInterval(b))
            {
                return (a<b);
            }
            return false;
        }

    private:
        const GenomeInterval&
        getInterval(const NodeAddressType& n) const
        {
            return (_set.getLocus(n.first).getNode(n.second).interval);
        }

        const SVLocusSet& _set;
    };

    typedef std::set<NodeAddressType, NodeAddressSorter> LocusSetIndexerType;

    friend
    std::ostream&
    operator<<(std::ostream& os, const NodeAddressType& a);

    SVLocus&
    getLocus(const LocusIndexType index)
    {
#ifdef DEBUG_SVL
        if (index>=_loci.size()) locusHurl(index,"non-const");
#endif

        assert(index<_loci.size());
        return _loci[index];
    }

    void
    locusHurl(const LocusIndexType index, const char* label) const;


    const SVLocusNode&
    getNode(const NodeAddressType n) const
    {
        return getLocus(n.first).getNode(n.second);
    }

    void
    clearLocus(const LocusIndexType locusIndex)
    {
        _loci[locusIndex].clear();
        _emptyLoci.insert(locusIndex);
        _source="UNKNOWN";
    }

    /// shared node intersection utility
    ///
    /// \param searchNodes the set of nodes to search for intersections in
    /// \param filterLocusIndex ignore intersections from this locus
    ///
    void
    getNodeIntersectCore(
        const LocusIndexType inputLocusIndex,
        const NodeIndexType inputNodeIndex,
        const LocusSetIndexerType& searchNodes,
        const LocusIndexType fitlerLocusIndex,
        std::set<NodeAddressType>& intersectNodes) const;

    /// get all nodes in this object which intersect with
    /// the inputNode
    void
    getNodeIntersect(
        const LocusIndexType locusIndex,
        const NodeIndexType nodeIndex,
        std::set<NodeAddressType>& intersectNodes) const
    {
        getNodeIntersectCore(locusIndex, nodeIndex, _inodes, locusIndex, intersectNodes);
    }

    void
    getIntersectingEdgeNodes(
        const LocusIndexType inputLocusIndex,
        const NodeIndexType inputRemoteNodeIndex,
        const EdgeMapType& remoteToLocal,
        const LocusSetIndexerType& remoteIntersect,
        std::vector<EdgeInfoType>& edges) const;

    ///
    /// \param isInputLocusMoved has the input locus been moved into the graph from an initial temporary locus?
    ///
    void
    getNodeMergeableIntersect(
        const LocusIndexType inputLocusIndex,
        const NodeIndexType inputNodeIndex,
        const bool isInputLocusMoved,
        std::set<NodeAddressType>& mergeIntersect) const;

    /// get all nodes in this object which intersect with
    /// a external node
    void
    getRegionIntersect(
        const GenomeInterval interval,
        std::set<NodeAddressType>& intersectNodes);

    /// assign all intersect clusters to the lowest index number that is not startLocusIndex
    ///
    void
    moveIntersectToLowIndex(
        const std::set<NodeAddressType>& intersectNodes,
        const LocusIndexType startLocusIndex,
        LocusIndexType& locusIndex);

    /// combine all content from 'from' locus into 'to' locus
    ///
    /// this is typically required when a node is merged
    /// which combines two loci
    void
    combineLoci(
        const LocusIndexType fromIndex,
        const LocusIndexType toIndex,
        const bool isClearSource = true);


    // add locus to this locusSet (intermediate step in merging)
    LocusIndexType
    insertLocus(
        const SVLocus& inputLocus);

    void
    removeNode(const NodeAddressType inputNodePtr)
    {
        LocusSetIndexerType::iterator iter(_inodes.find(inputNodePtr));
        if (iter == _inodes.end()) return;

        SVLocus& locus(getLocus(inputNodePtr.first));
        locus.eraseNode(inputNodePtr.second);
    }

    bool
    isOverlapAllowed() const
    {
        return (! _isFinalized);
    }

    void
    mergeNodePtr(NodeAddressType fromPtr,
                 NodeAddressType toPtr);

    // update index when nodes are moved:
    void
    recieve_notification(const notifier<SVLocusNodeMoveMessage>&,
                         const SVLocusNodeMoveMessage& msg)
    {

        if (msg.first)
        {
            // add
#ifdef DEBUG_SVL
            log_os << "SVLocusSetObserver: Adding node: " << msg.second.first << ":" << msg.second.second << "\n";
#endif
            _inodes.insert(msg.second);
            updateMaxRegionSize(getNode(msg.second).interval);
        }
        else
        {
            // delete
#ifdef DEBUG_SVL
            log_os << "SVLocusSetObserver: Deleting node: " << msg.second.first << ":" << msg.second.second << "\n";
#endif
            _inodes.erase(msg.second);
        }
    }

    void
    updateMaxRegionSize(const GenomeInterval& interval)
    {
        const unsigned tid(interval.tid);
        if (tid >= _maxRegionSize.size())
        {
            _maxRegionSize.resize((tid+1),0);
        }
        _maxRegionSize[tid] = std::max(_maxRegionSize[tid], interval.range.size());
    }


    void
    reconstructIndex();

    void
    clearIndex()
    {
        _emptyLoci.clear();
        _inodes.clear();
        _maxRegionSize.clear();
    }

    void
    dumpIndex(std::ostream& os) const;

    ///////////////////// data

public:
    bam_header_info header;
private:

    // contains the full set of loci
    locusset_type _loci;
    std::set<unsigned> _emptyLoci;

    // provides an intersection search of overlapping nodes given a bound node size:
    LocusSetIndexerType _inodes;

    // maximum region size per chromosome:
    std::vector<unsigned> _maxRegionSize;

    // simple debug string describing the source of this
    std::string _source;

    // to reduce noise in the graph, we only merge once shared edges (including self-edges), reach this count:
    unsigned _minMergeEdgeCount;

    // the graph has intermediate states (during build) when overlapping regions are allowed,
    // once complete, overlaps are not present and disallowed:
    bool _isFinalized;

    // total number of observations removed on edges with less than minMergeEdgeCount counts
    unsigned _totalCleaned;

    // total number of non-filtered anomalous reads scanned but not in graph:
    unsigned long _totalAnom;

    // total number of non-filtered non-anomalous reads scanned but not in graph:
    unsigned long _totalNonAnom;
};


std::ostream&
operator<<(std::ostream& os, const SVLocusSet::NodeAddressType& a);
