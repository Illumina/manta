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
#include "manta/SVLocus.hh"

#include <iosfwd>
#include <string>
#include <vector>



// A set of non-overlapping SVLocus objects
//
struct SVLocusSet : public observer<SVLocusNodeMoveMessage>
{
    typedef std::vector<SVLocus> locusset_type;
    typedef locusset_type::const_iterator const_iterator;

    SVLocusSet() :
        _inodes(NodeAddressSorter(*this)),
        _source("UNKNOWN")
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
        assert(index<_loci.size());
        return _loci[index];
    }

    /// merge locus into this:
    ///
    /// locus is destroyed in this process
    ///
    void
    merge(SVLocus& locus);

    /// merge locus set into this:
    ///
    /// locus set is destroyed in this process
    ///
    void
    merge(SVLocusSet& set);

    void
    clear()
    {
        _loci.clear();
        _emptyLoci.clear();
        _inodes.clear();
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
        const int32_t tid,
        const int32_t beginPos,
        const int32_t endPos);

    // dump stats on each locus in tsv format:
    void
    dumpStats(std::ostream& os) const;

    const std::string&
    getSource() const
    {
        return _source;
    }


private:

    typedef std::pair<LocusIndexType,NodeIndexType> NodeAddressType;

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
        assert(index<_loci.size());
        return _loci[index];
    }

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

    /// get all nodes in this object which intersect with
    /// the inputNode
    void
    getNodeIntersect(
        const LocusIndexType locusIndex,
        const NodeIndexType nodeIndex,
        LocusSetIndexerType& intersect);

    /// get all nodes in this object which intersect with
    /// a external node
    void
    getRegionIntersect(
        const int32_t tid,
        const int32_t beginPos,
        const int32_t endPos,
        LocusSetIndexerType& intersect);

    /// assign all intersect clusters to the lowest index number that is not startLocusIndex
    ///
    void
    moveIntersectToLowIndex(
        LocusSetIndexerType& intersect,
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
    reconstructIndex();

    void
    dumpIndex(std::ostream& os) const;


    /// check that internal data-structures are in
    /// a consistent state, throw on error
    void
    checkState(const bool isCheckOverlap = false) const;

    ///////////////////// data

public:
    bam_header_info header;
private:
    // contains the full set of loci
    locusset_type _loci;
    std::set<unsigned> _emptyLoci;

    // provides an intersection search of non-overlapping nodes:
    LocusSetIndexerType _inodes;

    // simple debug string describing the source of this
    std::string _source;
};


std::ostream&
operator<<(std::ostream& os, const SVLocusSet::NodeAddressType& a);
