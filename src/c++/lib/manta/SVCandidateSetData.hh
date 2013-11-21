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

#pragma once

#include "blt_util/bam_record.hh"
#include "svgraph/GenomeInterval.hh"
#include "manta/SVBreakend.hh"

#include "alignment/Alignment.hh"

#include <iosfwd>
#include <map>
#include <vector>

//#define DEBUG_SVDATA

#ifdef DEBUG_SVDATA
#include "blt_util/log.hh"
#include <iostream>
#endif


/// A read associated with an SV associated set of regions
///
/// note this read could be linked with zero to many specific SVCandidates
///
struct SVCandidateSetRead
{
    SVCandidateSetRead() :
        isNode1(true)
    {}

    bool
    isSet() const
    {
        return (! bamrec.empty());
    }

    //realignment info, etc...
    bam_record bamrec;
    bool isNode1; ///< used to link this read to node1 or node2 in the original graph ordering, note this is not the same as read1 and read2
};

std::ostream&
operator<<(std::ostream& os, const SVCandidateSetRead& svr);


/// capture details of the link between a read pair and an SV
struct SVPairAssociation
{
    typedef uint16_t index_t;

    explicit
    SVPairAssociation(
        const index_t initIndex = 0,
        const SVEvidenceType::index_t initEvtype = SVEvidenceType::UNKNOWN) :
        index(initIndex),
        evtype(initEvtype)
    {}

    index_t index;
    SVEvidenceType::index_t evtype;  // is this association drawn form a read pair or other evidence source? (ie. CIGAR)
};

std::ostream&
operator<<(std::ostream& os, const SVPairAssociation& sva);


/// A read associated with an SV associated set of regions
///
/// note this read could be linked with zero to many specific SVCandidates
///
struct SVCandidateSetReadPair
{
    SVCandidateSetReadPair()
    {}

    const char*
    qname() const
    {
        if      (read1.isSet()) return read1.bamrec.qname();
        else if (read2.isSet()) return read2.bamrec.qname();
        return NULL;
    }

    std::vector<SVPairAssociation> svLink; ///< which SVs from the set are this molecule associated with?
    SVCandidateSetRead read1;
    SVCandidateSetRead read2;
};

std::ostream&
operator<<(std::ostream& os, const SVCandidateSetReadPair& svp);



/// SVCandidateSet data associated with a specific bam-file/read-group
///
struct SVCandidateSetReadPairSampleGroup
{
    typedef std::vector<SVCandidateSetReadPair> pair_t;
    typedef pair_t::iterator iterator;
    typedef pair_t::const_iterator const_iterator;

    SVCandidateSetReadPairSampleGroup() :
        _isIncomplete(false)
    {}

    /// add a new bam record to the set:
    void
    add(const bam_record& bamRead,
        const bool isExpectRepeat,
        const bool isNode1);

    iterator
    begin()
    {
        return _pairs.begin();
    }

    iterator
    end()
    {
        return _pairs.end();
    }

    const_iterator
    begin() const
    {
        return _pairs.begin();
    }

    const_iterator
    end() const
    {
        return _pairs.end();
    }

    unsigned
    size() const
    {
        return _pairs.size();
    }

    bool
    isIncomplete() const
    {
        return _isIncomplete;
    }

    void
    setIncomplete()
    {
        _isIncomplete = true;
    }

private:
    typedef std::string bamqname_t;
    typedef std::map<bamqname_t,unsigned> pindex_t;

    SVCandidateSetReadPair&
    getReadPair(const pindex_t::key_type& key);

    pair_t _pairs;
    pindex_t _pairIndex;

    bool _isIncomplete; ///< this flag can be set if the object grows too large to insert more data into it
};


/// Data gathered from a set of regions implicated to contain one or more SVs
///
/// Note these data are used for initial hypothesis generation, therefore the
/// reads are potentially associated with zero to many specific SV candidates
/// (although we expect any one read to usually be associated with no more
/// one).
///
///
///
struct SVCandidateSetData
{
    SVCandidateSetData() :
        _isSkipped(false)
    {}

    /// get evidence associated with a specific sample group:
    SVCandidateSetReadPairSampleGroup&
    getDataGroup(const unsigned bamIndex)
    {
        data_t::iterator diter(_data.find(bamIndex));
        if (diter != _data.end()) return diter->second;

        std::pair<data_t::iterator,bool> diter2 = _data.insert(std::make_pair(bamIndex,SVCandidateSetReadPairSampleGroup()));
        return diter2.first->second;
    }

    /// get evidence associated with a specific sample group:
    const SVCandidateSetReadPairSampleGroup&
    getDataGroup(const unsigned bamIndex) const
    {
        data_t::const_iterator diter(_data.find(bamIndex));
        assert(diter != _data.end());
        return diter->second;
    }

    void
    clear()
    {
        _data.clear();
        _searchIntervals.clear();
        _isSkipped=false;
    }

    void
    setSkipped()
    {
        _isSkipped=true;
    }

    bool
    isSkipped() const
    {
        return _isSkipped;
    }

    /// return true if this search interval overlaps with any previous
    /// intervals
    ///
    /// if true, we know to expect repeated qnames
    ///
    bool
    setNewSearchInterval(const GenomeInterval& newSearch);

private:
    typedef std::map<unsigned,SVCandidateSetReadPairSampleGroup> data_t;
    data_t _data;

    std::vector<GenomeInterval> _searchIntervals;

    bool _isSkipped; // bam region was never scanned for this variant (b/c of a temporary state during model development, no fundamemntal reason)
};
