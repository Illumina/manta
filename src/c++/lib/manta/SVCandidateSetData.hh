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

#pragma once

#include "blt_util/bam_record.hh"
#include "svgraph/GenomeInterval.hh"

#include "alignment/Alignment.hh"
#include "alignment/GlobalJumpAligner.hh"

#include "SVLocusAssembler.hh"

#include <iosfwd>
#include <map>
#include <vector>

//#define DEBUG_SVDATA


/// A read associated with an SV associated set of regions
///
/// note this read could be linked with zero to many specific SVCandidates
///
struct SVCandidateSetRead
{
    bool
    isSet() const
    {
        return (! bamrec.empty());
    }

    //realignment info, etc...
    bam_record bamrec;
};

std::ostream&
operator<<(std::ostream& os, const SVCandidateSetRead& svr);



/// A read associated with an SV associated set of regions
///
/// note this read could be linked with zero to many specific SVCandidates
///
struct SVCandidateSetReadPair
{
    SVCandidateSetReadPair()
    {}

    typedef uint16_t index_t;
    std::vector<index_t> svIndex; ///< which SVs from the set are this molecule associated with?
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

    /// add a new bam record to the set:
    void
    add(const bam_record& bamRead,
        const bool isExpectRepeat);

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


private:
    typedef std::string bamqname_t;
    typedef std::map<bamqname_t,unsigned> pindex_t;

    SVCandidateSetReadPair&
    getReadPair(const pindex_t::key_type& key);

    pair_t _pairs;
    pindex_t _pairIndex;
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
    // should be a template...
    typedef JumpAlignmentResult<int> JumpAlignmentResultType;

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

    Assembly&
    getAssembly()
    {
        return ctgs;
    }

    const Assembly&
    getAssembly() const
    {
        return ctgs;
    }

    std::vector<JumpAlignmentResultType>&
    getAlignments()
    {
        return align;
    }

    const std::vector<JumpAlignmentResultType>&
    getAlignments() const
    {
        return align;
    }


    void
    clear()
    {
        _data.clear();
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

    // assembled contigs for both breakpoints
    Assembly ctgs;
    // contig alignments
    std::vector<JumpAlignmentResultType > align;

};
