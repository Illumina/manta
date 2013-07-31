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


struct SVCandidateRead
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
operator<<(std::ostream& os, const SVCandidateRead& svr);


struct SVCandidateReadPair
{
    SVCandidateReadPair()
        : svIndex(0)
    {}

    // which sv is this read associated with?
    unsigned short svIndex;
    SVCandidateRead read1;
    SVCandidateRead read2;
};


std::ostream&
operator<<(std::ostream& os, const SVCandidateReadPair& svp);



/// SVCandidateData associated with a specific bam-file/read-group
struct SVCandidateReadPairGroup
{
    typedef std::vector<SVCandidateReadPair> pair_t;
    typedef pair_t::iterator iterator;
    typedef pair_t::const_iterator const_iterator;

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

    SVCandidateReadPair&
    getReadPair(const pindex_t::key_type& key);

    pair_t _pairs;
    pindex_t _pairIndex;
};



/// any information which travels with the final SVCandidate
/// which would be useful at score time
///
/// example: realigned/reassembled reads
///
/// ideally scoring model should not have to touch bam again after hygen step..
///

struct SVCandidateData
{
	// should be a template...
	typedef JumpAlignmentResult<int> JumpAlignmentResultType;

	SVCandidateReadPairGroup&
    getDataGroup(const unsigned bamIndex)
    {
        data_t::iterator diter(_data.find(bamIndex));
        if (diter != _data.end()) return diter->second;

        std::pair<data_t::iterator,bool> diter2 = _data.insert(std::make_pair(bamIndex,SVCandidateReadPairGroup()));
        return diter2.first->second;
    }

    const SVCandidateReadPairGroup&
    getDataGroup(const unsigned bamIndex) const
    {
        data_t::const_iterator diter(_data.find(bamIndex));
        assert(diter != _data.end());
        return diter->second;
    }

    Assembly&
    getAssembly() { return ctgs; }

    const Assembly&
    getAssembly() const { return ctgs; }

    std::vector<JumpAlignmentResultType>&
    getAlignments() { return align; }

    const std::vector<JumpAlignmentResultType>&
    getAlignments() const { return align; }


    void
    clear()
    {
        _data.clear();
    }

    /// return true if this search interval overlaps with any previous:
    bool
    setNewSearchInterval(const GenomeInterval& newSearch);

private:
    typedef std::map<unsigned,SVCandidateReadPairGroup> data_t;
    data_t _data;

    std::vector<GenomeInterval> _searchIntervals;

    // assembled contigs for both breakpoints
    Assembly ctgs;
    // contig alignments
    std::vector<JumpAlignmentResultType > align;

};
