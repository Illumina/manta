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

#include <map>
#include <vector>


struct SVCandidateRead
{
    bool
    isSet() const
    {
        return (*(bamrec.qname()) != '\0');
    }

    //realignment info, etc...
    bam_record bamrec;
};


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



/// SVCandidateData associated with a specific bam-file/read-group
struct SVCandidateDataGroup
{
    typedef std::vector<SVCandidateReadPair> pair_t;
    typedef pair_t::iterator iterator;
    typedef pair_t::const_iterator const_iterator;

    void
    add(const bam_record& read);

    iterator
    begin()
    {
        return _pairs.begin();
    }

    iterator
    end()
    {
        return _pairs.begin();
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
    SVCandidateDataGroup&
    getDataGroup(const unsigned bamIndex)
    {
        data_t::iterator diter(_data.find(bamIndex));
        if(diter != _data.end()) return diter->second;

        std::pair<data_t::iterator,bool> diter2 = _data.insert(std::make_pair(bamIndex,SVCandidateDataGroup()));
        return diter2.first->second;
    }

    void
    clear()
    {
        _data.clear();
    }

private:
    typedef std::map<unsigned,SVCandidateDataGroup> data_t;
    data_t _data;
};
