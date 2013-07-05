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
    SVCandidateReadPair() :
        bamIndex(0)
    {}

    unsigned bamIndex;
    SVCandidateRead read1;
    SVCandidateRead read2;
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
    void
    add(const bam_record& read,
        const unsigned bamIndex);

private:
    typedef std::pair<unsigned,std::string> bamqname_t;
    typedef std::map<bamqname_t,unsigned> pindex_t;

    SVCandidateReadPair&
    getReadPair(const pindex_t::key_type& key);

    std::vector<SVCandidateReadPair> _pairs;
    pindex_t _pairIndex;
};
