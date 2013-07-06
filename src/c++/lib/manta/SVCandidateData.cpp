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

#include "blt_util/log.hh"
#include "common/Exceptions.hh"
#include "manta/SVCandidateData.hh"

#include <cassert>

#include <iostream>
#include <sstream>



std::ostream&
operator<<(std::ostream& os, const bam_record& br)
{
    os << br.qname() << "/" << br.read_no() << " tid:pos " << br.target_id() << ":" << (br.pos()-1);
    return os;
}



std::ostream&
operator<<(std::ostream& os, const SVCandidateRead& svr)
{
    os << "SVCandidateRead: " << svr.bamrec << "\n";
    return os;
}



std::ostream&
operator<<(std::ostream& os, const SVCandidateReadPair& svp)
{
    os << "SVCandidateReadPair svIndex: " << svp.svIndex << "\n";
    os << "\tread1: " << svp.read1;
    os << "\tread2: " << svp.read2;
    return os;
}



SVCandidateReadPair&
SVCandidateDataGroup::
getReadPair(const pindex_t::key_type& key)
{
    const pindex_t::const_iterator kiter(_pairIndex.find(key));

    if(kiter == _pairIndex.end())
    {
        _pairIndex[key] = _pairs.size();
        _pairs.push_back(SVCandidateReadPair());
        return _pairs.back();
    }
    else
    {
        return _pairs[kiter->second];
    }
}



void
SVCandidateDataGroup::
add(const bam_record& bamRead)
{
    using namespace illumina::common;

#ifdef DEBUG_SVDATA
    log_os << "SVDataGroup adding: " << bamRead << "\n";
#endif

    SVCandidateReadPair& pair(getReadPair(bamRead.qname()));

    SVCandidateRead* targetReadPtr(&(pair.read1));
    if(2 == bamRead.read_no())
    {
        targetReadPtr = (&(pair.read2));
    }
    if(targetReadPtr->isSet())
    {
        std::ostringstream oss;
        oss << "Unexpected read name collision.\n"
            << "\tExisting read: " << (*targetReadPtr) << "\n"
            << "\tNew read: " << bamRead << "\n";
        BOOST_THROW_EXCEPTION(PreConditionException(oss.str())); 
    }
    targetReadPtr->bamrec = bamRead;
}
