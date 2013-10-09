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

//#include "blt_util/log.hh"
#include "common/Exceptions.hh"
#include "manta/SVCandidateSetData.hh"

#include "boost/foreach.hpp"

#include <cassert>

#include <iostream>
#include <sstream>



std::ostream&
operator<<(std::ostream& os, const SVCandidateSetRead& svr)
{
    os << "SVCandidateSetRead: " << svr.bamrec << "\n";
    return os;
}



std::ostream&
operator<<(std::ostream& os, const SVCandidateSetReadPair& svp)
{
    os << "SVCandidateReadPair svIndices:";
    BOOST_FOREACH(const SVCandidateSetReadPair::index_t index, svp.svIndex)
    {
        os << " " << index;
    }
    os << "\n";
    os << "\tread1: " << svp.read1;
    os << "\tread2: " << svp.read2;
    return os;
}



SVCandidateSetReadPair&
SVCandidateSetReadPairSampleGroup::
getReadPair(const pindex_t::key_type& key)
{
    const pindex_t::const_iterator kiter(_pairIndex.find(key));

    if (kiter == _pairIndex.end())
    {
        _pairIndex[key] = _pairs.size();
        _pairs.push_back(SVCandidateSetReadPair());
        return _pairs.back();
    }
    else
    {
        return _pairs[kiter->second];
    }
}



void
SVCandidateSetReadPairSampleGroup::
add(const bam_record& bamRead,
    const bool isExpectRepeat)
{
    using namespace illumina::common;

#ifdef DEBUG_SVDATA
    log_os << "SVDataGroup adding: " << bamRead << "\n";
#endif

    SVCandidateSetReadPair& pair(getReadPair(bamRead.qname()));

    SVCandidateSetRead* targetReadPtr(&(pair.read1));
    if (2 == bamRead.read_no())
    {
        targetReadPtr = (&(pair.read2));
    }
    if (targetReadPtr->isSet())
    {
        //if (isExpectRepeat)
        if (isExpectRepeat || bamRead.is_supplement())
        {
            return;
        }
        /*
        if (bamRead.is_supplement()) {
          oss << "Found Supp: " << bamRead << std::endl;
          BOOST_THROW_EXCEPTION(LogicException(oss.str()));
        }
        */

        //if (! bamRead.is_supplement())
        //{
        std::ostringstream oss;
            oss << "Unexpected read name collision.\n"
                << "\tExisting read: " << (*targetReadPtr) << "\n"
                << "\tNew read: " << bamRead << "\n";
            //oss << "\tIs Supp 1: " << targetReadPtr->is_supplement() << std::endl;
            oss << "\tIs Supp 2: " << bamRead.is_supplement() << std::endl;
            BOOST_THROW_EXCEPTION(LogicException(oss.str()));
            //}
    }
    if (! bamRead.is_supplement())
    {
      targetReadPtr->bamrec = bamRead;
    }
}


bool
SVCandidateSetData::
setNewSearchInterval(const GenomeInterval& newSearch)
{
    bool retval(false);
    BOOST_FOREACH(const GenomeInterval& oldSearch, _searchIntervals)
    {
        if (oldSearch.isIntersect(newSearch))
        {
            retval=true;
            break;
        }
    }
    _searchIntervals.push_back(newSearch);
    return retval;
}

