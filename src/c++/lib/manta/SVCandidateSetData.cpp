// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2014 Illumina, Inc.
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
operator<<(std::ostream& os, const SVPairAssociation& sva)
{
    os << " svindex: " << sva.index << " evidenceType: " << SVEvidenceType::label(sva.evtype);
    return os;
}



std::ostream&
operator<<(std::ostream& os, const SVCandidateSetReadPair& svp)
{
    os << "SVCandidateReadPair svIndices:";
    BOOST_FOREACH(const SVPairAssociation& sva, svp.svLink)
    {
        os << sva << "\n";
    }
    os << "\n";
    os << "\tread1: " << svp.read1;
    os << "\tread2: " << svp.read2;
    return os;
}



SVCandidateSetReadPair*
SVCandidateSetReadPairSampleGroup::
getReadPair(const pindex_t::key_type& key)
{
    const pindex_t::const_iterator kiter(_pairIndex.find(key));

    if (kiter == _pairIndex.end())
    {
        /// don't add more pairs to the object once it's full:
        if (isFull()) return NULL;

        _pairIndex[key] = _pairs.size();
        _pairs.push_back(SVCandidateSetReadPair());
        return &(_pairs.back());
    }
    else
    {
        return &(_pairs[kiter->second]);
    }
}



void
SVCandidateSetReadPairSampleGroup::
add(const bam_record& bamRead,
    const bool isExpectRepeat,
    const bool isNode1,
    const bool isSubMapped)
{
    using namespace illumina::common;

#ifdef DEBUG_SVDATA
    log_os << "SVDataGroup adding: " << bamRead << "\n";
#endif

    SVCandidateSetReadPair* pairPtr(getReadPair(bamRead.qname()));
    if (NULL == pairPtr) return;

    SVCandidateSetReadPair& pair(*pairPtr);

    SVCandidateSetRead* targetReadPtr(&(pair.read1));
    if (2 == bamRead.read_no())
    {
        targetReadPtr = (&(pair.read2));
    }
    if (targetReadPtr->isSet())
    {
        if (isExpectRepeat) return;

        std::ostringstream oss;
        oss << "Unexpected read name collision.\n"
            << "\tExisting read: " << (*targetReadPtr) << "\n"
            << "\tNew read: " << bamRead << "\n";
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }

    targetReadPtr->bamrec = bamRead;
    targetReadPtr->isNode1 = isNode1;
    targetReadPtr->isSubMapped = isSubMapped;
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
