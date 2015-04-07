// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2015 Illumina, Inc.
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
operator<<(std::ostream& os, const SVSequenceFragmentAssociation& sva)
{
    os << " svindex: " << sva.index << " evidenceType: " << SVEvidenceType::label(sva.evtype);
    return os;
}



std::ostream&
operator<<(std::ostream& os, const SVCandidateSetSequenceFragment& svp)
{
    os << "SVCandidateReadPair svIndices:";
    for (const SVSequenceFragmentAssociation& sva : svp.svLink)
    {
        os << sva << "\n";
    }
    os << "\n";
    os << "\tread1: " << svp.read1;
    os << "\tread2: " << svp.read2;
    return os;
}



SVCandidateSetSequenceFragment*
SVCandidateSetSequenceFragmentSampleGroup::
getSequenceFragment(
    const pindex_t::key_type& key)
{
    const pindex_t::const_iterator kiter(_pairIndex.find(key));

    if (kiter == _pairIndex.end())
    {
        /// don't add more pairs to the object once it's full:
        if (isFull()) return nullptr;

        _pairIndex[key] = _pairs.size();
        _pairs.emplace_back();
        return &(_pairs.back());
    }
    else
    {
        return &(_pairs[kiter->second]);
    }
}



void
SVCandidateSetSequenceFragmentSampleGroup::
add(const bam_record& bamRead,
    const bool isExpectRepeat,
    const bool isNode1,
    const bool isSubMapped)
{
    using namespace illumina::common;

#ifdef DEBUG_SVDATA
    log_os << "SVDataGroup adding: " << bamRead << "\n";
#endif

    SVCandidateSetSequenceFragment* fragPtr(getSequenceFragment(bamRead.qname()));
    if (nullptr == fragPtr) return;

    SVCandidateSetSequenceFragment& fragment(*fragPtr);

    SVCandidateSetRead* targetReadPtr(nullptr);
    if (2 == bamRead.read_no())
    {
        if (bamRead.is_supplement())
        {
            fragment.read2Supplemental.emplace_back();
            targetReadPtr = (&(fragment.read2Supplemental.back()));
        }
        else
        {
            targetReadPtr = (&(fragment.read2));
        }
    }
    else
    {
        if (bamRead.is_supplement())
        {
            fragment.read1Supplemental.emplace_back();
            targetReadPtr = (&(fragment.read1Supplemental.back()));
        }
        else
        {
            targetReadPtr = (&(fragment.read1));
        }
    }

    SVCandidateSetRead& targetRead(*targetReadPtr);
    if (targetRead.isSet())
    {
        if (isExpectRepeat) return;

        std::ostringstream oss;
        oss << "Unexpected read name collision.\n"
            << "\tExisting read: " << targetRead << "\n"
            << "\tNew read: " << bamRead << "\n";
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }


    targetRead.bamrec = bamRead;
    targetRead.isNode1 = isNode1;
    targetRead.isSubMapped = isSubMapped;
    targetRead.mappedReadCount = _mappedReadCount;
    targetRead.subMappedReadCount = _subMappedReadCount;
}



bool
SVCandidateSetData::
setNewSearchInterval(const GenomeInterval& newSearch)
{
    bool retval(false);
    for (const GenomeInterval& oldSearch : _searchIntervals)
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
