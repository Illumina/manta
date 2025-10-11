//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2025 Illumina, Inc.
//
// This program is licensed under the terms of the Polyform strict license
//
// ***As far as the law allows, the software comes as is, without
// any warranty or condition, and the licensor will not be liable
// to you for any damages arising out of these terms or the use
// or nature of the software, under any kind of legal claim.***
//
// You should have received a copy of the PolyForm Strict License 1.0.0
// along with this program.  If not, see <https://polyformproject.org/licenses/strict/1.0.0>.
//
//

/// \file
/// \author Chris Saunders

#include "manta/SVCandidateSetData.hpp"
#include "common/Exceptions.hpp"
#include "htsapi/bam_record_util.hpp"

#include <cassert>

#include <iostream>
#include <sstream>

#ifdef DEBUG_SVDATA
#include "blt_util/log.hpp"
#endif

std::ostream& operator<<(std::ostream& os, const SVCandidateSetRead& svr)
{
  os << "SVCandidateSetRead: " << svr.bamrec << "\n";
  os << "Read_index: " << svr.readIndex << "\n";
  return os;
}

std::ostream& operator<<(std::ostream& os, const SVSequenceFragmentAssociation& sva)
{
  os << " svindex: " << sva.index << " evidenceType: " << SVEvidenceType::label(sva.evtype);
  return os;
}

std::ostream& operator<<(std::ostream& os, const SVCandidateSetSequenceFragment& svp)
{
  os << "SVCandidateReadPair svIndices:";
  for (const SVSequenceFragmentAssociation& sva : svp.svLink) {
    os << sva << "\n";
  }
  os << "\n";
  os << "\tread1: " << svp.read1;
  os << "\tread2: " << svp.read2;
  return os;
}

SVCandidateSetSequenceFragment* SVCandidateSetSequenceFragmentSampleGroup::getSequenceFragment(
    const pindex_t::key_type& key)
{
  const pindex_t::const_iterator kiter(_pairIndex.find(key));

  if (kiter == _pairIndex.end()) {
    /// don't add more pairs to the object once it's full:
    if (isFull()) return nullptr;

    _pairIndex[key] = _pairs.size();
    _pairs.emplace_back();
    return &(_pairs.back());
  } else {
    return &(_pairs[kiter->second]);
  }
}

void SVCandidateSetSequenceFragmentSampleGroup::add(
    const bam_header_info& bamHeader,
    const bam_record&      bamRead,
    const bool             isExpectRepeat,
    const bool             isSourcedFromGraphEdgeNode1,
    const bool             isSubMapped)
{
  using namespace illumina::common;

#ifdef DEBUG_SVDATA
  log_os << "SVDataGroup adding: " << bamRead << "\n";
#endif

  SVCandidateSetSequenceFragment* fragPtr(getSequenceFragment(bamRead.qname()));
  if (nullptr == fragPtr) return;

  SVCandidateSetSequenceFragment& fragment(*fragPtr);

  SVCandidateSetRead* targetReadPtr(nullptr);
  if (2 == bamRead.read_no()) {
    if (bamRead.isNonStrictSupplement()) {
      fragment.read2Supplemental.emplace_back();
      targetReadPtr = (&(fragment.read2Supplemental.back()));
    } else {
      targetReadPtr = (&(fragment.read2));
    }
  } else {
    if (bamRead.isNonStrictSupplement()) {
      fragment.read1Supplemental.emplace_back();
      targetReadPtr = (&(fragment.read1Supplemental.back()));
    } else {
      targetReadPtr = (&(fragment.read1));
    }
  }

  SVCandidateSetRead& targetRead(*targetReadPtr);
  if (targetRead.isSet()) {
    if (isExpectRepeat) return;

    std::ostringstream oss;
    oss << "Unexpected alignment name collision. Source: '" << dataSourceName << "'\n"
        << "\tExisting read: ";
    summarizeAlignmentRecord(bamHeader, targetRead.bamrec, oss);
    oss << "\n"
        << "\tNew read: ";
    summarizeAlignmentRecord(bamHeader, bamRead, oss);
    oss << "\n";
    BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
  }

  targetRead.bamrec                      = bamRead;
  targetRead.isSourcedFromGraphEdgeNode1 = isSourcedFromGraphEdgeNode1;
  targetRead.isSubMapped                 = isSubMapped;
  targetRead.readIndex                   = (isSubMapped ? _subMappedReadIndex : _mappedReadIndex);
}

bool SVCandidateSetData::setNewSearchInterval(const GenomeInterval& newSearch)
{
  bool retval(false);
  for (const GenomeInterval& oldSearch : _searchIntervals) {
    if (oldSearch.isIntersect(newSearch)) {
      retval = true;
      break;
    }
  }
  _searchIntervals.push_back(newSearch);
  return retval;
}
