//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2019 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

/// \file
/// \author Chris Saunders
///

#pragma once

#include <iosfwd>
#include <map>
#include <vector>

#include "alignment/Alignment.hpp"
#include "htsapi/bam_header_info.hpp"
#include "htsapi/bam_record.hpp"
#include "manta/SVBreakend.hpp"
#include "svgraph/GenomeInterval.hpp"

//#define DEBUG_SVDATA

#ifdef DEBUG_SVDATA
#include <iostream>
#include "blt_util/log.hpp"
#endif

/// \brief A read associated with an SV associated set of regions
///
/// Note this read could be linked with zero to many specific SVCandidates
///
struct SVCandidateSetRead {
  bool isSet() const { return (!bamrec.empty()); }

  bool isAnchored() const { return (isSet() && (!isSubMapped)); }

  // realignment info, etc...
  bam_record bamrec;

  /// True if this read is associated with Node1 in the corresponding SV locus graph node supporting the SV
  /// candidate.
  ///
  /// This is used to link this read to node1 or node2 in the original graph ordering,
  /// note this is not the same as read1 and read2
  bool isSourcedFromGraphEdgeNode1 = true;

  /// is mapq below the minimum normally required to use this read
  bool isSubMapped = false;

  /// relative index of this read compared to all reads with the same mapping status in bam-input order
  double readIndex = 0;
};

std::ostream& operator<<(std::ostream& os, const SVCandidateSetRead& svr);

/// capture details of the link between a sequence fragment and an SV
struct SVSequenceFragmentAssociation {
  typedef uint16_t index_t;

  explicit SVSequenceFragmentAssociation(
      const index_t initIndex = 0, const SVEvidenceType::index_t initEvtype = SVEvidenceType::UNKNOWN)
    : index(initIndex), evtype(initEvtype)
  {
  }

  index_t index;

  /// is the association from anom read pair, split read, CIGAR, etc?
  SVEvidenceType::index_t evtype;
};

std::ostream& operator<<(std::ostream& os, const SVSequenceFragmentAssociation& sva);

/// A DNA/RNA fragment associated with an SV-associated set of regions
///
/// note this read could be linked with zero to many specific SVCandidates
///
struct SVCandidateSetSequenceFragment {
  const char* qname() const
  {
    if (read1.isSet())
      return read1.bamrec.qname();
    else if (read2.isSet())
      return read2.bamrec.qname();
    else if ((!read1Supplemental.empty()) && read1Supplemental.front().isSet()) {
      return read1Supplemental.front().bamrec.qname();
    } else if ((!read2Supplemental.empty()) && read2Supplemental.front().isSet()) {
      return read2Supplemental.front().bamrec.qname();
    }
    return nullptr;
  }

  bool isAnchored() const { return (read1.isAnchored() || read2.isAnchored()); }

  bool checkReadPair() const
  {
    if (read1.isSet() && read2.isSet()) {
      if (read1.bamrec.target_id() != read2.bamrec.mate_target_id()) return false;
      if (read2.bamrec.target_id() != read1.bamrec.mate_target_id()) return false;
      if (read1.bamrec.pos() != read2.bamrec.mate_pos()) return false;
      if (read2.bamrec.pos() != read1.bamrec.mate_pos()) return false;
      if (read1.bamrec.is_fwd_strand() != read2.bamrec.is_mate_fwd_strand()) return false;
      if (read2.bamrec.is_fwd_strand() != read1.bamrec.is_mate_fwd_strand()) return false;
    }

    return true;
  }

  /// Which SVs from the set are this molecule associated with?
  std::vector<SVSequenceFragmentAssociation> svLink;

  SVCandidateSetRead              read1;
  std::vector<SVCandidateSetRead> read1Supplemental;
  SVCandidateSetRead              read2;
  std::vector<SVCandidateSetRead> read2Supplemental;
};

std::ostream& operator<<(std::ostream& os, const SVCandidateSetSequenceFragment& svp);

/// SVCandidateSet data associated with a specific bam-file/read-group
///
struct SVCandidateSetSequenceFragmentSampleGroup {
  typedef std::vector<SVCandidateSetSequenceFragment> pair_t;
  typedef pair_t::iterator                            iterator;
  typedef pair_t::const_iterator                      const_iterator;

  /// \brief Increment once for each eligible read considered
  ///
  /// This information is used to determine signal rate so that it can be compared to sample specific noise
  /// levels
  void increment(const bool isSubMapped)
  {
    if (isSubMapped) {
      _subMappedReadIndex++;
    } else {
      _mappedReadIndex++;
    }
  }

  /// Add a new bam record to the set
  ///
  /// \param[in] bamHeader Bam header information is (only) used to improve the detail of exception messages.
  ///
  /// \param[in] bamRead New bam record to add to the sample group set
  ///
  /// \param[in[ isExpectRepeat If false, raise an exception for a repeated BAM QNAME, otherwise skip all but
  /// the first repeated QNAME instance.
  ///
  /// \param[in] isSourcedFromGraphEdgeNode1 True if \p bamRead was discovered from SV locus graph edge node1
  ///
  /// \param[in] isSubMapped True if read is below default mapping quality threshold
  void add(
      const bam_header_info& bamHeader,
      const bam_record&      bamRead,
      const bool             isExpectRepeat,
      const bool             isSourcedFromGraphEdgeNode1,
      const bool             isSubMapped);

  iterator begin() { return _pairs.begin(); }

  iterator end() { return _pairs.end(); }

  const_iterator begin() const { return _pairs.begin(); }

  const_iterator end() const { return _pairs.end(); }

  unsigned size() const { return _pairs.size(); }

  bool isFull() const { return _isFull; }

  void setFull() { _isFull = true; }

private:
  typedef std::string                    bamqname_t;
  typedef std::map<bamqname_t, unsigned> pindex_t;

  /// get existing fragment or return pointer for a new fragment
  ///
  /// this will return null for new fragments when isFull() is true
  ///
  SVCandidateSetSequenceFragment* getSequenceFragment(const pindex_t::key_type& key);

public:
  /// Record a name for the data source to improve error messages:
  std::string dataSourceName = "UNKNOWN";

private:
  pair_t   _pairs;
  pindex_t _pairIndex;

  bool _isFull = false;  ///< this flag can be set if the object grows too large to insert more data into it

  /// Tracks the relative index of all mapped reads as read off of the input bam file. This is used to provide
  /// a relative index number for all reads supporting a particular SV candidate, so that supporting read
  /// density can be estimated. For instance, if 3 reads supporting a breakpoint have mapped read counts of
  /// {100,200,300}, we can roughly estimate the 1 of 100 reads support the breakpoint.
  double _mappedReadIndex = 0;

  /// same as above for reads with mapping quality below theshold
  double _subMappedReadIndex = 0;
};

/// \brief Data gathered from a set of regions implicated to contain one or more SVs
///
/// Note these data are used for initial hypothesis generation, therefore the
/// reads are potentially associated with zero to many specific SV candidates
/// (although we expect any one read to usually be associated with no more
/// one).
///
struct SVCandidateSetData {
  /// \brief Get evidence associated with a specific sample group:
  SVCandidateSetSequenceFragmentSampleGroup& getDataGroup(const unsigned bamIndex)
  {
    data_t::iterator diter(_data.find(bamIndex));
    if (diter != _data.end()) return diter->second;

    std::pair<data_t::iterator, bool> diter2 =
        _data.insert(std::make_pair(bamIndex, SVCandidateSetSequenceFragmentSampleGroup()));
    return diter2.first->second;
  }

  /// \brief Get evidence associated with a specific sample group:
  const SVCandidateSetSequenceFragmentSampleGroup& getDataGroup(const unsigned bamIndex) const
  {
    data_t::const_iterator diter(_data.find(bamIndex));
    assert(diter != _data.end());
    return diter->second;
  }

  void clear()
  {
    _data.clear();
    _searchIntervals.clear();
  }

  /// return true if this search interval overlaps with any previous
  /// intervals
  ///
  /// if true, we know to expect repeated qnames
  ///
  bool setNewSearchInterval(const GenomeInterval& newSearch);

private:
  typedef std::map<unsigned, SVCandidateSetSequenceFragmentSampleGroup> data_t;
  data_t                                                                _data;

  std::vector<GenomeInterval> _searchIntervals;
};
