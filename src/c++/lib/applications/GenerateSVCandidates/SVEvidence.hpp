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
/// \author Chris Saunders and Xiaoyu Chen
///
///
/// \brief classes required to build the SVEvidence object
///
/// Note that the SVEvidence object is a kind of mini-database, accumulating all information about the
/// relationship of each fragment with a specific SV candidate. This allows us to tie together:
///
/// 1. spanning read information from the full fragment
/// 2. split read information from read1
/// 3. split read information from read2
/// 4. fragment mapping properties from read1 & read2
///
/// ...for each fragment. The object stores this per-fragment information for all fragments impacting
/// all alleles at a given locus (But note at present we limit the alternate alleles to one).
///

#pragma once

#include <cassert>
#include <iosfwd>
#include <map>
#include <string>
#include <vector>

/// For a single read from a read pair, track all support data specific to an individual breakend of a single
/// allele
///
struct SVFragmentEvidenceAlleleBreakendPerRead {
  /// True if the read has been evaluated for support across the given allele breakend
  bool isSplitEvaluated = false;

  /// If isSplitEvaluated is true, this indicates if the read supports the given allele breakend
  bool isSplitSupport = false;

  /// If isSplitEvaluated is true, this indicates if the read supports the given allele breakend by a more
  /// permissive criteria
  bool isTier2SplitSupport = false;

  /// If isSplitEvaluated is true, this provides the read's evidence score
  float splitEvidence = 0;

  /// If isSplitEvaluated is true, this is the log likelihood the best alignment across the given breakend
  float splitLnLhood = 0;
};

std::ostream& operator<<(std::ostream& os, const SVFragmentEvidenceAlleleBreakendPerRead& svbpr);

/// Track all support data from an individual fragment specific to an individual breakend of a single allele
///
struct SVFragmentEvidenceAlleleBreakend {
  SVFragmentEvidenceAlleleBreakendPerRead& getRead(const bool isRead1) { return (isRead1 ? read1 : read2); }

  const SVFragmentEvidenceAlleleBreakendPerRead& getRead(const bool isRead1) const
  {
    return (isRead1 ? read1 : read2);
  }

  void clearPairSupport()
  {
    isFragmentSupport = false;
    fragLengthProb    = 0;
  }

  /// If true, paired-read analysis shows that this read pair fragment supports this allele on this breakend
  bool isFragmentSupport = false;

  /// If isFragmentSupport is true, this is the prob of the fragment size given the specified allele breakend
  float fragLengthProb = 0;

  SVFragmentEvidenceAlleleBreakendPerRead read1;  ///< read1 specific evidence
  SVFragmentEvidenceAlleleBreakendPerRead read2;  ///< read2 specific evidence
};

std::ostream& operator<<(std::ostream& os, const SVFragmentEvidenceAlleleBreakend& svbp);

/// Track all support data from an individual fragment specific to a single allele of an SV candidate
///
struct SVFragmentEvidenceAllele {
  SVFragmentEvidenceAlleleBreakend& getBp(const bool isBp1) { return (isBp1 ? bp1 : bp2); }

  /// \return A tuple of boolean values indicating if the read supports breakend1 or breakend2, respectively,
  /// of the given allele.
  std::pair<bool, bool> isAnySplitReadSupport(const bool isRead1) const
  {
    return std::make_pair(bp1.getRead(isRead1).isSplitSupport, bp2.getRead(isRead1).isSplitSupport);
  }

  std::pair<bool, bool> isAnyTier2SplitReadSupport(const bool isRead1) const
  {
    return std::make_pair(bp1.getRead(isRead1).isTier2SplitSupport, bp2.getRead(isRead1).isTier2SplitSupport);
  }

  void clearPairSupport()
  {
    bp1.clearPairSupport();
    bp2.clearPairSupport();
  }

  SVFragmentEvidenceAlleleBreakend bp1;
  SVFragmentEvidenceAlleleBreakend bp2;
};

std::ostream& operator<<(std::ostream& os, const SVFragmentEvidenceAllele& sval);

/// Store properties of the reads in a fragment which are not tightly coupled to any one allele/bp, etc....
///
struct SVFragmentEvidenceRead {
  /// TODO set anchor policy wrt shadow state!!!
  bool isAnchored(const bool isTier2) const { return (isTier2 ? _isTier2Anchored : _isAnchored); }

  bool isObservedAnchor(const bool isTier2) const { return (isScanned && isAnchored(isTier2)); }

  void setAnchored(const bool val) { _isAnchored = val; }

  void setTier2Anchored(const bool val) { _isTier2Anchored = val; }

  /// If true, this read's bam record has been scanned to fill in the remaining values in this object
  bool isScanned = false;

  /// Read was originally unmapped but had a mapped mate read, mapq is MAPQ of the mate in this case
  bool isShadow = false;

  unsigned mapq = 0;
  unsigned size = 0;

private:
  /// If true, the read is found and known to have a confident mapping wrt fragment support
  bool _isAnchored = false;

  /// If true, the read is found and known to have a confident mapping wrt fragment support at tier2
  bool _isTier2Anchored = false;
};

std::ostream& operator<<(std::ostream& os, const SVFragmentEvidenceRead& svr);

/// Track all support data from an individual fragment specific to an SV hypothesis
///
/// This is both to prevent double-counting of evidence and to consolidate different
/// sources of information (paired-split, etc).
///
struct SVFragmentEvidence {
  SVFragmentEvidenceRead& getRead(const bool isRead1) { return (isRead1 ? read1 : read2); }

  const SVFragmentEvidenceRead& getRead(const bool isRead1) const { return (isRead1 ? read1 : read2); }

  /// \return True if this fragment provides paired read support for either breakend of either allele
  bool isAnySpanningPairSupport() const
  {
    const bool isRefSupport(ref.bp1.isFragmentSupport || ref.bp2.isFragmentSupport);
    const bool isAltSupport(alt.bp1.isFragmentSupport || alt.bp2.isFragmentSupport);

    return (isRefSupport || isAltSupport);
  }

  /// does this fragment read provide any pair evidence for any bp of the ALT allele?
  bool isAltSpanningPairSupport() const
  {
    const bool isAltSupport(alt.bp1.isFragmentSupport || alt.bp2.isFragmentSupport);

    return (isAltSupport);
  }

  /// Indicates if this fragment read provides split evidence for any allele/bp combination
  ///
  /// \return A tuple of boolean values indicating if the read supports breakend1 or breakend2, respectively,
  /// of either the ref or alt alleles.
  std::pair<bool, bool> isAnySplitReadSupport(const bool isRead1) const
  {
    const std::pair<bool, bool> isAlt(alt.isAnySplitReadSupport(isRead1));
    const std::pair<bool, bool> isRef(ref.isAnySplitReadSupport(isRead1));

    return std::make_pair((isAlt.first || isRef.first), (isAlt.second || isRef.second));
  }

  /// does this fragment read provide any split evidence for any bp of the ALT allele?
  bool isAltSplitReadSupport(const bool isRead1) const
  {
    const std::pair<bool, bool> isAlt(alt.isAnySplitReadSupport(isRead1));

    return (isAlt.first || isAlt.second);
  }

  /// does this fragment read provide any split evidence for any allele/bp combination?
  std::pair<bool, bool> isAnyTier2SplitReadSupport(const bool isRead1) const
  {
    const std::pair<bool, bool> isAlt(alt.isAnyTier2SplitReadSupport(isRead1));
    const std::pair<bool, bool> isRef(ref.isAnyTier2SplitReadSupport(isRead1));

    return std::make_pair((isAlt.first || isRef.first), (isAlt.second || isRef.second));
  }

  /// does this fragment read provide any split evidence for any bp of the ALT allele?
  bool isAltTier2SplitReadSupport(const bool isRead1) const
  {
    const std::pair<bool, bool> isAlt(alt.isAnyTier2SplitReadSupport(isRead1));

    return (isAlt.first || isAlt.second);
  }

  void clearPairSupport()
  {
    ref.clearPairSupport();
    alt.clearPairSupport();
  }

  SVFragmentEvidenceRead read1;
  SVFragmentEvidenceRead read2;

  SVFragmentEvidenceAllele alt;
  SVFragmentEvidenceAllele ref;
};

std::ostream& operator<<(std::ostream& os, const SVFragmentEvidence& sve);

/// Track all support data for an SV hypothesis
///
/// Note how this object is different than SVScoreInfoSomatic -- it is highly detailed and meant to be
/// processed to create summary statistics and scores later. Those scores and summary statistics should go
/// into objects like SomaticSVSCoreInfo to be written out in whichever output format is selected.
///
struct SVEvidence {
  typedef std::map<std::string, SVFragmentEvidence> evidenceTrack_t;

  unsigned size() const { return samples.size(); }

  evidenceTrack_t& getSampleEvidence(const unsigned index)
  {
    assert(index < size());
    return samples[index];
  }

  const evidenceTrack_t& getSampleEvidence(const unsigned index) const
  {
    assert(index < size());
    return samples[index];
  }

  std::vector<evidenceTrack_t> samples;
};
