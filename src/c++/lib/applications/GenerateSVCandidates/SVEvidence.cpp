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
/// \author Chris Saunders and Xiaoyu Chen
///

#include "SVEvidence.hpp"

#include <iostream>

std::ostream& operator<<(std::ostream& os, const SVFragmentEvidenceAlleleBreakendPerRead& svbpr)
{
  os << "isEval: " << svbpr.isSplitEvaluated << " isSplitSupport: " << svbpr.isSplitSupport
     << " isTier2SplitSupport: " << svbpr.isTier2SplitSupport << " splitEvidence: " << svbpr.splitEvidence
     << " splitLnLhood: " << svbpr.splitLnLhood;

  return os;
}

std::ostream& operator<<(std::ostream& os, const SVFragmentEvidenceAlleleBreakend& svbp)
{
  os << "isFrag: " << svbp.isFragmentSupport << " fragProb: " << svbp.fragLengthProb << "\n";
  os << "read1ev: " << svbp.read1 << "\n";
  os << "read2ev: " << svbp.read2 << "\n";
  return os;
}

std::ostream& operator<<(std::ostream& os, const SVFragmentEvidenceAllele& sval)
{
  os << "----BP1: " << sval.bp1;
  os << "----BP2: " << sval.bp2;
  return os;
}

std::ostream& operator<<(std::ostream& os, const SVFragmentEvidenceRead& svr)
{
  os << "readinfo isScanned: " << svr.isScanned << " isAnchored: " << svr.isAnchored(false)
     << " isTier2Anchored: " << svr.isAnchored(true) << " isShadow: " << svr.isShadow << " mapq: " << svr.mapq
     << " size: " << svr.size;
  return os;
}

std::ostream& operator<<(std::ostream& os, const SVFragmentEvidence& sve)
{
  os << "FRAGMENT_START\n"
     << "read1: " << sve.read1 << "\n"
     << "read2: " << sve.read2 << "\n"
     << "+++++++++++ALT\n"
     << sve.alt << "+++++++++++REF\n"
     << sve.ref << "FRAGMENT_END\n";

  return os;
}
