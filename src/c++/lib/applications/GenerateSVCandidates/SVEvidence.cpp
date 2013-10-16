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
/// \author Chris Saunders and Xiaoyu Chen
///

#include "SVEvidence.hh"

#include <iostream>



std::ostream&
operator<<(
    std::ostream& os,
    const SVFragmentEvidenceAlleleBreakendPerRead& svbpr)
{
    os << "isEval: " << svbpr.isSplitEvaluated
       << " isSplitSupport: " << svbpr.isSplitSupport
       << " splitEvidence: " << svbpr.splitEvidence
       << " splitLnLhood: " << svbpr.splitLnLhood;

    return os;
}




std::ostream&
operator<<(
    std::ostream& os,
    const SVFragmentEvidenceAlleleBreakend& svbp)
{
    os << "isFrag: " << svbp.isFragmentSupport << " fragProb: " << svbp.fragLengthProb << "\n";
    os << "read1ev: " << svbp.read1 << "\n";
    os << "read2ev: " << svbp.read2 << "\n";
    return os;
}



std::ostream&
operator<<(
    std::ostream& os,
    const SVFragmentEvidenceAllele& sval)
{
    os << "----BP1: " << sval.bp1;
    os << "----BP2: " << sval.bp2;
    return os;
}



std::ostream&
operator<<(
    std::ostream& os,
    const SVFragmentEvidenceRead& svr)
{
    os << "readinfo isScanned: " << svr.isScanned << " isAnchored: " << svr.isAnchored << " mapq: " << svr.mapq;
    return os;
}



std::ostream&
operator<<(
    std::ostream& os,
    const SVFragmentEvidence& sve)
{
    os << "FRAGMENT_START\n"
       << "read1: " << sve.read1 << "\n"
       << "read2: " << sve.read2 << "\n"
       << "+++++++++++ALT\n" << sve.alt
       << "+++++++++++REF\n" << sve.ref
       << "FRAGMENT_END\n";

    return os;
}
