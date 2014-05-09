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
/// \author Chris Saunders and Xiaoyu Chen
///

#include "SVScorerShared.hh"



void
setReadEvidence(
    const unsigned minMapQ,
    const unsigned minTier2MapQ,
    const unsigned mapq,
    const unsigned readSize,
    const bool isShadow,
    SVFragmentEvidenceRead& read)
{
    if (read.isScanned) return;

    read.isScanned = true;
    read.mapq = mapq;
    read.isShadow = isShadow;
    read.setAnchored(read.mapq >= minMapQ);
    read.setTier2Anchored(read.mapq >= minTier2MapQ);
    read.size = readSize;
}
