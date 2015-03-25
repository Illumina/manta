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
/// \author Chris Saunders and Xiaoyu Chen
///

#pragma once

#include "htsapi/bam_record.hh"
#include "SVEvidence.hh"


void
setReadEvidence(
    const unsigned minMapQ,
    const unsigned minTier2MapQ,
    const unsigned mapq,
    const unsigned readSize,
    const bool isShadow,
    SVFragmentEvidenceRead& read);


inline
void
setReadEvidence(
    const unsigned minMapQ,
    const unsigned minTier2MapQ,
    const bam_record& bamRead,
    const bool isShadow,
    SVFragmentEvidenceRead& read)
{
    setReadEvidence(minMapQ, minTier2MapQ,
                    bamRead.map_qual(), bamRead.read_size(),
                    isShadow, read);
}
