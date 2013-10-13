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

#pragma once

#include "blt_util/bam_record.hh"
#include "SVEvidence.hh"


void
setReadEvidence(
    const unsigned minMapQ,
    const bam_record& bamRead,
    SVFragmentEvidenceRead& read);
