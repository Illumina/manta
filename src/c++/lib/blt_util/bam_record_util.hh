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
///
///

#pragma once

#include "blt_util/align_path.hh"
#include "blt_util/bam_record.hh"


/// is this read part of mapped pair with 'Innie' orientation?
///
/// Note this does not test MAPQ or fragment size, but could
/// be used as the core of a 'proper-pair' predicate
bool
is_innie_pair(
    const bam_record& bam_read);

/// return average basecall qscore for this read
unsigned
get_avg_quality(
    const bam_record& bam_read);


/// A struct which holds minimal alignment information processed
// from a BAM record or other source
struct SimpleAlignment
{
    SimpleAlignment() :
        is_fwd_strand(true),
        tid(0),
        pos(0)
    {}

    SimpleAlignment(const bam_record& bamRead);

    bool is_fwd_strand;
    int32_t tid;
    pos_t pos;
    ALIGNPATH::path_t path;
};

