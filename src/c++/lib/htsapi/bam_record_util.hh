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
///
///

#pragma once

#include "blt_util/align_path.hh"
#include "blt_util/SimpleAlignment.hh"
#include "htsapi/bam_record.hh"


/// is this read part of a pair where both members are mapped?
bool
is_mapped_pair(
    const bam_record& bam_read);


/// is this read part of a pair where both members are mapped to the same chrom?
///
bool
is_mapped_chrom_pair(
    const bam_record& bam_read);

/// is this read part of mapped pair with 'Innie' orientation?
///
/// Note this does not test MAPQ or fragment size, but could
/// be used as the core of a 'proper-pair' predicate
bool
is_innie_pair(
    const bam_record& bam_read);

/// detect cases where the dna fragment is likely to be shorter than read length
///
/// note this is an approximation because it's based on a single bam record, an
/// exact answer would require both records in the pair. In practice, this should
/// be good enough.
bool
is_possible_adapter_pair(
    const bam_record& bam_read,
    const SimpleAlignment& matchedAlignment);

/// return average basecall qscore for this read
unsigned
get_avg_quality(
    const bam_record& bam_read);

/// select 'first' read in pair such that you
/// consistently get only one read per-pair
/// (assuming the bam file is properly formated)
inline
bool
isFirstRead(
    const bam_record& bamRead)
{
    if (bamRead.pos() < bamRead.mate_pos()) return true;
    if ((bamRead.pos() == bamRead.mate_pos()) && bamRead.is_first()) return true;
    return false;
}

/// get BAM RG tag, return an empty string "" if no RG tag exists:
inline
const char*
getReadGroup(
    const bam_record& bamRead)
{
    static const char defaultRG[] = "";
    static const char rgTag[] = {'R','G'};

    const char* rgStr(bamRead.get_string_tag(rgTag));

    return ((NULL == rgStr) ? defaultRG : rgStr);
}
