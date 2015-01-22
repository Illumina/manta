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
///
///

#include "bam_record_util.hh"
#include "align_path_bam_util.hh"



bool
is_mapped_chrom_pair(
    const bam_record& bam_read)
{
    if (! bam_read.is_paired()) return false;
    if (bam_read.is_unmapped() || bam_read.is_mate_unmapped()) return false;
    if (bam_read.target_id() != bam_read.mate_target_id()) return false;
    return true;
}



// note this is designed to return true for the common case
// of pos == mate_pos occurring for short FFPE fragments
//
bool
is_innie_pair(
    const bam_record& bam_read)
{
    if (! is_mapped_chrom_pair(bam_read)) return false;

    if     (bam_read.pos() < bam_read.mate_pos())
    {
        if (! bam_read.is_fwd_strand()) return false;
        if (  bam_read.is_mate_fwd_strand()) return false;
    }
    else if (bam_read.pos() > bam_read.mate_pos())
    {
        if (  bam_read.is_fwd_strand()) return false;
        if (! bam_read.is_mate_fwd_strand()) return false;
    }
    else
    {
        if (bam_read.is_fwd_strand() == bam_read.is_mate_fwd_strand()) return false;
    }

    return true;
}



// detect cases where the dna fragment size is likely to be less
// than the read length (likely here because without access to local
// and remote records we can't say with complete certainty.
//
bool
is_possible_adapter_pair(
    const bam_record& bam_read)
{
    if (! is_mapped_chrom_pair(bam_read)) return false;
    if (bam_read.is_fwd_strand() == bam_read.is_mate_fwd_strand()) return false;

    // the threshold of 5 here allows for a few bases of spurious adapter-reference match:
    return (std::abs(bam_read.pos()-bam_read.mate_pos()) < 5);
}



unsigned
get_avg_quality(
    const bam_record& bam_read)
{
    const unsigned len(bam_read.read_size());
    if (0 == len) return 0;

    const uint8_t* qual(bam_read.qual());
    unsigned sum(0);
    for (unsigned i(0); i<len; ++i)
    {
        sum+=qual[i];
    }
    // this does not capture the decimal remainder but well...
    return (sum/len);
}
