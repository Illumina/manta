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

#include "SimpleAlignment_bam_util.hh"
#include "htsapi/align_path_bam_util.hh"



SimpleAlignment
getAlignment(
    const bam_record& bamRead)
{
    SimpleAlignment al;
    al.is_fwd_strand=bamRead.is_fwd_strand();
    al.tid=bamRead.target_id();
    al.pos=(bamRead.pos()-1);

    bam_cigar_to_apath(bamRead.raw_cigar(),bamRead.n_cigar(),al.path);
    return al;
}


SimpleAlignment
getFakeMateAlignment(
    const bam_record& bamRead)
{
    SimpleAlignment al;
    assert(! bamRead.is_mate_unmapped());
    al.is_fwd_strand=bamRead.is_mate_fwd_strand();
    al.tid=bamRead.mate_target_id();
    al.pos=(bamRead.mate_pos()-1);
    al.path.emplace_back(ALIGNPATH::MATCH, bamRead.read_size());
    return al;
}
