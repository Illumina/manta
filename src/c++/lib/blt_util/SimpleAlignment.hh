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

#pragma once

#include "blt_util/align_path.hh"
#include "blt_util/known_pos_range2.hh"


/// holds minimal alignment information processed
/// from a BAM record or other source
struct SimpleAlignment
{
    bool is_fwd_strand = true;
    int32_t tid = 0;
    pos_t pos = 0;
    ALIGNPATH::path_t path;
};


/// convert segment_type to match if the segment exists before or after all match segments currently in the alignment
///
SimpleAlignment
matchifyEdgeSegmentType(
    const SimpleAlignment& al,
    const ALIGNPATH::align_t segmentType,
    const bool isMatchLeadingEdge = true,
    const bool isMatchTrailingEdge = true);


/// transform an alignment such that any soft-clipped edge segments
/// become match.
///
/// segments are joined and start pos is adjusted appropriately
///
inline
SimpleAlignment
matchifyEdgeSoftClip(const SimpleAlignment& al)
{
    return matchifyEdgeSegmentType(al, ALIGNPATH::SOFT_CLIP);
}


/// get the range in reference coordinates if you did run matchifyEdgeSoftClip on an alignment:
known_pos_range2
matchifyEdgeSoftClipRefRange(const SimpleAlignment& al);
