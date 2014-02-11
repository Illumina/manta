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

#include "blt_util/SimpleAlignment.hh"
#include "blt_util/align_path_bam_util.hh"



SimpleAlignment::
SimpleAlignment(const bam_record& bamRead) :
    is_fwd_strand(bamRead.is_fwd_strand()),
    tid(bamRead.target_id()),
    pos(bamRead.pos()-1)
{
    bam_cigar_to_apath(bamRead.raw_cigar(),bamRead.n_cigar(),path);
}



/// convert segment_type to match if the segment exists before or after all match segments currently in the alignment
///
SimpleAlignment
matchifyEdgeSegmentType(
    const SimpleAlignment& al,
    const ALIGNPATH::align_t segmentType,
    const bool isMatchLeadingEdge,
    const bool isMatchTrailingEdge)
{
    using namespace ALIGNPATH;

    assert(is_segment_type_read_length(segmentType));

    SimpleAlignment al2;
    al2.is_fwd_strand=al.is_fwd_strand;
    al2.tid=al.tid;
    al2.pos=al.pos;

    const std::pair<unsigned,unsigned> ends(get_match_edge_segments(al.path));
    const unsigned as(al.path.size());
    for (unsigned i(0); i<as; ++i)
    {
        const path_segment& ps(al.path[i]);
        const bool isLeadingEdgeSegment(i<ends.first);
        const bool isTrailingEdgeSegment(i>ends.second);
        const bool isTargetType(ps.type==segmentType);
        const bool isCandidateEdge((isMatchLeadingEdge && isLeadingEdgeSegment) ||
                                   (isMatchTrailingEdge && isTrailingEdgeSegment));
        const bool isEdgeTarget(isCandidateEdge && isTargetType);
        if (isEdgeTarget && isLeadingEdgeSegment) al2.pos-=ps.length;
        if (isEdgeTarget || (ps.type==MATCH))
        {
            if ((! al2.path.empty()) && (al2.path.back().type == MATCH))
            {
                al2.path.back().length += ps.length;
            }
            else
            {
                al2.path.push_back(ps);
                al2.path.back().type = MATCH;
            }
        }
        else
        {
            al2.path.push_back(ps);
        }
    }

    return al2;
}
