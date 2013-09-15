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
/// \author Chris Saunders
///


#include "depth_buffer_util.hh"

#include "boost/foreach.hpp"



void
add_alignment_to_depth_buffer(
    const pos_t& pos,
    const ALIGNPATH::path_t& apath,
    depth_buffer& db)
{
    using namespace ALIGNPATH;

    pos_t ref_head_pos(pos);

    BOOST_FOREACH(const path_segment& ps, apath)
    {
        if ( is_segment_align_match(ps.type) )
        {
            for (unsigned j(0); j<ps.length; ++j) db.inc(ref_head_pos+static_cast<pos_t>(j));
        }

        if ( is_segment_type_ref_length(ps.type) ) ref_head_pos += ps.length;
    }
}

