// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
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


#include "bam_seq_read_util.hh"



static
void
get_read_align_strand_end_skip(const bam_seq& bseq,
                               unsigned& end_skip)
{
    unsigned read_end(bseq.size());

    while (read_end>0)
    {
        if (bseq.get_char(read_end-1)=='N') read_end--;
        else break;
    }

    end_skip=bseq.size()-read_end;
}



void
get_read_fwd_strand_skip(const bam_seq& bseq,
                         const bool is_fwd_strand,
                         unsigned& begin_skip,
                         unsigned& end_skip)
{
    begin_skip=0;
    if (is_fwd_strand)
    {
        get_read_align_strand_end_skip(bseq,end_skip);
    }
    else
    {
        end_skip=0;
        const unsigned bsize(bseq.size());
        while (begin_skip<bsize)
        {
            if (bseq.get_char(begin_skip)=='N') begin_skip++;
            else break;
        }
    }
}
