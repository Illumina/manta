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

/// \brief bam record manipulation functions
///

#pragma once


extern "C" {
#define __STDC_LIMIT_MACROS
#include "bam.h"
}


namespace BAM_FLAG
{
enum index_t {
    PAIRED        = 0x001,
    PROPER_PAIR   = 0x002,
    UNMAPPED      = 0x004,
    MATE_UNMAPPED = 0x008,
    STRAND        = 0x010,
    MATE_STRAND   = 0x020,
    FIRST_READ    = 0x040,
    SECOND_READ   = 0x080,
    SECONDARY     = 0x100,
    FILTER        = 0x200,
    DUPLICATE     = 0x400,
    SUPPLEMENT    = 0x800
};
}


// insert new qname
//
void
edit_bam_qname(const char* name,
               bam1_t& br);

// Read length is taken from the read string. Assumes offset has
// already been removed from qual
void
edit_bam_read_and_quality(const char* read,
                          const uint8_t* qual,
                          bam1_t& br);

// remove all copies of optional field "tag"
//
void
nuke_bam_aux_field(bam1_t& br,
                   const char* tag);


// store an unsigned int to optional field "tag", optimize storage
// based on size of x
void
bam_aux_append_unsigned(bam1_t& br,
                        const char* tag,
                        uint32_t x);

// change the size of a subsegment of the bam data, 'end' identifies
// the byte offset of the end of the segment and 'delta' is the change
// to the segment size
//
void
change_bam_data_segment_len(const int end,
                            const int delta,
                            bam1_t& br);

// Update bam record bin value -- call after updating pos and/or
// cigar fields.
//
inline
void
bam_update_bin(bam1_t& br)
{
    // set bin value:
    //
    // Test for position rather than looking at the unmapped flag
    // because we want to index shadow reads.
    //
    bam1_core_t& brc(br.core);
    if (brc.pos>=0)
    {
        if (brc.n_cigar!=0)
        {
            // normal case:
            brc.bin = bam_reg2bin(brc.pos, bam_calend(&brc, bam1_cigar(&br)));
        }
        else
        {
            // shadow case: (match logic from samtools)
            brc.bin = bam_reg2bin(brc.pos, brc.pos+1);
        }
    }
    else
    {
        // unmapped, non-shadow reads:
        brc.bin = 0;
    }
}
