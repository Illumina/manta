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
/// \author Chris Saunders
///

#pragma once

#include "hts_streamer.hh"
#include "vcf_record.hh"


struct vcf_streamer : private hts_streamer
{
    // optionally provide a BAM header to validate vcf chromosome names against
    //
    vcf_streamer(
        const char* filename,
        const char* region,
        const bam_header_t* bh = nullptr);

    ~vcf_streamer();

    // advance to next vcf record
    //
    // is_indel_only - if set, skip all records except indels
    //
    bool next(const bool is_indel_only=false);

    const vcf_record*
    get_record_ptr() const
    {
        if (_is_record_set) return &_vcfrec;
        else                return nullptr;
    }

    void report_state(std::ostream& os) const;

private:
    bcf_hdr_t* _hdr;
    vcf_record _vcfrec;
};
