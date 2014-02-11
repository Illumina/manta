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


#include "blt_util/bam_util.hh"
#include "blt_util/tabix_util.hh"
#include "blt_util/vcf_record.hh"


struct vcf_streamer
{

    // optionally provide a BAM header to validate vcf chromosome names against
    //
    explicit
    vcf_streamer(const char* filename,
                 const char* region = NULL,
                 const bam_header_t* bh = NULL);

    ~vcf_streamer();

    // advance to next vcf record
    //
    // is_indel_only - if set, skip all records except indels
    //
    bool next(const bool is_indel_only=false);

    const vcf_record* get_record_ptr() const
    {
        if (_is_record_set) return &_vcfrec;
        else               return NULL;
    }

    const char* name() const
    {
        return _stream_name.c_str();
    }

    unsigned record_no() const
    {
        return _record_no;
    }

    void report_state(std::ostream& os) const;

    //const bam_header_t*
    //get_header() const { return _bfp->header; }

private:
    bool _is_record_set;
    bool _is_stream_end;
    unsigned _record_no;
    std::string _stream_name;

    tabix_t* _tfp;
    ti_iter_t _titer;

    vcf_record _vcfrec;
};

