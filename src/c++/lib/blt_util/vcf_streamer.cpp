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

#include "blt_util/blt_exception.hh"
#include "blt_util/log.hh"
#include "blt_util/seq_util.hh"
#include "blt_util/vcf_streamer.hh"

#include <cassert>
#include <cstdlib>
#include <sys/stat.h>

#include <iostream>
#include <set>
#include <string>



// return true only if all chromosomes in the vcf exist in the
// bam header
static
void
check_vcf_header_compatability(const char* vcf_filename,
                               const ti_index_t* vh,
                               const bam_header_t* bh)
{

    assert(NULL != bh);
    assert(NULL != vh);

    // build set of chrom labels from BAM:
    std::set<std::string> bamlabels;
    for (int32_t i(0); i<bh->n_targets; ++i)
    {
        bamlabels.insert(std::string(bh->target_name[i]));
    }
    int n_vcf_labels(0);
    const char** vcf_labels = ti_seqname(vh, &n_vcf_labels);

    for (int i(0); i<n_vcf_labels; ++i)
    {
        if (bamlabels.find(std::string(vcf_labels[i])) == bamlabels.end())
        {
            log_os << "ERROR: Chromosome label '" << vcf_labels[i] << "' in VCF file '" << vcf_filename << "' does not exist in the BAM header\n";
            exit(EXIT_FAILURE);
        }
    }

    free(vcf_labels);
}



vcf_streamer::
vcf_streamer(const char* filename,
             const char* region,
             const bam_header_t* bh)
    : _is_record_set(false), _is_stream_end(false), _record_no(0), _stream_name(filename),
      _tfp(NULL), _titer(NULL)
{

    if (NULL == filename)
    {
        throw blt_exception("vcf filename is null ptr");
    }

    if ('\0' == *filename)
    {
        throw blt_exception("vcf filename is empty string");
    }

    _tfp = ti_open(filename, 0);

    if (NULL == _tfp)
    {
        log_os << "ERROR: Failed to open VCF file: '" << filename << "'\n";
        exit(EXIT_FAILURE);
    }

    // read from a specific region:
    if (ti_lazy_index_load(_tfp) < 0)
    {
        log_os << "ERROR: Failed to load index for VCF file: '" << filename << "'\n";
        exit(EXIT_FAILURE);
    }

    if (NULL != bh)
    {
        check_vcf_header_compatability(filename,_tfp->idx,bh);
    }

    if (NULL == region)
    {
        // read the whole VCF file:
        _titer = ti_query(_tfp, 0, 0, 0);
        return;
    }

    int tid,beg,end;
    if (ti_parse_region(_tfp->idx, region, &tid, &beg, &end) == 0)
    {
        _titer = ti_queryi(_tfp, tid, beg, end);
    }
    else
    {
        _is_stream_end=true;
    }
}



vcf_streamer::
~vcf_streamer()
{
    if (NULL != _titer) ti_iter_destroy(_titer);
    if (NULL != _tfp) ti_close(_tfp);
}



bool
vcf_streamer::
next(const bool is_indel_only)
{
    if (_is_stream_end || (NULL==_tfp) || (NULL==_titer)) return false;

    while (true)
    {
        int len;
        const char* vcf_record_string(ti_read(_tfp, _titer, &len));

        _is_stream_end=(NULL == vcf_record_string);
        _is_record_set=(! _is_stream_end);
        if (! _is_record_set) break;
        _record_no++;

        if (! _vcfrec.set(vcf_record_string))
        {
            log_os << "ERROR: Can't parse vcf record: '" << vcf_record_string << "'\n";
            exit(EXIT_FAILURE);
        }
        if (! _vcfrec.is_valid()) continue;
        if (is_indel_only && (! _vcfrec.is_indel())) continue;

        break; // found expected vcf record type
    }

    return _is_record_set;
}



void
vcf_streamer::
report_state(std::ostream& os) const
{

    const vcf_record* vcfp(get_record_ptr());

    os << "\tvcf_stream_label: " << name() << "\n";
    if (NULL != vcfp)
    {
        os << "\tvcf_stream_record_no: " << record_no() << "\n"
           << "\tvcf_record: " << *(vcfp) << "\n";
    }
    else
    {
        os << "\tno vcf record currently set\n";
    }
}
