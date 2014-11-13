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

#include "blt_util/blt_exception.hh"
#include "blt_util/log.hh"
#include "htsapi/bam_streamer.hh"

#include <cassert>
#include <cstdlib>

#include <fstream>
#include <iostream>



bam_streamer::
bam_streamer(
    const char* filename,
    const char* region)
    : _is_record_set(false),
      _bfp(nullptr),
      _hidx(nullptr),
      _hitr(nullptr),
      _record_no(0),
      _stream_name(filename),
      _is_region(false)
{
    assert(nullptr != filename);
    if ('\0' == *filename)
    {
        throw blt_exception("Can't initialize bam_streamer with empty filename\n");
    }

    _bfp = samopen(filename, "rb", 0);

    if (nullptr == _bfp)
    {
        log_os << "ERROR: Failed to open SAM/BAM/CRAM file: " << filename << "\n";
        exit(EXIT_FAILURE);
    }

    if (nullptr == region)
    {
        // read the whole BAM file:

        if (_bfp->header->n_targets)
        {
            // parse a fake region so that header->hash is created
            std::string fake_region(target_id_to_name(0));
            fake_region += ":1-1";
            int ref,beg,end;
            bam_parse_region(_bfp->header, fake_region.c_str(), &ref, &beg, &end);
        }
        return;
    }

    // read a specific region of the bam file:
    set_new_region(region);
}



bam_streamer::
~bam_streamer()
{
    if (nullptr != _hitr) hts_itr_destroy(_hitr);
    if (nullptr != _hidx) hts_idx_destroy(_hidx);
    if (nullptr != _bfp) samclose(_bfp);
}



static
bool
fexists(const char* filename)
{
    std::ifstream ifile(filename);
    return ifile;
}



static
bool
hasEnding(
    const std::string& fullString,
    const std::string& ending)
{
    if (fullString.length() < ending.length()) return false;
    return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
}



// load index if it hasn't been set already:
void
bam_streamer::
_load_index()
{
    /// TODO: Find out whether _hidx can be destroyed after the HTS
    /// iterator is created, in which case this could be a local
    /// variable. Until we know, _hidx should persist for the lifetime
    /// of _hiter
    if (nullptr != _hidx) return;

    std::string index_base(name());

    // hack to allow GATK/Picard bai name convention:
    if ((! fexists((index_base+".bai").c_str())) &&
        (! fexists((index_base+".csa").c_str())) &&
        (! fexists((index_base+".crai").c_str())))
    {
        static const std::string bamext(".bam");
        if (hasEnding(index_base,bamext))
        {
            index_base=index_base.substr(0,index_base.length()-bamext.length());
        }
    }

    _hidx = sam_index_load(_bfp->file, index_base.c_str());
    if (nullptr == _hidx)
    {
        log_os << "ERROR: BAM/CRAM index is not available for file: " << name() << "\n";
        exit(EXIT_FAILURE);
    }
}



void
bam_streamer::
set_new_region(const char* region)
{
    int ref,beg,end;
    bam_parse_region(_bfp->header, region, &ref, &beg, &end); // parse the region

    set_new_region(ref,beg,end);
    _region=region;
}



void
bam_streamer::
set_new_region(const int ref, const int beg, const int end)
{
    if (nullptr != _hitr) hts_itr_destroy(_hitr);

    _load_index();

    if (ref < 0)
    {
        log_os << "ERROR: Invalid region specified for BAM file: " << name() << "\n";
        exit(EXIT_FAILURE);
    }

    _hitr = sam_itr_queryi(_hidx,ref,beg,end);
    _is_region = true;
    _region.clear();

    _is_record_set = false;
    _record_no = 0;
}


bool
bam_streamer::
next()
{
    if (nullptr == _bfp) return false;

    int ret;
    if (nullptr == _hitr)
    {
        ret = samread(_bfp, _brec._bp);
    }
    else
    {
        ret = sam_itr_next(_bfp->file, _hitr, _brec._bp);
    }

    _is_record_set=(ret >= 0);
    if (_is_record_set) _record_no++;

    return _is_record_set;
}



const char*
bam_streamer::
target_id_to_name(const int32_t tid) const
{
    // assert(tid < _bfp->header->n_targets);
    if (tid<0)
    {
        static const char unmapped[] = "*";
        return unmapped;
    }
    return _bfp->header->target_name[tid];
}



int32_t
bam_streamer::
target_name_to_id(const char* seq_name) const
{
    return bam_get_tid(_bfp->header,seq_name);
}



void
bam_streamer::
report_state(std::ostream& os) const
{
    const bam_record* bamp(get_record_ptr());

    os << "\tbam_stream_label: " << name() << "\n";
    if (_is_region && (! _region.empty()))
    {
        os << "\tbam_stream_selected_region: " << _region << "\n";
    }
    if (nullptr != bamp)
    {
        os << "\tbam_stream_record_no: " << record_no() << "\n";
        os << "\tbam_record QNAME/read_number: " << bamp->qname() << "/" << bamp->read_no() << "\n";
        const char* chrom_name(target_id_to_name(bamp->target_id()));
        os << "\tbam record RNAME: " << chrom_name << "\n";
        os << "\tbam record POS: " << bamp->pos() << "\n";

    }
    else
    {
        os << "\tno bam record currently set\n";
    }
}
