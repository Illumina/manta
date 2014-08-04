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

/// \file

/// \author Chris Saunders
///

#include "blt_util/bam_streamer.hh"
#include "blt_util/blt_exception.hh"
#include "blt_util/log.hh"

#include <cstdlib>

#include <fstream>
#include <iostream>



bam_streamer::
bam_streamer(const char* filename,
             const char* region)
    : _is_record_set(false), _bfp(NULL), _bidx(NULL), _biter(NULL),
      _record_no(0), _stream_name(filename), _is_region(false)
{

    assert(NULL != filename);
    if ('\0' == *filename)
    {
        throw blt_exception("Can't initialize bam_streamer with empty filename\n");
    }


    _bfp = samopen(filename, "rb", 0);

    if (NULL == _bfp)
    {
        log_os << "ERROR: Failed to open SAM/BAM file: " << filename << "\n";
        exit(EXIT_FAILURE);
    }


    if (NULL == region)
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
    if (NULL != _biter) bam_iter_destroy(_biter);
    if (NULL != _bidx) bam_index_destroy(_bidx);
    if (NULL != _bfp) samclose(_bfp);
}



static
bool fexists(const char* filename)
{
    std::ifstream ifile(filename);
    return ifile;
}



static
bool hasEnding(
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

    if (NULL != _bidx) return;

    // use the BAM index to read a region of the BAM file
    if (! (_bfp->type&0x01))
    {
        log_os << "ERROR: file must be in BAM format for region lookup: " << name() << "\n";
        exit(EXIT_FAILURE);
    }

    /// TODO: Find out whether _bidx can be destroyed after the BAM
    /// iterator is created, in which case this could be a local
    /// variable. Until we know, _bidx should persist for the lifetime
    /// of _biter
    std::string index_base(name());
    if (! fexists((index_base+".bai").c_str()))
    {
        static const std::string bamext(".bam");
        if (hasEnding(index_base,bamext))
        {
            index_base=index_base.substr(0,index_base.length()-bamext.length());
        }
    }

    _bidx = bam_index_load(index_base.c_str()); // load BAM index
    if (NULL == _bidx)
    {
        log_os << "ERROR: BAM index is not available for file: " << name() << "\n";
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

    if (NULL != _biter) bam_iter_destroy(_biter);

    _load_index();

    if (ref < 0)
    {
        log_os << "ERROR: Invalid region specified for BAM file: " << name() << "\n";
        exit(EXIT_FAILURE);
    }

    _biter = bam_iter_query(_bidx,ref,beg,end);
    _is_region = true;
    _region.clear();

    _is_record_set = false;
    _record_no = 0;
}


bool
bam_streamer::
next()
{
    if (NULL==_bfp) return false;

    int ret;
    if (NULL == _biter)
    {
        ret = samread(_bfp, _brec._bp);
    }
    else
    {
        ret = bam_iter_read(_bfp->x.bam, _biter, _brec._bp);
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
    if (NULL != bamp)
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
