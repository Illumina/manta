//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2019 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

/// \file
/// \author Chris Saunders
///

#include "htsapi/bam_streamer.hpp"
#include "blt_util/blt_exception.hpp"
#include "blt_util/log.hpp"
#include "htsapi/bam_header_util.hpp"

#include <cassert>
#include <cstdlib>

#include <fstream>
#include <iostream>
#include <sstream>

stream_state_reporter::~stream_state_reporter() {}

bam_streamer::bam_streamer(const char* filename, const char* referenceFilename, const char* region)
  : _is_record_set(false),
    _hfp(nullptr),
    _hdr(nullptr),
    _hidx(nullptr),
    _hitr(nullptr),
    _record_no(0),
    _stream_name(filename),
    _is_region(false)
{
  assert(nullptr != filename);
  if ('\0' == *filename) {
    throw blt_exception("Can't initialize bam_streamer with empty filename\n");
  }

  _hfp = hts_open(filename, "rb");

  if (nullptr == _hfp) {
    std::ostringstream oss;
    oss << "Failed to open SAM/BAM/CRAM file for reading: '" << name() << "'";
    throw blt_exception(oss.str().c_str());
  }

  if (nullptr != referenceFilename) {
    const std::string referenceFilenameIndex(std::string(referenceFilename) + ".fai");
    int               ret = hts_set_fai_filename(_hfp, referenceFilenameIndex.c_str());
    if (ret != 0) {
      std::ostringstream oss;
      oss << "Failed to use reference: '" << referenceFilename << "' for BAM/CRAM file: '" << name() << "'";
      throw blt_exception(oss.str().c_str());
    }
  }

  _hdr = sam_hdr_read(_hfp);

  if (nullptr == _hdr) {
    std::ostringstream oss;
    oss << "Failed to parse header from SAM/BAM/CRAM file: " << name();
    throw blt_exception(oss.str().c_str());
  }

  if (nullptr == region) {
    // setup to read the whole BAM file by default if resetRegion() is not called:
    if (_hdr->n_targets) {
      // parse any contig name so that header->hash is created
      // ignore returned tid value, so doesn't matter if fake name
      // exists
      target_name_to_id("fake_name");
    }
  } else {
    // read a specific region of the bam file:
    resetRegion(region);
  }
}

bam_streamer::~bam_streamer()
{
  if (nullptr != _hitr) hts_itr_destroy(_hitr);
  if (nullptr != _hidx) hts_idx_destroy(_hidx);
  if (nullptr != _hdr) bam_hdr_destroy(_hdr);
  if (nullptr != _hfp) {
    const int retval = hts_close(_hfp);
    if (retval != 0) {
      log_os << "ERROR: Failed to close BAM/CRAM file: '" << name() << "'\n";
      std::exit(EXIT_FAILURE);
    }
  }
}

static bool fexists(const char* filename)
{
  std::ifstream ifile(filename);
  return (!ifile.fail());
}

static bool hasEnding(const std::string& fullString, const std::string& ending)
{
  if (fullString.length() < ending.length()) return false;
  return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
}

// load index if it hasn't been set already:
void bam_streamer::_load_index()
{
  /// TODO: Find out whether _hidx can be destroyed after the HTS
  /// iterator is created, in which case this could be a local
  /// variable. Until we know, _hidx should persist for the lifetime
  /// of _hiter
  if (nullptr != _hidx) return;

  std::string index_base(name());

  // hack to allow GATK/Picard bai name convention:
  if ((!fexists((index_base + ".bai").c_str())) && (!fexists((index_base + ".csi").c_str())) &&
      (!fexists((index_base + ".crai").c_str()))) {
    static const std::string bamext(".bam");
    if (hasEnding(index_base, bamext)) {
      index_base = index_base.substr(0, index_base.length() - bamext.length());
    }
  }

  _hidx = sam_index_load(_hfp, index_base.c_str());
  if (nullptr == _hidx) {
    std::ostringstream oss;
    oss << "BAM/CRAM index is not available for file: '" << name() << "'";
    throw blt_exception(oss.str().c_str());
  }
}

void bam_streamer::resetRegion(const char* region)
{
  int32_t referenceContigId, beginPos, endPos;
  parse_bam_region_from_hdr(_hdr, region, referenceContigId, beginPos, endPos);

  try {
    resetRegion(referenceContigId, beginPos, endPos);
    _region = region;
  } catch (const std::exception& /*e*/) {
    log_os << "ERROR: exception while fetching BAM/CRAM region: '" << region << "' from file '" << name()
           << "'\n";
    throw;
  }
}

void bam_streamer::resetRegion(int referenceContigId, int beginPos, int endPos)
{
  if (nullptr != _hitr) hts_itr_destroy(_hitr);

  _load_index();

  if (referenceContigId < 0) {
    std::ostringstream oss;
    oss << "Invalid region (contig index: " << referenceContigId << ") specified for BAM/CRAM file: '"
        << name() << "'";
    throw blt_exception(oss.str().c_str());
  }

  _hitr = sam_itr_queryi(_hidx, referenceContigId, beginPos, endPos);
  if (_hitr == nullptr) {
    std::ostringstream oss;
    oss << "Failed to fetch region: #" << referenceContigId << ":" << beginPos << "-" << endPos
        << " specified for BAM/CRAM file: '" << name() << "'";
    throw blt_exception(oss.str().c_str());
  }
  _is_region = true;
  _region.clear();

  _is_record_set = false;
  _record_no     = 0;
}

bool bam_streamer::next()
{
  if (nullptr == _hfp) return false;

  int ret;
  if (nullptr == _hitr) {
    ret = sam_read1(_hfp, _hdr, _brec._bp);

    // Semi-documented sam_read1 API: -1 is expected read failure at end of stream, any other negative value
    // is an error
    if (ret < -1) {
      std::ostringstream oss;
      oss << "Unexpected return value from htslib sam_read1 function '" << ret
          << "' while attempting to read BAM/CRAM file:\n";
      report_state(oss);
      throw blt_exception(oss.str().c_str());
    }
  } else {
    ret = sam_itr_next(_hfp, _hitr, _brec._bp);

    // Re sam_itr_next API: -1 is expected read failure at end of stream. As of htslib v1.5 errors also give a
    // return value of -1. If PR #575 is accepted then errors should return a value less than -1.
    if (ret < -1) {
      std::ostringstream oss;
      oss << "Unexpected return value from htslib sam_itr_next function '" << ret
          << "' while attempting to read BAM/CRAM file:\n";
      report_state(oss);
      throw blt_exception(oss.str().c_str());
    }
  }

  _is_record_set = (ret >= 0);
  if (_is_record_set) _record_no++;

  return _is_record_set;
}

const char* bam_streamer::target_id_to_name(const int32_t tid) const
{
  // assert(tid < _bfp->header->n_targets);
  if (tid < 0) {
    static const char unmapped[] = "*";
    return unmapped;
  }
  return _hdr->target_name[tid];
}

int32_t bam_streamer::target_name_to_id(const char* seq_name) const
{
  return bam_name2id(_hdr, seq_name);
}

void bam_streamer::report_state(std::ostream& os) const
{
  const bam_record* bamp(get_record_ptr());

  os << "\tbam_stream_label: '" << name() << "'\n";
  if (_is_region && (!_region.empty())) {
    os << "\tbam_stream_selected_region: " << _region << "\n";
  }
  if (nullptr != bamp) {
    os << "\tbam_stream_record_no: " << record_no() << "\n";
    os << "\tbam_record QNAME/read_number: " << bamp->qname() << "/" << bamp->read_no() << "\n";
    const char* chrom_name(target_id_to_name(bamp->target_id()));
    os << "\tbam record RNAME: " << chrom_name << "\n";
    os << "\tbam record POS: " << bamp->pos() << "\n";
  } else {
    os << "\tno bam record currently set\n";
  }
}
