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

#include "vcf_streamer.hpp"

#include "blt_util/log.hpp"
#include "blt_util/seq_util.hpp"

#include "common/Exceptions.hpp"

#include <sys/stat.h>
#include <cassert>

#include <iostream>
#include <set>
#include <sstream>
#include <string>

/// return true only if all chromosomes in the bcf/vcf exist in the
/// bam header
static void check_bam_bcf_header_compatibility(
    const char* bcf_filename, const bcf_hdr_t* bcf_header, const bam_hdr_t& bam_header)
{
  assert(bcf_header);

  // build set of chrom labels from BAM:
  std::set<std::string> bamLabels;
  for (int32_t i(0); i < bam_header.n_targets; ++i) {
    bamLabels.insert(std::string(bam_header.target_name[i]));
  }
  int labelCount(0);

  const char** bcfLabels = bcf_hdr_seqnames(bcf_header, &labelCount);

  for (int labelIndex(0); labelIndex < labelCount; ++labelIndex) {
    if (bamLabels.find(std::string(bcfLabels[labelIndex])) != bamLabels.end()) continue;
    std::ostringstream oss;
    oss << "Chromosome label '" << bcfLabels[labelIndex] << "' in BCF/VCF file '" << bcf_filename
        << "' does not exist in the BAM header";
    BOOST_THROW_EXCEPTION(illumina::common::GeneralException(oss.str()));
  }

  free(bcfLabels);
}

vcf_streamer::vcf_streamer(const char* filename, const char* region, const bool requireNormalized)
  : hts_streamer(filename, region), _hdr(nullptr), _requireNormalized(requireNormalized)
{
  //
  // note with the switch to samtools 1.X vcf/bcf still involve predominantly separate
  // apis -- no bcf support added here, but a shared function has been chosen where possible
  // ... for_instance hts_open/bcf_hdr_read should work with both vcf and bcf
  //

  _hdr = bcf_hdr_read(_hfp);
  if (!_hdr) {
    std::ostringstream oss;
    oss << "Failed to load header for VCF file: '" << filename << "'";
    BOOST_THROW_EXCEPTION(illumina::common::GeneralException(oss.str()));
  }
  _sampleCount = bcf_hdr_nsamples(_hdr);
}

vcf_streamer::~vcf_streamer()
{
  if (_hdr) bcf_hdr_destroy(_hdr);
}

bool vcf_streamer::next()
{
  if (_is_stream_end || (!_hfp) || (!_titr)) return false;

  while (true) {
    if (tbx_itr_next(_hfp, _tidx, _titr, &_kstr) < 0) {
      _is_stream_end = true;
    } else {
      _is_stream_end = (nullptr == _kstr.s);
    }
    _is_record_set = (!_is_stream_end);
    if (!_is_record_set) break;

    // filter out header for whole file access case:
    if (_kstr.s[0] == '#') continue;

    _record_no++;

    if (!_vcfrec.set(_kstr.s)) {
      std::ostringstream oss;
      oss << "Can't parse VCF record: '" << _kstr.s << "'";
      report_state(oss);
      BOOST_THROW_EXCEPTION(illumina::common::GeneralException(oss.str()));
    }

    if (_vcfrec.isSimpleVariantLocus()) {
      if (!_vcfrec.is_normalized()) {
        std::ostringstream oss;
        oss << "Input VCF record is not normalized:\n";
        report_state(oss);

        if (_requireNormalized) {
          oss << "Please normalize all records in this VCF with a tool such as vt, then resubmit\n";
          BOOST_THROW_EXCEPTION(illumina::common::GeneralException(oss.str()));
        } else {
          log_os << "WARNING: " << oss.str();
          log_os << "All alleles in this VCF record have been skipped.\n";
          continue;
        }
      }
    }

    break;  // found expected vcf record type
  }

  return _is_record_set;
}

void vcf_streamer::report_state(std::ostream& os) const
{
  const vcf_record* vcfp(get_record_ptr());

  os << "\tvcf_stream_label: " << name() << "\n";
  if (vcfp) {
    os << "\tvcf_stream_record_no: " << record_no() << "\n"
       << "\tvcf_record: " << *(vcfp) << "\n";
  } else {
    os << "\tno vcf record currently set\n";
  }
}

void vcf_streamer::validateBamHeaderChromSync(const bam_hdr_t& header) const
{
  check_bam_bcf_header_compatibility(name(), _hdr, header);
}
