//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2025 Illumina, Inc.
//
// This program is licensed under the terms of the Polyform strict license
//
// ***As far as the law allows, the software comes as is, without
// any warranty or condition, and the licensor will not be liable
// to you for any damages arising out of these terms or the use
// or nature of the software, under any kind of legal claim.***
//
// You should have received a copy of the PolyForm Strict License 1.0.0
// along with this program.  If not, see <https://polyformproject.org/licenses/strict/1.0.0>.
//
//

/// \file
///
/// \brief VCF utilities
///
/// \author Chris Saunders

#pragma once

#include <cstring>
#include <iosfwd>
#include <vector>

namespace VCFID {
enum index_t { CHROM, POS, ID, REF, ALT, QUAL, FILT, INFO, FORMAT, SAMPLE, SIZE };
}

inline const char* vcf_col_label()
{
  static const char h[] = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
  return h;
}

std::ostream& vcf_fileDate(std::ostream& os);

void write_vcf_filter(std::ostream& os, const char* id, const char* desc);

/// look for 'key' in vcf FORMAT field, provide index of key or return
/// false
///
inline bool get_format_key_index(const char* format, const char* key, unsigned& index)
{
  index = 0;
  do {
    if (index > 0) format++;
    if (0 == strncmp(format, key, strlen(key))) return true;
    index++;
  } while (nullptr != (format = strchr(format, ':')));
  return false;
}

// return pointer to
//
inline const char* get_format_string_nocopy(const char* const* word, const char* key)
{
  unsigned keynum(0);
  if (!get_format_key_index(word[VCFID::FORMAT], key, keynum)) return nullptr;

  const char* sample(word[VCFID::SAMPLE]);
  for (; keynum > 0; sample++) {
    if (!*sample) return nullptr;
    if ((*sample) == ':') keynum--;
  }
  return sample;
}

/// returns -1 for '.' alleles
void parse_gt(const char* gt, std::vector<int>& gti, const bool is_allow_bad_end_char = false);
