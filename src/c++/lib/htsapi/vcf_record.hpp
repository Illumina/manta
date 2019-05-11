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

#pragma once

#include "blt_util/seq_util.hpp"

#include <iosfwd>
#include <string>
#include <vector>

struct vcf_record {
  vcf_record() { clear(); }

  /// set the vcf record from record string s, return false on error
  bool set(const char* s);

  void clear()
  {
    chrom.clear();
    pos = 0;
    ref.clear();
    alt.clear();
    line = nullptr;
  }

  /// test if the variant represents a "simple" SNV or small indel
  /// where all REF/ALT entries are composed of only [ACGTN] sequences
  //
  // N.B. the current implementation of is_valid does not
  // allow symbolic ALTs, e.g. representing a deletion using
  // the <DEL> symbolic allele and the END INFO field.  This
  // is probably fine, but worth noting
  bool isSimpleVariantLocus() const
  {
    if (ref.empty()) return false;
    if (alt.empty()) return false;
    if (!is_valid_seq(ref.c_str())) return false;
    for (const auto& alt_allele : alt) {
      if (!is_valid_seq(alt_allele.c_str())) return false;
    }
    return true;
  }

  /// check for REF or ALT alleles with a size > 1, alleles with equal length REF and ALT sequences will
  /// count as an indel
  bool is_indel() const
  {
    if (!isSimpleVariantLocus()) return false;
    if ((ref.size() > 1) && (!alt.empty())) return true;
    for (const auto& alt_allele : alt) {
      if (alt_allele.size() > 1) return true;
    }
    return false;
  }

  bool is_snv() const
  {
    if (!isSimpleVariantLocus()) return false;
    if (1 != ref.size()) return false;
    for (const auto& alt_allele : alt) {
      if (1 != alt_allele.size()) return false;
    }
    return true;
  }

  /// complements is_snv() by taking case of the form: REF="A", ALT="."
  bool is_ref_site() const
  {
    if (ref.empty()) return false;
    if (!is_valid_seq(ref.c_str())) return false;
    if (1 != ref.size()) return false;
    return (alt.empty());
  }

  bool is_normalized() const;

  std::string              chrom;
  int                      pos = 0;
  std::string              ref;
  std::vector<std::string> alt;
  const char*              line = nullptr;
};

std::ostream& operator<<(std::ostream& os, const vcf_record& vcfr);
