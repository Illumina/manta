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

#include "vcf_record.hpp"
#include "blt_util/parse_util.hpp"

#include <cassert>
#include <cctype>

#include <algorithm>
#include <iostream>

struct convert {
  void operator()(char& c) const { c = toupper((unsigned char)c); }
};

static void stoupper(std::string& s)
{
  std::for_each(s.begin(), s.end(), convert());
}

bool vcf_record::set(const char* s)
{
  static const char     sep('\t');
  static const unsigned maxword(5);

  clear();

  line = s;

  // simple tab parse:
  const char* start(s);
  const char* p(start);

  unsigned wordindex(0);
  while (wordindex < maxword) {
    if ((*p == sep) || (*p == '\n') || (*p == '\0')) {
      switch (wordindex) {
      case 0:
        chrom = std::string(start, p - start);
        break;
      case 1:
        pos = illumina::blt_util::parse_int(start);
        assert(start == p);
        break;
      case 2:
        // skip this field...
        break;
      case 3:
        ref = std::string(start, p - start);
        stoupper(ref);
        break;
      case 4:
        // additional parse loop for ',' character:
        {
          const char* p2(start);

          while (p2 <= p) {
            if ((*p2 == ',') || (p2 == p)) {
              // recognize '.' value and leave alt empty in this case:
              const bool isEmpty((p2 == p) and alt.empty() and ((p2 - start) == 1) and (*start == '.'));
              if (not isEmpty) {
                alt.emplace_back(start, p2 - start);
                stoupper(alt.back());
              }
              start = p2 + 1;
            }
            p2++;
          }
        }
        break;
      default:
        assert(0);
        break;
      }
      start = p + 1;
      wordindex++;
    }
    if ((*p == '\n') || (*p == '\0')) break;
    ++p;
  }

  return (wordindex >= maxword);
}

bool vcf_record::is_normalized() const
{
  // normalized indels are left-aligned, reference-padded, and parsimonious
  // normalized SNVs are a single differing base
  // normalized MNVs (and complex alleles) have differing bases at the beginning
  // and end of the ref and alt alleles.  However, many input VCFs have
  // reference-padded MNVs, which should not affect Strelka's calling, so for
  // now, we're allowing variants to violate left-parsimony in MNVs and complex
  // alleles
  // see http://genome.sph.umich.edu/wiki/Variant_Normalization
  unsigned ref_length = ref.size();
  assert(ref_length != 0);

  for (const auto& alt_allele : alt) {
    unsigned alt_length = alt_allele.size();
    assert(alt_length != 0);

    // all normalized variants with the same length ref and alt
    // must differ at the last ref and alt base.  Any indel that
    // has more than one base in both the ref and the alt must
    // also differ at the last base (this should only happen at
    // complex indels).  This checks for right-padding, i.e. parsimony
    if ((alt_length > 1 && ref_length > 1) || alt_length == ref_length) {
      if ((*alt_allele.rbegin()) == (*ref.rbegin())) {
        return false;
      }
    }

    if (alt_length != ref_length) {
      // this checks that indels are reference-padded
      if ((*alt_allele.begin()) == (*ref.begin())) {
        // this checks that they're left-shifted
        for (unsigned i = ref_length - 1, j = alt_length - 1;; ++i, ++j) {
          if (ref[i] != alt_allele[j]) {
            break;
          } else if (i == 0 || j == 0) {
            return false;
          }
        }
      }
      // if the first and last bases of alleles with two differing lengths
      // do not match, the record represents a complex allele, and fulfills
      // normalization requirements
    }
  }
  return true;
}

std::ostream& operator<<(std::ostream& os, const vcf_record& vcfr)
{
  os << vcfr.chrom << '\t' << vcfr.pos << '\t' << '.' << '\t' << vcfr.ref << '\t';

  const unsigned nalt(vcfr.alt.size());
  for (unsigned a(0); a < nalt; ++a) {
    if (a > 0) os << ',';
    os << vcfr.alt[a];
  }
  os << '\t' << '.' << '\t' << '.' << '\t' << '.' << '\t' << '.' << '\n';

  return os;
}
