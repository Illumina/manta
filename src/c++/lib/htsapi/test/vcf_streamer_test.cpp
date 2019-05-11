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

#include "testConfig.h"

#include "htsapi/vcf_streamer.hpp"

#include "common/Exceptions.hpp"

#include "boost/test/unit_test.hpp"

BOOST_AUTO_TEST_SUITE(test_vcf_streamer)

static const char* getTestpath()
{
  static const std::string testPath(std::string(TEST_DATA_PATH) + "/vcf_streamer_test.vcf.gz");
  return testPath.c_str();
}

BOOST_AUTO_TEST_CASE(test_vcf_streamer_region)
{
  vcf_streamer vcfs(getTestpath(), "chr1:750000-822000");

  const vcf_record* vptr(nullptr);

  // first check that valid record does not throw an exception
  BOOST_REQUIRE_NO_THROW(vcfs.next());
  vptr = vcfs.get_record_ptr();
  assert(vptr != nullptr);

  // VCF record test for
  // chr1    757807  MantaDEL:44:0:0:0:1:0   CCCTGGCCAGCAGATCCACCCTGTCTATACTACCTG    C       ...
  // testing variant assignment, normalize and validity checks, position, and ref
  BOOST_REQUIRE(vptr->isSimpleVariantLocus());
  BOOST_REQUIRE(vptr->is_indel());
  BOOST_REQUIRE(vptr->is_normalized());
  BOOST_REQUIRE(!vptr->is_snv());

  BOOST_REQUIRE_EQUAL(vptr->pos, 757807);
  BOOST_REQUIRE_EQUAL(vptr->ref, "CCCTGGCCAGCAGATCCACCCTGTCTATACTACCTG");

  // also check that a valid record returns true
  BOOST_REQUIRE(vcfs.next());
  vptr = vcfs.get_record_ptr();
  assert(vptr != nullptr);

  // VCF record test for
  // chr1    758807  FAKED   C       T,A     ...
  // testing variant assignment, normalize and validity checks, position, alt size, and alt
  BOOST_REQUIRE(vptr->isSimpleVariantLocus());
  BOOST_REQUIRE(!vptr->is_indel());
  BOOST_REQUIRE(vptr->is_snv());
  BOOST_REQUIRE(vptr->is_normalized());

  BOOST_REQUIRE_EQUAL(vptr->pos, 758807);
  BOOST_REQUIRE_EQUAL(vptr->alt.size(), 2u);
  BOOST_REQUIRE_EQUAL(vptr->alt[0], "T");

  BOOST_REQUIRE(vcfs.next());
  vptr = vcfs.get_record_ptr();
  assert(vptr != nullptr);

  // VCF record test for
  // chr1    821604  MantaINS:53:0:0:0:1:0   T       TGCCCTTTGGCAGAGCAGGTGTGCTGTGCTG ...
  // testing variant assignment, normalize and validity checks, position, alt size, and alt
  BOOST_REQUIRE(vptr->isSimpleVariantLocus());
  BOOST_REQUIRE(vptr->is_indel());
  BOOST_REQUIRE(!vptr->is_snv());
  BOOST_REQUIRE(vptr->is_normalized());

  BOOST_REQUIRE_EQUAL(vptr->pos, 821604);
  BOOST_REQUIRE_EQUAL(vptr->alt.size(), 1u);
  BOOST_REQUIRE_EQUAL(vptr->alt[0], "TGCCCTTTGGCAGAGCAGGTGTGCTGTGCTG");

  vcfs.resetRegion("chr10:89717700-89717810");

  // tests that we can reset VCF regions
  BOOST_REQUIRE(vcfs.next());
  vptr = vcfs.get_record_ptr();
  assert(vptr != nullptr);

  // VCF record test for
  // chr10   89717769        COSM30622       TA      T       ...
  // testing that deletion is reported as normalized
  BOOST_REQUIRE(vptr->isSimpleVariantLocus());
  BOOST_REQUIRE(vptr->is_indel());
  BOOST_REQUIRE(!vptr->is_snv());
  BOOST_REQUIRE(vptr->is_normalized());
  BOOST_REQUIRE_EQUAL(vptr->pos, 89717769);
  BOOST_REQUIRE_EQUAL(vptr->alt.size(), 1u);
  BOOST_REQUIRE_EQUAL(vptr->alt[0], "T");

  BOOST_REQUIRE_THROW(vcfs.next(), illumina::common::GeneralException);
  vptr = vcfs.get_record_ptr();
  assert(vptr != nullptr);

  // VCF record test for
  // chr10   89717774        COSM5809        AA      A       ...
  // testing that deletion is not reported as left-shifted
  BOOST_REQUIRE(vptr->isSimpleVariantLocus());
  BOOST_REQUIRE(vptr->is_indel());
  BOOST_REQUIRE(!vptr->is_snv());
  BOOST_REQUIRE(!vptr->is_normalized());
  BOOST_REQUIRE_EQUAL(vptr->pos, 89717774);
  BOOST_REQUIRE_EQUAL(vptr->alt.size(), 1u);
  BOOST_REQUIRE_EQUAL(vptr->alt[0], "A");

  BOOST_REQUIRE_THROW(vcfs.next(), illumina::common::GeneralException);
  vptr = vcfs.get_record_ptr();
  assert(vptr != nullptr);

  // VCF record test for
  // chr10   89717775        COSM5012        A       AA      ...
  // testing that insertion is not reported as left-shifted
  BOOST_REQUIRE(vptr->isSimpleVariantLocus());
  BOOST_REQUIRE(vptr->is_indel());
  BOOST_REQUIRE(!vptr->is_snv());
  BOOST_REQUIRE(!vptr->is_normalized());
  BOOST_REQUIRE_EQUAL(vptr->pos, 89717775);
  BOOST_REQUIRE_EQUAL(vptr->alt.size(), 1u);
  BOOST_REQUIRE_EQUAL(vptr->alt[0], "AA");

  BOOST_REQUIRE(vcfs.next());
  vptr = vcfs.get_record_ptr();
  assert(vptr != nullptr);

  // VCF record test for
  // chr10   89717776        VALID_MNV_1     TCG     AGT     ...
  // testing that MNV is reported as normalized
  BOOST_REQUIRE(vptr->isSimpleVariantLocus());
  BOOST_REQUIRE(vptr->is_indel());
  BOOST_REQUIRE(!vptr->is_snv());
  BOOST_REQUIRE(vptr->is_normalized());
  BOOST_REQUIRE_EQUAL(vptr->pos, 89717776);
  BOOST_REQUIRE_EQUAL(vptr->alt.size(), 1u);
  BOOST_REQUIRE_EQUAL(vptr->alt[0], "AGT");

  BOOST_REQUIRE_THROW(vcfs.next(), illumina::common::GeneralException);
  vptr = vcfs.get_record_ptr();
  assert(vptr != nullptr);

  // VCF record test for
  // chr10   89717779        INVALID_MNV_1   TCG     AGG     ...
  // testing that MNV is not reported as normalized
  BOOST_REQUIRE(vptr->isSimpleVariantLocus());
  BOOST_REQUIRE(vptr->is_indel());
  BOOST_REQUIRE(!vptr->is_snv());
  BOOST_REQUIRE(!vptr->is_normalized());
  BOOST_REQUIRE_EQUAL(vptr->pos, 89717779);
  BOOST_REQUIRE_EQUAL(vptr->alt.size(), 1u);
  BOOST_REQUIRE_EQUAL(vptr->alt[0], "AGG");

  BOOST_REQUIRE(vcfs.next());
  vptr = vcfs.get_record_ptr();
  assert(vptr != nullptr);

  // VCF record test for
  // chr10   89717782        INVALID_MNV_2   ACG     AGT     ...
  // testing that MNV is reported as normalized, despite first base
  // matching (Manta returns MNV candidates with reference-padding, so
  // this is to maintain consistency with its inputs)
  BOOST_REQUIRE(vptr->isSimpleVariantLocus());
  BOOST_REQUIRE(vptr->is_indel());
  BOOST_REQUIRE(!vptr->is_snv());
  BOOST_REQUIRE(vptr->is_normalized());
  BOOST_REQUIRE_EQUAL(vptr->pos, 89717782);
  BOOST_REQUIRE_EQUAL(vptr->alt.size(), 1u);
  BOOST_REQUIRE_EQUAL(vptr->alt[0], "AGT");

  BOOST_REQUIRE_THROW(vcfs.next(), illumina::common::GeneralException);
  vptr = vcfs.get_record_ptr();
  assert(vptr != nullptr);

  // VCF record test for
  // chr10   89717785        INVALID_SNV   A     A     ...
  // testing that SNV is not reported as normalized
  BOOST_REQUIRE(vptr->isSimpleVariantLocus());
  BOOST_REQUIRE(!vptr->is_indel());
  BOOST_REQUIRE(vptr->is_snv());
  BOOST_REQUIRE(!vptr->is_normalized());
  BOOST_REQUIRE_EQUAL(vptr->pos, 89717785);
  BOOST_REQUIRE_EQUAL(vptr->alt.size(), 1u);
  BOOST_REQUIRE_EQUAL(vptr->alt[0], "A");

  BOOST_REQUIRE_THROW(vcfs.next(), illumina::common::GeneralException);
  vptr = vcfs.get_record_ptr();
  assert(vptr != nullptr);

  // VCF record test for
  // chr10   89717785        RIGHT_PAD_INDEL ACCC    AC      ...
  // testing that indel is not reported as normalized
  BOOST_REQUIRE(vptr->isSimpleVariantLocus());
  BOOST_REQUIRE(vptr->is_indel());
  BOOST_REQUIRE(!vptr->is_snv());
  BOOST_REQUIRE(!vptr->is_normalized());
  BOOST_REQUIRE_EQUAL(vptr->pos, 89717785);
  BOOST_REQUIRE_EQUAL(vptr->alt.size(), 1u);
  BOOST_REQUIRE_EQUAL(vptr->alt[0], "AC");

  BOOST_REQUIRE(vcfs.next());
  vptr = vcfs.get_record_ptr();
  assert(vptr != nullptr);

  // VCF record test for
  // chr10   89717790        COMPLEX GAGCTGTG   AGCT ...
  // testing that complex allele is reported as normalized
  BOOST_REQUIRE(vptr->isSimpleVariantLocus());
  BOOST_REQUIRE(vptr->is_indel());
  BOOST_REQUIRE(!vptr->is_snv());
  BOOST_REQUIRE(vptr->is_normalized());
  BOOST_REQUIRE_EQUAL(vptr->pos, 89717790);
  BOOST_REQUIRE_EQUAL(vptr->alt.size(), 1u);
  BOOST_REQUIRE_EQUAL(vptr->alt[0], "AGCT");

  BOOST_REQUIRE(vcfs.next());
  vptr = vcfs.get_record_ptr();
  assert(vptr != nullptr);

  // VCF record test for complex val
  // chr10   89717799    INVALID_ALT A   <DEL>   .   .   .
  // this should not be parsed as an SNV or indel
  BOOST_REQUIRE(not vptr->isSimpleVariantLocus());

  // testing that next returns false after last record
  BOOST_REQUIRE(!vcfs.next());
  vptr = vcfs.get_record_ptr();
  assert(vptr == nullptr);
}

BOOST_AUTO_TEST_SUITE_END()
