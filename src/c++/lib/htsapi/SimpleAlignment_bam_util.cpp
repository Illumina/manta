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

///
///
///

#include "SimpleAlignment_bam_util.hpp"
#include "htsapi/align_path_bam_util.hpp"

void getAlignment(const bam_record& bamRead, SimpleAlignment& al)
{
  al.is_fwd_strand = bamRead.is_fwd_strand();
  al.tid           = bamRead.target_id();
  al.pos           = (bamRead.pos() - 1);

  bam_cigar_to_apath(bamRead.raw_cigar(), bamRead.n_cigar(), al.path);
}

SimpleAlignment getAlignment(const bam_record& bamRead)
{
  SimpleAlignment al;
  getAlignment(bamRead, al);
  return al;
}

SimpleAlignment getKnownOrFakedMateAlignment(const bam_record& bamRead)
{
  using namespace ALIGNPATH;

  SimpleAlignment al;
  assert(!bamRead.is_mate_unmapped());
  al.is_fwd_strand = bamRead.is_mate_fwd_strand();
  al.tid           = bamRead.mate_target_id();
  al.pos           = (bamRead.mate_pos() - 1);

  static const char mateCigarTag[] = {'M', 'C'};
  const char*       mateCigarString(bamRead.get_string_tag(mateCigarTag));
  if (nullptr != mateCigarString) {
    cigar_to_apath(mateCigarString, al.path);
  } else {
    al.path.emplace_back(MATCH, bamRead.read_size());
  }
  return al;
}
