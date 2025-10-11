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
/// \author Chris Saunders
///

#include "htsapi/align_path_bam_util.hpp"

#include <cassert>

using namespace ALIGNPATH;

void bam_cigar_to_apath(const uint32_t* bam_cigar, const unsigned n_cigar, path_t& apath)
{
  // this assertion isn't really required...
  //    assert(n_cigar>0);
  apath.resize(n_cigar);
  for (unsigned i(0); i < n_cigar; ++i) {
    apath[i].length = (bam_cigar[i] >> BAM_CIGAR_SHIFT);
    apath[i].type   = static_cast<align_t>(1 + (bam_cigar[i] & BAM_CIGAR_MASK));
  }
}

void apath_to_bam_cigar(const path_t& apath, uint32_t* bam_cigar)
{
  const unsigned as(apath.size());
  for (unsigned i(0); i < as; ++i) {
    const path_segment& ps(apath[i]);
    assert(ps.type != NONE);
    bam_cigar[i] = (ps.length << BAM_CIGAR_SHIFT | (static_cast<uint32_t>(ps.type) - 1));
  }
}

void edit_bam_cigar(const path_t& apath, bam1_t& br)
{
  bam1_core_t& bc(br.core);

  const int old_n_cigar(bc.n_cigar);
  const int new_n_cigar(apath.size());
  const int delta(4 * (new_n_cigar - old_n_cigar));

  if (0 != delta) {
    const int end(bc.l_qname + (4 * old_n_cigar));
    change_bam_data_segment_len(end, delta, br);
    bc.n_cigar = new_n_cigar;
  }

  // update content of cigar array:
  apath_to_bam_cigar(apath, bam_get_cigar(&br));
}
