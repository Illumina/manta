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
/// \author Ole Schulz-Trieglaff
///

#pragma once

#include <iosfwd>
#include <set>
#include <string>
#include <vector>

#include "blt_util/known_pos_range2.hpp"

/// \brief data pertaining to a de-novo assembly contig
///
/// stores for each contig the sequence and the number of reads
/// containing its seeding k-mer
///
struct AssembledContig {
  std::string seq;  ///< contigsequence

  // reads used for assembly of contig <read_no,mapping position to contig>
  //std::map<std::string,int> contigReads;

  unsigned seedReadCount = 0;  ///< no of reads containing the seeding kmer

  std::set<unsigned> supportReads;
  std::set<unsigned> rejectReads;

  known_pos_range2 conservativeRange;  ///< subsection of the contig with conservative coverage
};

std::ostream& operator<<(std::ostream& os, const AssembledContig& contig);

typedef std::vector<AssembledContig> Assembly;
