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

#include "AssembledContig.hpp"
#include "blt_util/seq_printer.hpp"

#include <iostream>

std::ostream& operator<<(std::ostream& os, const AssembledContig& contig)
{
  os << "CONTIG size: " << contig.seq.size() << " seedCount: " << contig.seedReadCount
     << " supportReads: " << contig.supportReads.size() << " seq:\n";
  printSeq(contig.seq, os);
  os << "\n";

  return os;
}
