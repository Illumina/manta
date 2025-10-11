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

#include "blt_util/chrom_depth_map.hpp"
#include "blt_util/blt_exception.hpp"
#include "blt_util/log.hpp"
#include "blt_util/parse_util.hpp"

#include <cstdlib>
#include <cstring>

#include <fstream>
#include <iostream>

// parse the chrom depth file
void parse_chrom_depth(const std::string& chrom_depth_file, cdmap_t& chrom_depth)
{
  if (chrom_depth_file.empty()) return;

  std::ifstream depth_is(chrom_depth_file.c_str());
  if (!depth_is) {
    log_os << "ERROR: Failed to open chrom depth file '" << chrom_depth_file << "'\n";
    exit(EXIT_FAILURE);
  }

  static const unsigned buff_size(1024);
  char                  buff[buff_size];

  unsigned line_no(0);

  while (true) {
    depth_is.getline(buff, buff_size);
    if (!depth_is) {
      if (depth_is.eof())
        break;
      else {
        log_os << "ERROR: unexpected failure while attempting to read chrom depth file line " << (line_no + 1)
               << "\n";
        exit(EXIT_FAILURE);
      }
    } else {
      ++line_no;
    }

    char* word2(strchr(buff, '\t'));
    if (nullptr == word2) {
      log_os << "ERROR: unexpected format in read chrom depth file line " << (line_no) << "\n";
      exit(EXIT_FAILURE);
    }
    *(word2++) = '\0';
    try {
      const char* s(word2);
      chrom_depth[buff] = illumina::blt_util::parse_double(s);
    } catch (const blt_exception&) {
      log_os << "ERROR: unexpected format in read chrom depth file line " << (line_no) << "\n";
      throw;
    }
    if (chrom_depth[buff] < 0) {
      log_os << "ERROR: Chromosome depth estimate is negative. Chromosome: '" << buff
             << "' Depth: " << chrom_depth[buff] << "\n";
      exit(EXIT_FAILURE);
    }
  }
}
