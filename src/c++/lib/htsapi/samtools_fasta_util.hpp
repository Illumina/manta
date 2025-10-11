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
/// \author Bret Barnes
///

#pragma once

#include <map>
#include <string>

/// retrieve a map of chromosome sizes from the fasta index
///
void get_chrom_sizes(const std::string& fai_file, std::map<std::string, unsigned>& chrom_sizes);

/// retrieve size of specific chromosome from the fasta index
///
unsigned get_chrom_length(const std::string& fai_file, const std::string& chrom_name);

/// get reference sequence from region
void get_region_seq(const std::string& ref_file, const std::string& fa_region, std::string& ref_seq);

/// get reference sequence from decomposed region
void get_region_seq(
    const std::string& ref_file,
    const std::string& chrom,
    const int          begin_pos,
    const int          end_pos,
    std::string&       ref_seq);

/// get reference sequence from decomposed region and run standardization result
///
/// \param begin_pos begin position (zero-indexed, closed)
/// \param end_pos end position (zero-indexed, closed)
void get_standardized_region_seq(
    const std::string& ref_file,
    const std::string& chrom,
    const int          begin_pos,
    const int          end_pos,
    std::string&       ref_seq);
