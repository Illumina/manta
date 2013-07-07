
///
/// \author Bret Barnes
///

#pragma once

#include <map>
#include <string>


/// retrieve a map of chromosome sizes from the fasta index
///
void
get_chrom_sizes(
        const std::string& fai_file,
        std::map<std::string,unsigned>& chrom_sizes);


/// retrieve size of specific chromosome from the fasta index
///
unsigned
get_chrom_length(
        const std::string& fai_file,
        const std::string& chrom_name);


/// get reference sequence from region
void
get_region_seq(
        const std::string& ref_file,
        const std::string& fa_region,
        std::string& ref_seq);

/// get reference sequence from decomposed region
void
get_region_seq(
        const std::string& ref_file,
        const std::string& chrom,
        const int begin_pos,
        const int end_pos,
        std::string& ref_seq);


/// get reference sequence from decomposed region and run standardization result
///
/// \param begin_pos begin position (zero-indexed)
/// \param end_pos end position (zero-indexed)
void
get_standardized_region_seq(
        const std::string& ref_file,
        const std::string& chrom,
        const int begin_pos,
        const int end_pos,
        std::string& ref_seq);
