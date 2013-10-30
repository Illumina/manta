// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

/// \author Chris Saunders
///

#include "blt_util/bam_header_util.hh"
#include "blt_util/blt_exception.hh"
#include "blt_util/parse_util.hh"
#include "blt_util/string_util.hh"

#include <sstream>



void
parse_bam_region(
    const bam_header_info& header,
    const std::string& region,
    int32_t& tid,
    int32_t& begin_pos,
    int32_t& end_pos)
{
    static const char region_sep1(':');
    std::vector<std::string> words;
    split_string(region,region_sep1,words);

    if (words.empty() || words[0].empty() || (words.size() > 2))
    {
        std::ostringstream oss;
        oss << "ERROR: can't parse bam_region " << region << "\n";
        throw blt_exception(oss.str().c_str());
    }

    bool isFound(false);
    const unsigned n_chroms(header.chrom_data.size());
    for (unsigned i(0); i<n_chroms; ++i)
    {
        if (words[0]==header.chrom_data[i].label)
        {
            tid=i;
            isFound=true;
            break;
        }
    }

    if (! isFound)
    {
        std::ostringstream oss;
        oss << "ERROR: can't parse bam_region " << region << "\n"
            << "\tchromosome: '" << words[0] << "' not found in header\n";
        throw blt_exception(oss.str().c_str());
    }

    begin_pos = 0;
    end_pos = header.chrom_data[tid].length;
    if (1 == words.size()) return;

    static const char region_sep2('-');
    std::vector<std::string> words2;
    split_string(words[1],region_sep2,words2);

    if (words2.empty() || (words2.size() > 2))
    {
        std::ostringstream oss;
        oss << "ERROR: can't parse bam_region " << region << "\n";
        throw blt_exception(oss.str().c_str());
    }

    begin_pos = (illumina::blt_util::parse_int_str(words2[0]))-1;

    if (1 == words2.size()) return;
    end_pos = (illumina::blt_util::parse_int_str(words2[1]));
}

/*****************************************************************************/

void
make_chrom_tid_map(
    const bam_header_info& header,
    std::map<std::string, int32_t>& chromNameTidMap)
{
    const unsigned n_chroms(header.chrom_data.size());

    for (unsigned i(0); i<n_chroms; ++i)
    {
        chromNameTidMap[header.chrom_data[i].label] = i;
    }
}

/*****************************************************************************/
