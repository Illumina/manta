// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#include "htsapi/bam_header_util.hh"
#include "blt_util/blt_exception.hh"
#include "blt_util/parse_util.hh"
#include "blt_util/string_util.hh"

#include <sstream>

#include "boost/tokenizer.hpp"



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
        oss << "ERROR: can't parse bam_region [err 1] " << region << "\n";
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
        oss << "ERROR: can't parse bam_region [err 2] " << region << "\n"
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
        oss << "ERROR: can't parse bam_region [err 3] " << region << "\n";
        throw blt_exception(oss.str().c_str());
    }

    begin_pos = (illumina::blt_util::parse_int_str(words2[0]))-1;

    if (1 == words2.size()) return;
    end_pos = (illumina::blt_util::parse_int_str(words2[1]));
}



bool
check_header_compatibility(
    const bam_header_t* h1,
    const bam_header_t* h2)
{
    if (h1->n_targets != h2->n_targets)
    {
        return false;
    }

    for (int32_t i(0); i<h1->n_targets; ++i)
    {
        if (h1->target_len[i] != h2->target_len[i]) return false;
        if (0 != strcmp(h1->target_name[i],h2->target_name[i])) return false;
    }
    return true;
}



std::string
get_bam_header_sample_name(
    const std::string& bam_header_text,
    const char* default_sample_name)
{
    using namespace boost;
    char_separator<char> sep("\t\n");
    tokenizer< char_separator<char>> tokens(bam_header_text, sep);
    for (const auto& t : tokens)
    {
        const auto sm = t.find("SM:");
        if (std::string::npos == sm) continue;
        return t.substr(sm+3);
    }

    return default_sample_name;
}
