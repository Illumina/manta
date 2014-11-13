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

/// random vcf utilities
///
/// \author Chris Saunders
///

#pragma once

#include <cstring>
#include <iosfwd>
#include <vector>


namespace VCFID
{
enum index_t
{
    CHROM,
    POS,
    ID,
    REF,
    ALT,
    QUAL,
    FILT,
    INFO,
    FORMAT,
    SAMPLE,
    SIZE
};
}



inline
const char*
vcf_col_label()
{
    static const char h[] = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
    return h;
}


std::ostream&
vcf_fileDate(std::ostream& os);


void
write_vcf_filter(
    std::ostream& os,
    const char* id,
    const char* desc);


// look for 'key' in vcf FORMAT field, provide index of key or return
// false
//
inline
bool
get_format_key_index(const char* format,
                     const char* key,
                     unsigned& index)
{
    index=0;
    do
    {
        if (index) format++;
        if (0==strncmp(format,key,strlen(key))) return true;
        index++;
    }
    while (NULL != (format=strchr(format,':')));
    return false;
}



// return pointer to
//
inline
const char*
get_format_string_nocopy(const char* const* word,
                         const char* key)
{
    unsigned keynum(0);
    if (! get_format_key_index(word[VCFID::FORMAT],key,keynum)) return NULL;

    const char* sample(word[VCFID::SAMPLE]);
    for (; keynum; sample++)
    {
        if (! *sample) return NULL;
        if ((*sample)==':') keynum--;
    }
    return sample;
}



// returns -1 for '.' alleles
void
parse_gt(const char* gt,
         std::vector<int>& gti,
         const bool is_allow_bad_end_char=false);
