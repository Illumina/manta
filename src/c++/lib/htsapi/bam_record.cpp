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

/// \file

/// \author Chris Saunders
///

#include "blt_util/blt_exception.hh"
#include "htsapi/align_path_bam_util.hh"
#include "htsapi/bam_record.hh"

#include <iostream>
#include <sstream>



std::ostream&
operator<<(std::ostream& os, const bam_record& br)
{
    if (br.empty())
    {
        os << "NONE";
    }
    else
    {
        os << br.qname() << "/" << br.read_no()
           << " tid:pos:strand " << br.target_id() << ":" << (br.pos()-1) << ":" << (br.is_fwd_strand() ? '+' : '-');

        ALIGNPATH::path_t apath;
        bam_cigar_to_apath(br.raw_cigar(),br.n_cigar(),apath);
        os << " cigar: " << apath;

        /// print SAtag if present:
        static const char satag[] = {'S','A'};
        const char* saStr(br.get_string_tag(satag));
        if (NULL != saStr)
        {
            os  << " sa: " << saStr;
        }

        if (br.is_paired())
        {
            os  << " mate_tid:pos:strand " << br.mate_target_id() << ":" << (br.mate_pos()-1) << ":" << (br.is_mate_fwd_strand() ? '+' : '-');
        }

    }
    return os;
}



unsigned
bam_record::
alt_map_qual(const char* tag) const
{
    uint8_t* alt_ptr(bam_aux_get(_bp,tag));
    if ((NULL != alt_ptr) && is_int_code(alt_ptr[0]))
    {
        const int alt_map(bam_aux2i(alt_ptr));
        if (alt_map<0)
        {
            std::ostringstream oss;
            oss << "ERROR: Unexpected negative value in optional BAM tag: '" << std::string(tag,2) << "'\n"
                << "\tRemove the --eland-compatibility flag to stop using this tag.\n";
            throw blt_exception(oss.str().c_str());
        }
        return static_cast<unsigned>(alt_map);
    }
    else
    {
        return map_qual();
    }
}



const char*
bam_record::
get_string_tag(const char* tag) const
{

    // retrieve the BAM tag
    uint8_t* pTag = bam_aux_get(_bp, tag);
    if (!pTag) return NULL;

    // skip tags that are not encoded as a null-terminated string
    if (pTag[0] != 'Z') return NULL;
    ++pTag;

    return (const char*)pTag;
}



bool
bam_record::
get_num_tag(const char* tag, int32_t& num) const
{

    // retrieve the BAM tag
    uint8_t* pTag = bam_aux_get(_bp, tag);
    if (!pTag) return false;

    // skip tags that are not encoded as integers
    if (!is_int_code(pTag[0])) return false;
    num = bam_aux2i(pTag);

    return true;
}
