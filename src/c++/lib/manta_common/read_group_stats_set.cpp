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
// <https://github.com/downloads/sequencing/licenses/>.
//
/**
 ** Copyright (c) 2007-2010 Illumina, Inc.
 **
 ** This software is covered by the "Illumina Genome Analyzer Software
 ** License Agreement" and the "Illumina Source Code License Agreement",
 ** and certain third party copyright/licenses, and any user of this
 ** source file is bound by the terms therein (see accompanying files
 ** Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
 ** Illumina_Source_Code_License_Agreement.pdf and third party
 ** copyright/license notices).
 **
 ** This file is part of the Guided Reassembly Of Unaligned Paired-End Reads
 ** (GROUPER) software package.
 **
 **/

#include "read_group_stats_set.hh"

#include "blt_util/log.hh"
#include "blt_util/parse_util.hh"
#include "blt_util/string_util.hh"

#include <iostream>



// Stats file data format
const unsigned HEAD_FILL_IDX = 0;
const unsigned HEAD_SRC_IDX  = 1;
const unsigned HEAD_NAME_IDX = 2;

// Stats file data format
const unsigned STAT_SOURCE_IDX             = 0;
const unsigned STAT_INS_SIZE_MEAN_IDX      = 1;
const unsigned STAT_INS_SIZE_SD_IDX        = 2;
const unsigned STAT_INS_SIZE_MEDIAN_IDX    = 3;

const unsigned STAT_QUALITY_MEAN_IDX       = 4;
const unsigned STAT_QUALITY_SD_IDX         = 5;
const unsigned STAT_QUALITY_MEDIAN_IDX     = 6;

const unsigned STAT_SINGLES_MEAN_IDX       = 7;
const unsigned STAT_SINGLES_SD_IDX         = 8;
const unsigned STAT_SINGLES_MEDIAN_IDX     = 9;

const unsigned STAT_READ1_LEN_IDX          = 10;
const unsigned STAT_READ2_LEN_IDX          = 11;
const unsigned STAT_REL_ORIENT_IDX         = 12;



void
read_group_stats_set::
load(std::istream& is) {

    using namespace illumina::blt_util;

    clear();
    std::map<int,std::string> gmap;

    std::string line;
    while (! is.eof()) {
        std::getline(is, line);
        if (line.length() == 0) continue;

        std::vector<std::string> data;
        split_string(line,'\t', data);
        if (((data[0] == "#") && (data[1] == "index"))
            || (data[0] == "# Bam_Size_To_Use")
            || (data[0] == "# Bam_Orig_Path")) {
            continue;
        }

        if (data[0] == "# Bam_Path") {
            gmap[parse_int_str(data[HEAD_SRC_IDX])] = data[HEAD_NAME_IDX];
            continue;
        }
        // Get key string
        const int32_t key = parse_int_str(data[STAT_SOURCE_IDX]);

        // Make sure we have a BAM file source mapping!
        if (gmap.count(key) == 0) {
            log_os << "[ERROR]: While loading stats file unmapped data found for index: " << key
                   << ", on line: " << line << '\n';
            exit(EXIT_FAILURE);
        }

        const read_group_stats rps(data);
        set_stats(gmap[key],rps);
    }
}



void
read_group_stats_set::
store(std::ostream& os) const {
    const unsigned n_groups(_group.size());
    for (unsigned i(0); i<n_groups; ++i)
    {
        os << "# Bam_Path\t" << i << "\t" << _group.get_key(i) << '\n';
    }
    // write column header for better readability
    os << "#\tindex"
       << "\tmeanInsertSize\tsdInsSize\tmedianInsSize"
       << "\tmeanQuality\tsdQuality\tmedianQuality"
       << "\tmeanSingleAln\tsdSingleAln\tmedianSingleAln"
       << "\treadLen1\treadLen2\treadOrientation"
       << '\n';

    for (unsigned i(0); i<n_groups; ++i)
    {
        os << i << '\t';
        get_stats(i).store(os);
        os << '\n';
    }
}

