//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2018 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

/// \file
/// \author Trevor Ramsay
///

#include "testAlignmentDataUtil.hh"

#include "blt_util/align_path.hh"
#include "htsapi/align_path_bam_util.hh"



bam_header_info
buildTestBamHeader()
{
    bam_header_info bamHeader;
    bamHeader.chrom_data.emplace_back("chrM",1000000);
    bamHeader.chrom_data.emplace_back("chrT",1000000);

    int32_t chromIndex(0); 
    for (const auto& chromData : bamHeader.chrom_data)
    {
        bamHeader.chrom_to_index.insert(std::make_pair(chromData.label, chromIndex));
        chromIndex++;
    }
    return bamHeader;
}



void
buildTestBamRecord(
    bam_record& bamRead,
    int targetID,
    int pos,
    int mateTargetID,
    int matePos,
    int fragmentSize,
    int mapQ,
    std::string cigarString,
    std::string querySeq)
{
    bam1_t& bamData(*(bamRead.get_data()));

    // set qname
    {
        edit_bam_qname("buildTestBamRecord", bamData);
    }

    // set CIGAR
    {
        if (cigarString.empty())
        {
            cigarString = std::to_string(fragmentSize) + "M";
        }

        ALIGNPATH::path_t inputPath;
        cigar_to_apath(cigarString.c_str(), inputPath);
        edit_bam_cigar(inputPath, bamData);
    }

    // set read and qual
    {
        if ( querySeq.empty() )
        {
            querySeq = std::string(fragmentSize,'A');
        }
        const unsigned querySize(querySeq.length());
        // initialize test qual array to all Q30's:
        std::unique_ptr<uint8_t[]> qual(new uint8_t[querySize]);
        for (unsigned i(0); i<querySize; ++i)
        {
            qual[i] = 30;
        }
        edit_bam_read_and_quality(querySeq.c_str(), qual.get(), bamData);
    }

    // Set some defaults for the read
    bamRead.toggle_is_paired();
    bamRead.toggle_is_mate_fwd_strand();
    bamData.core.pos = pos;
    bamData.core.isize = fragmentSize;
    bamData.core.qual = mapQ;
    bamRead.set_target_id(targetID);

    // Set mate info
    bamData.core.mtid = mateTargetID;
    bamData.core.mpos = matePos;

    static const char nhTag[] = {'N','H'};
    static const char nmTag[] = {'N','M'};
    static const char rgTag[] = {'R','G'};
    bam_aux_append_unsigned(bamData, nhTag, 1);
    bam_aux_append_unsigned(bamData, nmTag, 1);
    bam_aux_append_unsigned(bamData, rgTag, 1);
}



void
addSupplementaryAlignmentEvidence(
        bam_record& bamRead,
        const std::string& svStr)
{
    static const char svtag[] = {'S','A'};
    bam_aux_append(bamRead.get_data(),svtag,'Z',(svStr.size()+1),
                   (uint8_t*)(svStr.c_str()));
}
