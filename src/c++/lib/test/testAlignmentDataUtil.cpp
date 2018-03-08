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
    // Add default chrom options to handle most testing conditions.
    bamHeader.chrom_to_index.insert(std::pair<std::string,int32_t>("chrM",0));
    bamHeader.chrom_to_index.insert(std::pair<std::string,int32_t>("chrT",1));
    bamHeader.chrom_data.emplace_back("chrM",1000000);
    bamHeader.chrom_data.emplace_back("chrT",1000000);

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
    if ( querySeq == "" )
    {
        for ( int i = 0; i < fragmentSize; ++i)
        {
            querySeq += "A";
        }
    }
    if ( cigarString == "" )
    {
        cigarString = std::to_string(fragmentSize) + "M";
    }

    const unsigned querySize(querySeq.length());

    // Get bam_record ptr
    ALIGNPATH::path_t inputPath;
    inputPath.push_back(ALIGNPATH::path_segment(ALIGNPATH::MATCH,querySize));
    cigar_to_apath(cigarString.c_str(),inputPath);

    bam1_t* bamDataPtr(bamRead.get_data());
    edit_bam_cigar(inputPath,*bamDataPtr);

    // initialize test qual array to all Q30's:
    std::unique_ptr<uint8_t[]> qual(new uint8_t[querySize]);
    for (unsigned i(0); i<querySize; ++i)
    {
        qual[i] = 30;
    }

    edit_bam_read_and_quality(querySeq.c_str(), qual.get(), *bamDataPtr);

    // Set some defaults for the read
    bamRead.toggle_is_paired();
    bamRead.set_target_id(0);
    bamRead.toggle_is_mate_fwd_strand();
    bamDataPtr->core.pos = pos;
    bamDataPtr->core.isize = fragmentSize;
    //bamDataPtr->core.l_qseq = fragmentSize;
    bamDataPtr->core.qual = mapQ;
    bamRead.set_target_id(targetID);

    // Set mate info
    bamDataPtr->core.mtid = mateTargetID;
    bamDataPtr->core.mpos = matePos;

    const char nhTag[] = {'N','H'};
    const char nmTag[] = {'N','M'};
    const char rgTag[] = {'R','G'};

    //bam1_t& bamRef(bamRead.get_data());
    bam_aux_append_unsigned(*bamDataPtr, nhTag, 1);
    bam_aux_append_unsigned(*bamDataPtr, nmTag, 1);
    bam_aux_append_unsigned(*bamDataPtr, rgTag, 1);
}



void
addSupplementaryAlignmentEvidence(
        bam_record& bamRead,
        const std::string& svStr)
{
    const char svtag[] = {'S','A'};
    bam_aux_append(bamRead.get_data(),svtag,'Z',(svStr.size()+1),
                   (uint8_t*)(svStr.c_str()));
}
