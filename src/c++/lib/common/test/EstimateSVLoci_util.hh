//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2017 Illumina, Inc.
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


///
/// \file
/// \author Trevor Ramsay
///

#pragma once

#include "htsapi/align_path_bam_util.hh"
#include "htsapi/SimpleAlignment_bam_util.hh"
#include "htsapi/bam_record_util.hh"
#include "htsapi/bam_record.hh"
#include "manta/SVLocusScanner.hh"

#include "boost/scoped_array.hpp"
#include "boost/test/unit_test.hpp"
#include "boost/algorithm/string.hpp"

#include <fstream>
#include <iostream>

/// Wrap string for handling input streams
class line : public std::string {};

std::istream &operator>>(std::istream &is, line &l)
{
    std::getline(is, l);
    return is;
}

/// \brief Create an alignment file based on the bam_header_info input
void
createAlignFile(
        const bam_header_info& info,
        const std::string& filename)
{
    if (FILE* file = fopen(filename.c_str(), "r")) {
        fclose(file);
        std::remove(filename.c_str());
    }
    std::ofstream outfile(filename);
    outfile << info;
    outfile.close();
}

/// \brief Create a stats file for testing.
///        Default has 4 elements with 250 observations each for 50, 75, 100, 125
void
createStatsFile(
        const std::string& filename,
        std::string fileInput = "")
{
    if ( fileInput.length() < 1) {
        fileInput = "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\" ?><!DOCTYPE boost_serialization><boost_serialization signature=\"serialization::archive\" version=\"11\">"
                "<numGroups>1</numGroups>\n"
                "<groupStats_0>\n"
                "<groupLabel>tempStatsGroup</groupLabel>\n"
                "<groupStats>\n"
                "<fragmentSizeDistribution>\n"
                "<totalObservationCount>1000</totalObservationCount>\n"
                "<elementCount>4</elementCount>\n"
                "<element>\n"
                "<size>50</size>\n"
                "<count>250</count>\n"
                "</element>"
                "<element>\n"
                "<size>75</size>\n"
                "<count>250</count>\n"
                "</element>"
                "<element>\n"
                "<size>100</size>\n"
                "<count>250</count>\n"
                "</element>"
                "<element>\n"
                "<size>125</size>\n"
                "<count>250</count>\n"
                "</element>"
                "</fragmentSizeDistribution>\n"
                "<pairOrientation>\n"
                "<pairOrientation>Rp</pairOrientation>\n"
                "</pairOrientation>\n"
                "<readCount>\n"
                "<totalReadCount>1000000</totalReadCount>\n"
                "<totalPairedReadCount>1000000</totalPairedReadCount>\n"
                "<totalUnpairedReadCount>0</totalUnpairedReadCount>\n"
                "<totalPairedLowMapqReadCount>0</totalPairedLowMapqReadCount>\n"
                "<totalHighConfidenceReadPairCount>1000000</totalHighConfidenceReadPairCount>\n"
                "</readCount>\n"
                "</groupStats>\n"
                "</groupStats_0>\n"
                "</boost_serialization>";
    }

    if (FILE* file = fopen(filename.c_str(), "r")) {
        fclose(file);
        std::remove(filename.c_str());
    }
    std::ofstream outfile(filename);
    if ( fileInput.length() > 0 ) {
        outfile << fileInput;
    }
    else
    {
        outfile << "This is a temporary file for unit testing" << std::endl;
    }
    outfile.close();
}

/// \brief Read the given stats file looking for a specific row of information.
std::string
getResultFromStatsFile(
        std::string statsFile,
        std::string category)
{
    std::ifstream is(statsFile);
    std::istream_iterator<line> start(is), end;
    std::vector<line> stats(start, end);

    // iterate through the file
    for ( std::vector<line>::iterator it=stats.begin(); it!=stats.end(); ++it) {
        // find all rows that contain the substring
        if ((*it).find(category) != std::string::npos ) {
            std::vector<std::string> strs;
            boost::split(strs, (*it), boost::is_any_of("\t"));
            // If the substring is exactly correct (avoid issues like NoFilter vs NoFilterAndIgnored)
            if ( strs[0].length() ==  category.length()) {
                return strs[1];
            }
        }
    }
    return "";
}

/// \brief Change the templateSize of the bam record.
void
changeTemplateSize(
        bam_record& bamRead,
        int newSize)
{
    bam1_t* bamDataPtr(bamRead.get_data());
    bamDataPtr->core.isize = newSize;
}

/// \brief Add supplementary evidence to a given read, defaulting to chrT from the default bam_header_info.
void addSupplementaryEvidence(
        bam_record &bamRead,
        std::string svStr = "chrT,6348,-,54H22M,50,0;")
{
    const char svtag[] = {'S','A'};
    bam_aux_append(bamRead.get_data(),svtag,'Z',(svStr.size()+1),
                   (uint8_t*)(svStr.c_str()));
}

/// \brief Create the BamAlignment for a given read using the SimpleAlignment object.
void
getBamAlignment(
        const bam_record& bamRead,
        SimpleAlignment& al)
{
    al.is_fwd_strand=bamRead.is_fwd_strand();
    al.tid=bamRead.target_id();
    al.pos=(bamRead.pos()-1);

    bam_cigar_to_apath(bamRead.raw_cigar(),bamRead.n_cigar(),al.path);
}

/// \brief Generic function to create a file with the given input
void
createFile(
        const std::string& filename,
        const std::string& fileInput = "")
{
    if (FILE* file = fopen(filename.c_str(), "r"))
    {
        fclose(file);
        std::remove(filename.c_str());
    }
    std::ofstream outfile(filename);
    if ( fileInput.length() > 0 )
    {
        outfile << fileInput;
    }
    else
    {
        outfile << "This is a temporary file for unit testing" << std::endl;
    }
    outfile.close();
}

/// \brief Return a bam_header_info object with references 0 and 1: "chrM" and "chrT", both size 1000000
bam_header_info
buildBamHeader()
{
    bam_header_info bamHeader = bam_header_info();
    // Add default chrom options to handle most testing conditions.
    bamHeader.chrom_to_index.insert(std::pair<std::string,int32_t>("chrM",0));
    bamHeader.chrom_to_index.insert(std::pair<std::string,int32_t>("chrT",1));
    bamHeader.chrom_data.emplace_back("chrM",1000000);
    bamHeader.chrom_data.emplace_back("chrT",1000000);

    return bamHeader;
}

/// \brief Build a bam_record based on the input parameters.
///        A default bam_record is a proper paired reference sequence.
void
buildBamRecord(
        bam_record& bamRead,
        int targetID = 0,
        int pos = 100,
        int mateTargetID = 0,
        int matePos = 200,
        int fragmentSize = 100,
        int mapQ = 15,
        std::string cigarString = "",
        std::string querySeq = "")
{
    if ( querySeq == "")
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
    boost::scoped_array<uint8_t> qual(new uint8_t[querySize]);
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
}

/// \brief SVLocusScanner builder
///        Default: Creates a DNA SVLocusScanner with minCandidateVariantSize = 8 and simple
///        stats and align files (Refer to createStatsFile() and createAlignFile respectively).
std::unique_ptr<SVLocusScanner>
buildSVLocusScanner(
        const bam_header_info& bamHeaderInfo,
        const std::string& statsFilename = std::string("tempStatsFile.txt"),
        const std::string& alignFilename = std::string("tempAlignFile.txt"),
        bool isRNA = false,
        int minCandidateVariantSizeInput = 8)
{
    ReadScannerOptions opts = ReadScannerOptions();
    opts.minCandidateVariantSize = minCandidateVariantSizeInput;
    const ReadScannerOptions& constRefOpts(opts);

    createStatsFile(statsFilename);
    createAlignFile(bamHeaderInfo, alignFilename);
    const std::vector<std::string> alignFilenameVector(1, alignFilename);
    bool& refIsRna(isRNA);

    std::unique_ptr<SVLocusScanner> newSVLS(new SVLocusScanner(constRefOpts, statsFilename, alignFilenameVector, refIsRna));

    return newSVLS;
}
