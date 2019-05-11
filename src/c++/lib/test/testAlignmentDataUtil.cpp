//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2019 Illumina, Inc.
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

#include "testAlignmentDataUtil.hpp"

#include "testConfig.h"

#include "blt_util/align_path.hpp"
#include "common/Exceptions.hpp"
#include "htsapi/align_path_bam_util.hpp"
#include "htsapi/bam_dumper.hpp"

#include <fstream>
#include <sstream>

const std::string& getTestReferenceFilename()
{
  static const std::string refPath(std::string(SHARED_TEST_DATA_PATH) + "/testGenome.fa");
  return refPath;
}

bam_header_info buildTestBamHeader()
{
  bam_header_info bamHeader;
  bamHeader.chrom_data.emplace_back("chrFoo", 500);
  bamHeader.chrom_data.emplace_back("chrBar", 500);

  int32_t chromIndex(0);
  for (const auto& chromData : bamHeader.chrom_data) {
    bamHeader.chrom_to_index.insert(std::make_pair(chromData.label, chromIndex));
    chromIndex++;
  }
  return bamHeader;
}

HtslibBamHeaderManager::HtslibBamHeaderManager(const std::vector<bam_header_info::chrom_info>& chromData)
  : _header(bam_hdr_init())
{
  _header->n_targets   = chromData.size();
  _header->target_len  = (uint32_t*)calloc(_header->n_targets, sizeof(uint32_t));
  _header->target_name = (char**)calloc(_header->n_targets, sizeof(char*));
  for (int i = 0; i < _header->n_targets; ++i) {
    _header->target_len[i]  = chromData[i].length;
    _header->target_name[i] = strdup(chromData[i].label.c_str());
  }
}

HtslibBamHeaderManager::~HtslibBamHeaderManager()
{
  bam_hdr_destroy(_header);
}

void buildTestBamFile(
    const bam_header_info&         bamHeader,
    const std::vector<bam_record>& readsToAdd,
    const std::string&             bamFilename)
{
  const HtslibBamHeaderManager bamHeaderManager(bamHeader.chrom_data);
  bam_dumper                   bamDumper(bamFilename.c_str(), bamHeaderManager.get());
  for (const bam_record& bamRecord : readsToAdd) {
    bamDumper.put_record(bamRecord.get_data());
  }
  bamDumper.close();

  const int indexStatus = bam_index_build(bamFilename.c_str(), 0);
  if (indexStatus < 0) {
    std::ostringstream oss;
    oss << "Failed to build index for bam file. bam_index_build return code: " << indexStatus
        << " bam filename: '" << bamFilename << "'\n";
    BOOST_THROW_EXCEPTION(illumina::common::GeneralException(oss.str().c_str()));
  }
}

void buildTestBamRecord(
    bam_record& bamRead,
    int         targetID,
    int         pos,
    int         mateTargetID,
    int         matePos,
    int         readLength,
    int         mapQ,
    std::string cigarString,
    std::string querySeq,
    int         fragmentSize)
{
  bam1_t& bamData(*(bamRead.get_data()));

  // set qname
  {
    edit_bam_qname("buildTestBamRecord", bamData);
  }

  // set CIGAR
  {
    if (cigarString.empty()) {
      cigarString = std::to_string(readLength) + "M";
    }

    ALIGNPATH::path_t inputPath;
    cigar_to_apath(cigarString.c_str(), inputPath);
    edit_bam_cigar(inputPath, bamData);
  }

  // set read and qual
  {
    if (querySeq.empty()) {
      querySeq = std::string(readLength, 'A');
    }
    const unsigned querySize(querySeq.length());
    // initialize test qual array to all Q30's:
    std::unique_ptr<uint8_t[]> qual(new uint8_t[querySize]);
    for (unsigned i(0); i < querySize; ++i) {
      qual[i] = 30;
    }
    edit_bam_read_and_quality(querySeq.c_str(), qual.get(), bamData);
  }
  // Set some defaults for the read
  bamRead.toggle_is_paired();
  bamRead.toggle_is_mate_fwd_strand();
  bamData.core.pos   = pos;
  bamData.core.isize = fragmentSize;
  bamData.core.qual  = mapQ;
  bamRead.set_target_id(targetID);

  // Set mate info
  bamData.core.mtid = mateTargetID;
  bamData.core.mpos = matePos;

  static const char nhTag[] = {'N', 'H'};
  static const char nmTag[] = {'N', 'M'};
  static const char rgTag[] = {'R', 'G'};
  bam_aux_append_unsigned(bamData, nhTag, 1);
  bam_aux_append_unsigned(bamData, nmTag, 1);
  bam_aux_append_unsigned(bamData, rgTag, 1);
}

void addSupplementaryAlignmentEvidence(bam_record& bamRead, const std::string& svStr)
{
  static const char svtag[] = {'S', 'A'};
  bam_aux_append(bamRead.get_data(), svtag, 'Z', (svStr.size() + 1), (const uint8_t*)(svStr.c_str()));
}

// Write two chromosomes with their depth
void buildTestChromosomeDepthFile(const std::string& depthFileName)
{
  std::ofstream outfile;
  outfile.open(depthFileName.c_str());
  outfile << "chrFoo"
          << "\t58.000" << std::endl;
  outfile << "chrBar"
          << "\t16.000" << std::endl;
  outfile.close();
}
