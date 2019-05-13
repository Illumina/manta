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
/// \brief Utility functions for mocking bam alignment data during unit testing
/// \author Trevor Ramsay
///

#pragma once

#include "htsapi/bam_header_info.hpp"
#include "htsapi/bam_record.hpp"

#include <string>
#include <vector>

/// \brief Return path to a file containing a simplified test reference
///
/// The test reference contains "chrFoo" and "chrBar", both with size 500
const std::string& getTestReferenceFilename();

/// \brief Return a test bam_header_info object
///
/// The reference structure matches the that returned by getTestReferenceFilename
bam_header_info buildTestBamHeader();

/// \brief Build and manage the lifetime of the htslib bam header structure, bam_hdr_t
///
/// Note that the htslib header struct created by this object has only the minimum information filled in
/// required for testing.
struct HtslibBamHeaderManager {
  /// \brief Initialize a bam header object with chromosome name, size and order as specified by \p chromData
  explicit HtslibBamHeaderManager(const std::vector<bam_header_info::chrom_info>& chromData = {});

  ~HtslibBamHeaderManager();

  const bam_hdr_t& get() const { return *_header; }

private:
  bam_hdr_t* _header;
};

/// \brief Write bamHeader and readsToAdd into new bam file
///
/// \param bamHeader Reference contig order and length
/// \param readsToAdd Bam records
/// \param bamFilename Name of file to write bam output to
void buildTestBamFile(
    const bam_header_info&         bamHeader,
    const std::vector<bam_record>& readsToAdd,
    const std::string&             bamFilename);

/// \brief Build a bam_record based on the input parameters.
///
/// A default bam_record is a proper paired reference sequence.
void buildTestBamRecord(
    bam_record& bamRead,
    int         targetID     = 0,
    int         pos          = 100,
    int         mateTargetID = 0,
    int         matePos      = 200,
    int         readLength   = 100,
    int         mapQ         = 15,
    std::string cigarString  = "",
    std::string querySeq     = "",
    int         fragmentSize = 400);

/// \brief Add supplementary alignment evidence to a bam record
///
/// This fills in an auxiliary "SA" tag, defaulting to a supplementary connections to "chrT"
/// from the default bam_header_info.
void addSupplementaryAlignmentEvidence(
    bam_record& bamRead, const std::string& svStr = "chrFoo,300,-,54H22M,50,0;");

/// \brief Change the templateSize of the bam record.
inline void changeTemplateSize(bam_record& bamRead, int newSize)
{
  bam1_t* bamDataPtr(bamRead.get_data());
  bamDataPtr->core.isize = newSize;
}

/// \brief Write chromosome depth to a file
///
/// \param depthFileName Name of file to write depth output to
void buildTestChromosomeDepthFile(const std::string& depthFileName);
