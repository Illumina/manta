//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2025 Illumina, Inc.
//
// This program is licensed under the terms of the Polyform strict license
//
// ***As far as the law allows, the software comes as is, without
// any warranty or condition, and the licensor will not be liable
// to you for any damages arising out of these terms or the use
// or nature of the software, under any kind of legal claim.***
//
// You should have received a copy of the PolyForm Strict License 1.0.0
// along with this program.  If not, see <https://polyformproject.org/licenses/strict/1.0.0>.
//
//

/// \file
/// \brief Utility functions to create temporary files for unit testing
/// \author Trevor Ramsay
///

#pragma once

#include <string>

struct bam_header_info;
struct SVLocusSet;

/// Provide access to temporary file and limit lifetime of the temporary file to the object filetime
///
/// This object should not be directly instantiated by client code.
struct TestFileMakerBase {
  ~TestFileMakerBase();

  const std::string& getFilename() const { return _tempFilename; }

protected:
  // don't allow this to be directly instantiated:
  TestFileMakerBase() = default;

  std::string _tempFilename;
};

/// \brief Get a non-existing temporary file name
///
/// If the temporary file name exists in the filesystem, it will be deleted when this object goes out of
/// scope.
struct TestFilenameMaker : public TestFileMakerBase {
  TestFilenameMaker();
};

/// \brief Get a non-existing temporary file name
///
/// If the temporary file name exists in the filesystem, it will be deleted when this object goes out of
/// scope. if the temporary file name with the suffix ".bai" exists in the filesystem, this will be deleted
/// when the object goes out of scope as well.
struct BamFilenameMaker : public TestFileMakerBase {
  BamFilenameMaker();
  ~BamFilenameMaker();
};

/// \brief Given a bam_header_info, serialize this content to a temporary file for test purposes.
///
/// The temporary file will be deleted when this object goes out of scope.
struct TestAlignHeaderFileMaker : public TestFileMakerBase {
  explicit TestAlignHeaderFileMaker(const bam_header_info& info);
};

/// \brief Mock a GetAlignmentStats output file for test purposes.
///
/// Mocked stats output has 4 elements with 250 observations each for 50, 75, 100, 125
///
/// The temporary file will be deleted when this object goes out of scope.
struct TestStatsFileMaker : public TestFileMakerBase {
  TestStatsFileMaker();
};

/// \brief Given an SVLocusSet, dump its stats to a temporary file.
///
/// The temporary file will be deleted when this object goes out of scope.
struct SVLocusSetStatsFileMaker : public TestFileMakerBase {
  explicit SVLocusSetStatsFileMaker(const SVLocusSet& svLocusSet);
};

/// \brief Get a non-existing temporary file name
///
/// If the temporary file name exists in the filesystem, it will be deleted when this object goes out of
/// scope.
struct TestChromosomeDepthFileMaker : public TestFileMakerBase {
  TestChromosomeDepthFileMaker();
};
