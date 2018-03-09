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
/// \brief Utility functions to create temporary files for unit testing
/// \author Trevor Ramsay
///

#pragma once

#include "htsapi/bam_header_info.hh"

#include <string>

/// Provide access to temporary file and limit lifetime of the temporary file to the object filetime
struct TestFileMakerBase
{
    ~TestFileMakerBase();

    const std::string&
    getFilename() const
    {
        return _tempFilename;
    }

protected:
    std::string _tempFilename;
};

/// \brief Get a non-existing temporary file name
///
/// If anything is written into the temporary file name, it will be deleted when this object goes out of scope.
struct TestFilenameMaker : public TestFileMakerBase
{
    TestFilenameMaker();
};


/// \brief Given a bam_header_info, serialize this content to a temporary file for test purposes.
///
/// The temporary file will be deleted when this object goes out of scope.
struct TestAlignHeaderFileMaker : public TestFileMakerBase
{
    explicit
    TestAlignHeaderFileMaker(const bam_header_info& info);
};

/// \brief Mock a GetAlignmentStats output file for test purposes.
///
/// The temporary file will be deleted when this object goes out of scope.
struct TestStatsFileMaker : public TestFileMakerBase
{
    TestStatsFileMaker();
};
