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

#include "testFileMakers.hh"

#include "manta/ReadGroupStatsSet.hh"
#include "test/testUtil.hh"

#include "boost/filesystem.hpp"

#include <fstream>



TestFileMakerBase::
~TestFileMakerBase()
{
    using namespace boost::filesystem;
    if (exists(_tempFilename))
    {
        remove(_tempFilename);
    }
}


TestFilenameMaker::
TestFilenameMaker()
{
    _tempFilename = getNewTempFile();
}



/// \brief Create an alignment file based on the bam_header_info input
/// \return Name of temporary file containing alignment fileoutput
static
std::string
createTestAlignHeaderFile(
    const bam_header_info& info)
{
    const std::string tempFile(getNewTempFile());
    {
        std::ofstream outfile(tempFile);
        assert(outfile);
        outfile << info;
    }
    return tempFile;
}



TestAlignHeaderFileMaker::
TestAlignHeaderFileMaker(const bam_header_info& info)
{
    _tempFilename = createTestAlignHeaderFile(info);
}



/// \brief Mock a GetAlignmentStats output file for testing.
///
/// Default has 4 elements with 250 observations each for 50, 75, 100, 125
///
/// \return Name of temporary file containing mock stats output
static
std::string
createTestStatsFile()
{
    const std::string tempFile(getNewTempFile());

    ReadGroupLabel rgKey("tempStatsGroup", "");
    ReadGroupStats rgStats;
    for (unsigned i(0);i<250; ++i)
    {
        rgStats.fragStats.addObservation(50);
        rgStats.fragStats.addObservation(75);
        rgStats.fragStats.addObservation(100);
        rgStats.fragStats.addObservation(125);
    }

    ReadGroupStatsSet rstats;
    rstats.setStats(rgKey, rgStats);

    rstats.save(tempFile.c_str());

    return tempFile;
}



TestStatsFileMaker::
TestStatsFileMaker()
{
    _tempFilename = createTestStatsFile();
}
