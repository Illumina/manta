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

#include "testFileMakers.hpp"

#include "htsapi/bam_header_info.hpp"
#include "manta/ReadGroupStatsSet.hpp"
#include "svgraph/SVLocusSet.hpp"
#include "test/testUtil.hpp"

#include "boost/filesystem.hpp"

#include <cassert>
#include <fstream>

TestFileMakerBase::~TestFileMakerBase()
{
  using namespace boost::filesystem;
  if (exists(_tempFilename)) {
    remove(_tempFilename);
  }
}

TestFilenameMaker::TestFilenameMaker()
{
  _tempFilename = getNewTempFile();
}

BamFilenameMaker::BamFilenameMaker()
{
  _tempFilename = getNewTempFile() + ".bam";
}

BamFilenameMaker::~BamFilenameMaker()
{
  using namespace boost::filesystem;
  const std::string& indexFilename(_tempFilename + ".bai");
  if (exists(indexFilename)) {
    remove(indexFilename);
  }
}

TestAlignHeaderFileMaker::TestAlignHeaderFileMaker(const bam_header_info& info)
{
  _tempFilename = getNewTempFile();
  std::ofstream os(_tempFilename);
  assert(os);
  os << info;
}

TestStatsFileMaker::TestStatsFileMaker()
{
  _tempFilename = getNewTempFile();

  ReadGroupLabel rgKey("tempStatsGroup", "");
  ReadGroupStats rgStats;
  for (unsigned i(0); i < 250; ++i) {
    rgStats.fragStats.addObservation(50);
    rgStats.fragStats.addObservation(75);
    rgStats.fragStats.addObservation(100);
    rgStats.fragStats.addObservation(125);
  }

  ReadGroupStatsSet rstats;
  rstats.setStats(rgKey, rgStats);

  rstats.save(_tempFilename.c_str());
}

SVLocusSetStatsFileMaker::SVLocusSetStatsFileMaker(const SVLocusSet& svLocusSet)
{
  _tempFilename = getNewTempFile();
  std::ofstream os(_tempFilename);
  assert(os);
  svLocusSet.dumpStats(os);
}

TestChromosomeDepthFileMaker::TestChromosomeDepthFileMaker()
{
  _tempFilename = getNewTempFile() + ".txt";
}
