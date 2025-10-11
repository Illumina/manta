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
