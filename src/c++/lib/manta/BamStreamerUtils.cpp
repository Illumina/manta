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

#include "BamStreamerUtils.hpp"

#include "common/Exceptions.hpp"
#include "htsapi/bam_header_util.hpp"

#include <sstream>

void openBamStreams(
    const std::string&                          referenceFilename,
    const std::vector<std::string>&             bamFilenames,
    std::vector<std::shared_ptr<bam_streamer>>& bamStreams)
{
  bamStreams.clear();
  for (const std::string& bamFilename : bamFilenames) {
    bamStreams.emplace_back(new bam_streamer(bamFilename.c_str(), referenceFilename.c_str()));
  }
}

void resetBamStreamsRegion(const std::string& region, std::vector<std::shared_ptr<bam_streamer>>& bamStreams)
{
  if (region.empty()) return;
  for (auto& bamStream : bamStreams) {
    bamStream->resetRegion(region.c_str());
  }
}

void assertCompatibleBamStreams(
    const std::vector<std::string>&                   bamFilenames,
    const std::vector<std::shared_ptr<bam_streamer>>& bamStreams)
{
  const unsigned bamCount(bamStreams.size());

  assert(0 != bamCount);

  if (bamCount < 2) return;

  const bam_hdr_t& compareHeader(bamStreams[0]->get_header());
  for (unsigned bamIndex(1); bamIndex < bamCount; ++bamIndex) {
    const bam_hdr_t& indexHeader(bamStreams[bamIndex]->get_header());
    if (!check_header_compatibility(compareHeader, indexHeader)) {
      std::ostringstream oss;
      oss << "Incompatible headers between alignment files:\n"
          << "\t" << bamFilenames[0] << "\n"
          << "\t" << bamFilenames[bamIndex] << "\n";
      BOOST_THROW_EXCEPTION(illumina::common::GeneralException(oss.str()));
    }
  }
}

input_stream_data mergeBamStreams(std::vector<std::shared_ptr<bam_streamer>>& bamStreams)
{
  int               bamIndex(0);
  input_stream_data sdata;
  for (auto& bamStream : bamStreams) {
    sdata.register_reads(*bamStream, bamIndex);
    bamIndex++;
  }
  return sdata;
}
