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
/// \author Chris Saunders
///

#include "EstimateSVLociRunner.hpp"
#include "SVLocusSetFinder.hpp"

#include "blt_util/input_stream_handler.hpp"
#include "blt_util/log.hpp"
#include "blt_util/time_util.hpp"
#include "htsapi/bam_header_util.hpp"
#include "manta/BamStreamerUtils.hpp"
#include "manta/SVReferenceUtil.hpp"
#include "svgraph/GenomeIntervalUtil.hpp"

#include <iostream>
#include <vector>

//#define DEBUG_ESL

EstimateSVLociRunner::EstimateSVLociRunner(const ESLOptions& opt) : _opt(opt)
{
  openBamStreams(_opt.referenceFilename, _opt.alignFileOpt.alignmentFilenames, _bamStreams);
  assertCompatibleBamStreams(opt.alignFileOpt.alignmentFilenames, _bamStreams);

  // assume bam headers are compatible after running assertCompatibleBamStreams
  const bam_hdr_t&      htslibBamHeaderInfoPtr(_bamStreams[0]->get_header());
  const bam_header_info bamHeaderInfo(htslibBamHeaderInfoPtr);

  _mergedSetPtr =
      std::make_shared<SVLocusSet>(_opt.graphOpt, bamHeaderInfo, _opt.alignFileOpt.alignmentFilenames);
}

void EstimateSVLociRunner::estimateSVLociForSingleRegion(const std::string& region)
{
  TimeTracker regionSVLocusSetBuildTimer;
  regionSVLocusSetBuildTimer.resume();

  resetBamStreamsRegion(region, _bamStreams);

  const GenomeInterval scanRegion(
      convertSamtoolsRegionToGenomeInterval(_mergedSetPtr->getBamHeader(), region));

#ifdef DEBUG_ESL
  static const std::string log_tag("EstimateSVLoci");
  log_os << log_tag << " scanRegion= " << scanRegion << "\n";
  log_os << log_tag << " startLociCount: " << _mergedSetPtr->size() << "\n";
#endif

  // grab the reference for segment we're estimating plus a buffer around the segment edges:
  static const unsigned refEdgeBufferSize(500);

  auto refSegmentPtr(std::make_shared<reference_contig_segment>());
  getIntervalReferenceSegment(
      _opt.referenceFilename,
      _mergedSetPtr->getBamHeader(),
      refEdgeBufferSize,
      scanRegion,
      (*refSegmentPtr.get()));

  SVLocusSetFinder locusFinder(_opt, scanRegion, refSegmentPtr, _mergedSetPtr);

  // loop through alignments from all samples:
  input_stream_handler sinput(mergeBamStreams(_bamStreams));
  while (sinput.next()) {
    const input_record_info current(sinput.get_current());

    if (current.itype != INPUT_TYPE::READ) {
      log_os << "ERROR: invalid input condition.\n";
      exit(EXIT_FAILURE);
    }

    const bam_streamer& readStream(*_bamStreams[current.sample_no]);
    const bam_record&   read(*(readStream.get_record_ptr()));

    locusFinder.update(readStream, read, current.sample_no);
  }

  // finished updating:
  locusFinder.flush();
  regionSVLocusSetBuildTimer.stop();
  const CpuTimes regionSVLocusSetBuildTimes(regionSVLocusSetBuildTimer.getTimes());
  _mergedSetPtr->addBuildTime(regionSVLocusSetBuildTimes);

#ifdef DEBUG_ESL
  log_os << log_tag << " endLociCount: " << _mergedSetPtr->size() << "\n";
  log_os << log_tag << " regionSVLocusSetBuildTimes: ";
  regionSVLocusSetBuildTimes.reportHr(log_os);
  log_os << "\n";
#endif
}
