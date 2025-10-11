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
