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
/// \author Chris Saunders
///

#include "estimateSVLociForSingleRegion.hh"
#include "SVLocusSetFinder.hh"

#include "blt_util/input_stream_handler.hh"
#include "blt_util/log.hh"
#include "blt_util/time_util.hh"
#include "htsapi/bam_header_util.hh"
#include "manta/BamStreamerUtils.hh"
#include "manta/SVReferenceUtil.hh"
#include "svgraph/GenomeIntervalUtil.hh"

#include <iostream>
#include <vector>


//#define DEBUG_ESL



void
estimateSVLociForSingleRegion(
    const ESLOptions& opt,
    const std::string& region,
    SVLocusSet& mergedSet)
{
    TimeTracker regionSVLocusSetBuildTimer;
    regionSVLocusSetBuildTimer.resume();

    std::vector<std::shared_ptr<bam_streamer>> bamStreams;
    openBamStreams(opt.referenceFilename, opt.alignFileOpt.alignmentFilenames, bamStreams);
    resetBamStreamsRegion(region, bamStreams);
    assertCompatibleBamStreams(opt.alignFileOpt.alignmentFilenames, bamStreams);

    // assume headers compatible after this point....

    const bam_hdr_t& header(bamStreams[0]->get_header());
    const auto bamHeaderPtr(std::make_shared<bam_header_info>(header));

    // Use the const reference for older objects which have an implicit lifetime contract with the caller
    const bam_header_info& bamHeader(*(bamHeaderPtr.get()));
    const GenomeInterval scanRegion(convertSamtoolsRegionToGenomeInterval(bamHeader, region));

#ifdef DEBUG_ESL
    static const std::string log_tag("EstimateSVLoci");
    log_os << log_tag << " scanRegion= " << scanRegion << "\n";
    log_os << log_tag << " startLociCount: " << mergedSet.size() << "\n";
#endif

    // grab the reference for segment we're estimating plus a buffer around the segment edges:
    static const unsigned refEdgeBufferSize(500);

    auto refSegmentPtr(std::make_shared<reference_contig_segment>());
    getIntervalReferenceSegment(opt.referenceFilename, bamHeader, refEdgeBufferSize, scanRegion, (*refSegmentPtr.get()));

    SVLocusSetFinder locusFinder(opt, scanRegion, bamHeaderPtr, refSegmentPtr, mergedSet);

    // loop through alignments from all samples:
    input_stream_handler sinput(mergeBamStreams(bamStreams));
    while (sinput.next())
    {
        const input_record_info current(sinput.get_current());

        if (current.itype != INPUT_TYPE::READ)
        {
            log_os << "ERROR: invalid input condition.\n";
            exit(EXIT_FAILURE);
        }

        const bam_streamer& readStream(*bamStreams[current.sample_no]);
        const bam_record& read(*(readStream.get_record_ptr()));

        locusFinder.update(readStream, read, current.sample_no);
    }

    // finished updating:
    locusFinder.flush();
    regionSVLocusSetBuildTimer.stop();
    const CpuTimes regionSVLocusSetBuildTimes(regionSVLocusSetBuildTimer.getTimes());
    mergedSet.addBuildTime(regionSVLocusSetBuildTimes);

#ifdef DEBUG_ESL
    log_os << log_tag << " endLociCount: " << mergedSet.size() << "\n";
    log_os << log_tag << " regionSVLocusSetBuildTimes: ";
    regionSVLocusSetBuildTimes.reportHr(log_os);
    log_os << "\n";
#endif
}
