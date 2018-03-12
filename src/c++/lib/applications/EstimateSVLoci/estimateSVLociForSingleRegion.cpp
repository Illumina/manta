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
#include "htsapi/bam_streamer.hh"
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

    typedef std::shared_ptr<bam_streamer> stream_ptr;
    std::vector<stream_ptr> bamStreams;

    // setup all data for main alignment loop:
    for (const std::string& alignmentFilename : opt.alignFileOpt.alignmentFilenames)
    {
        stream_ptr tmp(new bam_streamer(alignmentFilename.c_str(), opt.referenceFilename.c_str()));
        if (! region.empty())
        {
            tmp->resetRegion(region.c_str());
        }
        bamStreams.push_back(tmp);
    }

    const unsigned bamCount(bamStreams.size());

    assert(0 != bamCount);

    // check bam header compatibility:
    if (bamCount > 1)
    {
        /// TODO: provide a better error exception for failed bam header check:
        const bam_hdr_t& compareHeader(bamStreams[0]->get_header());
        for (unsigned bamIndex(1); bamIndex<bamCount; ++bamIndex)
        {
            const bam_hdr_t& indexHeader(bamStreams[bamIndex]->get_header());
            if (! check_header_compatibility(compareHeader,indexHeader))
            {
                log_os << "ERROR: incompatible bam headers between files:\n"
                       << "\t" << opt.alignFileOpt.alignmentFilenames[0] << "\n"
                       << "\t" << opt.alignFileOpt.alignmentFilenames[bamIndex] << "\n";
                exit(EXIT_FAILURE);
            }
        }
    }

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

    input_stream_data sdata;
    for (unsigned bamIndex(0); bamIndex<bamCount; ++bamIndex)
    {
        sdata.register_reads(*bamStreams[bamIndex],bamIndex);
    }

    // loop through alignments:
    input_stream_handler sinput(sdata);
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
