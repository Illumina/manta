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

#include "BamStreamerUtils.hh"



void
openBamStreams(
    const std::string& referenceFilename,
    const std::vector<std::string>& bamFilenames,
    std::vector<std::shared_ptr<bam_streamer>>& bamStreams)
{
    bamStreams.clear();
    for (const std::string& bamFilename : bamFilenames)
    {
        bamStreams.emplace_back(new bam_streamer(bamFilename.c_str(), referenceFilename.c_str()));
    }
}



void
resetBamStreamsRegion(
    const std::string& region,
    std::vector<std::shared_ptr<bam_streamer>>& bamStreams)
{
    if (region.empty()) return;
    for (auto& bamStream : bamStreams)
    {
        bamStream->resetRegion(region.c_str());
    }
}



input_stream_data
mergeBamStreams(
    std::vector<std::shared_ptr<bam_streamer>>& bamStreams)
{
    int bamIndex(0);
    input_stream_data sdata;
    for (auto& bamStream : bamStreams)
    {
        sdata.register_reads(*bamStream, bamIndex);
        bamIndex++;
    }
    return sdata;
}