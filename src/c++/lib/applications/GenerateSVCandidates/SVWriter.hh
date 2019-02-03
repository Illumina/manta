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

#pragma once

#include "GSCOptions.hh"

#include "common/OutStream.hh"
#include "manta/JunctionIdGenerator.hh"
#include "manta/SVCandidateAssemblyData.hh"
#include "manta/SVCandidateSetData.hh"
#include "manta/SVMultiJunctionCandidate.hh"
#include "format/VcfWriterCandidateSV.hh"
#include "format/VcfWriterDiploidSV.hh"
#include "format/VcfWriterSomaticSV.hh"
#include "format/VcfWriterTumorSV.hh"
#include "format/VcfWriterRnaSV.hh"
#include "SVSupports.hh"


struct SVWriter
{
    SVWriter(
        const GSCOptions& initOpt,
        const bam_header_info& bamHeaderInfo,
        const char* progName,
        const char* progVersion,
        const std::vector<std::string>& sampleNames);

    void
    writeSV(
        const EdgeInfo& edge,
        const SVCandidateSetData& svData,
        const std::vector<SVCandidateAssemblyData>& mjAssemblyData,
        const SVMultiJunctionCandidate& mjSV,
        const std::vector<bool>& isCandidateJunctionFiltered,
        const std::vector<bool>& isScoredJunctionFiltered,
        const std::vector<SVId>& junctionSVId,
        const std::vector<SVModelScoreInfo>& mjModelScoreInfo,
        const SVModelScoreInfo& mjJointModelScoreInfo,
        const bool isMJEvent,
        SupportSamples& svSupports) const;

    ///////////////////////// data:
    const GSCOptions& opt;
    unsigned diploidSampleCount;

    OutStream candfs;
    OutStream dipfs;
    OutStream somfs;
    OutStream tumfs;
    OutStream rnafs;

    VcfWriterCandidateSV candWriter;
    VcfWriterDiploidSV diploidWriter;
    VcfWriterSomaticSV somWriter;
    VcfWriterTumorSV tumorWriter;
    VcfWriterRnaSV rnaWriter;
};
