// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/downloads/sequencing/licenses/>.
//

///
/// \author Chris Saunders
///

#pragma once

#include "manta/SomaticSVScoreInfo.hh"
#include "format/VcfWriterSV.hh"
#include "options/SomaticCallOptions.hh"


struct VcfWriterSomaticSV : public VcfWriterSV
{
    VcfWriterSomaticSV(
        const SomaticCallOptions& somaticOpt,
        const bool isMaxDepthFilter,
        const std::string& referenceFilename,
        const SVLocusSet& set,
        std::ostream& os) :
        VcfWriterSV(referenceFilename,set,os),
        _somaticOpt(somaticOpt),
        _isMaxDepthFilter(isMaxDepthFilter),
        _ssInfoPtr(NULL)
    {}

    void
    writeSV(
        const EdgeInfo& edge,
        const SVCandidateSetData& svData,
        const SVCandidateAssemblyData& adata,
        const unsigned svIndex,
        const SVCandidate& sv,
        const SomaticSVScoreInfo& ssInfo);

private:

    void
    addHeaderInfo() const;

    void
    addHeaderFilters() const;

    void
    modifyInfo(
        const SVBreakend& bp1,
        const SVBreakend& bp2,
        const bool isFirstOfPair,
        const SVCandidateSetData& svData,
        const SVCandidateAssemblyData& adata,
        std::vector<std::string>& infotags) const;

    std::string
    getFilter() const;

    const SomaticCallOptions& _somaticOpt;
    const bool _isMaxDepthFilter;
    const SomaticSVScoreInfo* _ssInfoPtr;
};

