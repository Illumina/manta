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
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#pragma once

#include "manta/SVModelScoreInfo.hh"
#include "format/VcfWriterSV.hh"
#include "options/CallOptionsDiploid.hh"


struct VcfWriterDiploidSV : public VcfWriterSV
{
    VcfWriterDiploidSV(
        const CallOptionsDiploid& diploidOpt,
        const bool isMaxDepthFilter,
        const std::string& referenceFilename,
        const SVLocusSet& set,
        std::ostream& os) :
        VcfWriterSV(referenceFilename,set,os),
        _diploidOpt(diploidOpt),
        _isMaxDepthFilter(isMaxDepthFilter),
        _modelScorePtr(NULL)
    {}

    void
    writeSV(
        const EdgeInfo& edge,
        const SVCandidateSetData& svData,
        const SVCandidateAssemblyData& adata,
        const SVCandidate& sv,
        const SVModelScoreInfo& ssInfo);

private:

    void
    addHeaderFormatSampleKey() const;

    void
    addHeaderInfo() const;

    void
    addHeaderFormat() const;

    void
    addHeaderFilters() const;

    void
    addSplitReadInfo(
        InfoTag_t& infotags) const;

    void
    modifyInfo(
        InfoTag_t& infotags) const;

    void
    modifySample(
        SampleTag_t& sampletags) const;

    void
    modifyTranslocInfo(
        const bool isFirstOfPair,
        InfoTag_t& infotags) const;

    void
    writeQual() const;

    void
    writeFilter() const;


    const CallOptionsDiploid& _diploidOpt;
    const bool _isMaxDepthFilter;
    const SVModelScoreInfo* _modelScorePtr;
};
