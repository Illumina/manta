// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2015 Illumina, Inc.
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
#include "format/VcfWriterScoredSV.hh"
#include "options/CallOptionsSomatic.hh"


struct VcfWriterSomaticSV : public VcfWriterSV, VcfWriterScoredSV
{
    static const bool isRNA = false;

    VcfWriterSomaticSV(
        const CallOptionsSomatic& somaticOpt,
        const bool isMaxDepthFilter,
        const std::string& referenceFilename,
        const SVLocusSet& set,
        std::ostream& os) :
        VcfWriterSV(referenceFilename, isRNA, set,os),
        _somaticOpt(somaticOpt),
        _isMaxDepthFilter(isMaxDepthFilter),
        _somaticInfoPtr(nullptr),
        _singleJunctionSomaticInfoPtr(nullptr)
    {}

    void
    writeSV(
        const SVCandidateSetData& svData,
        const SVCandidateAssemblyData& adata,
        const SVCandidate& sv,
        const SVId& svId,
        const SVScoreInfo& baseInfo,
        const SVScoreInfoSomatic& somaticInfo,
        const EventInfo& event,
        const SVScoreInfoSomatic& singleJunctionSomaticInfo);

private:

    void
    addHeaderFormatSampleKey() const override;

    void
    addHeaderInfo() const override;

    void
    addHeaderFormat() const override;

    void
    addHeaderFilters() const override;

    void
    modifyInfo(
        const EventInfo& event,
        std::vector<std::string>& infotags) const override;

    void
    modifyTranslocInfo(
        const SVCandidate& sv,
        const bool isFirstOfPair,
        std::vector<std::string>& infotags) const override;

    void
    modifySample(
        const SVCandidate& sv,
        SampleTag_t& sampletags) const override;

    void
    writeFilter() const override;

    const SVScoreInfoSomatic&
    getSomaticInfo() const
    {
        assert(NULL != _somaticInfoPtr);
        return *_somaticInfoPtr;
    }

    const SVScoreInfoSomatic&
    getSingleJunctionSomaticInfo() const
    {
        assert(NULL != _singleJunctionSomaticInfoPtr);
        return *_singleJunctionSomaticInfoPtr;
    }


    const CallOptionsSomatic& _somaticOpt;
    const bool _isMaxDepthFilter;
    const SVScoreInfoSomatic* _somaticInfoPtr;
    const SVScoreInfoSomatic* _singleJunctionSomaticInfoPtr;
};

