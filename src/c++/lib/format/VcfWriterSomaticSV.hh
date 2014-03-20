// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2014 Illumina, Inc.
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
#include "options/CallOptionsSomatic.hh"


struct VcfWriterSomaticSV : public VcfWriterSV
{
    VcfWriterSomaticSV(
        const CallOptionsSomatic& somaticOpt,
        const bool isMaxDepthFilter,
        const std::string& referenceFilename,
        const SVLocusSet& set,
        std::ostream& os) :
        VcfWriterSV(referenceFilename,set,os),
        _somaticOpt(somaticOpt),
        _isMaxDepthFilter(isMaxDepthFilter),
        _baseInfoPtr(NULL),
        _somaticInfoPtr(NULL)
    {}

    void
    writeSV(
        const SVCandidateSetData& svData,
        const SVCandidateAssemblyData& adata,
        const SVCandidate& sv,
        const SVId& svId,
        const SVScoreInfo& baseInfo,
        const SVScoreInfoSomatic& somaticInfo);

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
    modifyInfo(
        std::vector<std::string>& infotags) const;

    void
    modifyTranslocInfo(
        const bool isFirstOfPair,
        std::vector<std::string>& infotags) const;

    void
    modifySample(
        const SVCandidate& sv,
        SampleTag_t& sampletags) const;

    void
    writeFilter() const;

    const SVScoreInfo&
    getBaseInfo() const
    {
        assert(NULL != _baseInfoPtr);
        return *_baseInfoPtr;
    }

    const SVScoreInfoSomatic&
    getSomaticInfo() const
    {
        assert(NULL != _somaticInfoPtr);
        return *_somaticInfoPtr;
    }


    const CallOptionsSomatic& _somaticOpt;
    const bool _isMaxDepthFilter;
    const SVScoreInfo* _baseInfoPtr;
    const SVScoreInfoSomatic* _somaticInfoPtr;
};

