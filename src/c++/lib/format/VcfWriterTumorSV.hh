
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
/// \author Xiaoyu Chen
///

#pragma once

#include "format/VcfWriterSV.hh"
#include "format/VcfWriterScoredSV.hh"
#include "options/CallOptionsTumor.hh"


struct VcfWriterTumorSV : public VcfWriterSV, VcfWriterScoredSV
{
	static const bool isRNA = false;

    VcfWriterTumorSV(
        const CallOptionsTumor& tumorOpt,
        const bool isMaxDepthFilter,
        const std::string& referenceFilename,
        const SVLocusSet& set,
        std::ostream& os) :
        VcfWriterSV(referenceFilename, isRNA, set, os),
	    _tumorOpt(tumorOpt),
	    _isMaxDepthFilter(isMaxDepthFilter),
	    _tumorInfoPtr(nullptr)
    {}

    void
    writeSV(
        const SVCandidateSetData& svData,
        const SVCandidateAssemblyData& adata,
        const SVCandidate& sv,
        const SVId& svId,
        const SVScoreInfo& baseInfo,
        const SVScoreInfoTumor& tumorInfo,
        const EventInfo& event);

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
	writeFilter() const override;

    void
    modifySample(
        const SVCandidate& sv,
        SampleTag_t& sampletags) const override;

    void
    modifyTranslocInfo(
        const SVCandidate& sv,
        const bool isFirstOfPair,
        InfoTag_t& infotags) const override;



    const SVScoreInfoTumor&
    getTumorInfo() const
    {
        assert(NULL != _tumorInfoPtr);
        return *_tumorInfoPtr;
    }

    const CallOptionsTumor& _tumorOpt;
    const bool _isMaxDepthFilter;
    const SVScoreInfoTumor* _tumorInfoPtr;
};

