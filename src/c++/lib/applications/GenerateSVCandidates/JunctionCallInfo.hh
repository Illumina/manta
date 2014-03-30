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
/// \author Chris Saunders and Xiaoyu Chen
///

#pragma once

#include "SVEvidence.hh"
#include "manta/SVScoreInfo.hh"
#include "manta/SVCandidate.hh"


/// manage all per-junction information consumed by an SV calling model
///
/// using this object facilities multi-breakend event scoring, but clearly
/// separating out per-junction input info from junction-independent info
///
struct JunctionCallInfo
{
    JunctionCallInfo() :
        _sv(NULL),
        _evidence(NULL),
        _baseInfo(NULL),
        _spanningPairWeight(0)
    {}

    const SVCandidate&
    getSV() const
    {
        assert(NULL != _sv);
        return *_sv;
    }

    const SVEvidence&
    getEvidence() const
    {
        assert(NULL != _evidence);
        return *_evidence;
    }

    const SVScoreInfo&
    getBaseInfo() const
    {
        assert(NULL != _baseInfo);
        return *_baseInfo;
    }

    float
    getSpanningWeight() const
    {
        return _spanningPairWeight;
    }

    void
    init(
        const SVCandidate& sv,
        const SVEvidence& evidence,
        const SVScoreInfo& baseInfo,
        const float spanningPairWeight)
    {
        _sv=&sv;
        _evidence=&evidence;
        _baseInfo=&baseInfo;
        _spanningPairWeight=spanningPairWeight;
    }

private:
    const SVCandidate* _sv;
    const SVEvidence* _evidence;
    const SVScoreInfo* _baseInfo;
    float _spanningPairWeight;
};
