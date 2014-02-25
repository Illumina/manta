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
/// \brief align a contig across two breakend regions
///


#pragma once

#include "Alignment.hh"
#include "blt_util/align_path.hh"


struct AlignerUtil
{
    static
    void
    updatePath(
        ALIGNPATH::path_t& path,
        ALIGNPATH::path_segment& ps,
        ALIGNPATH::align_t atype)
    {
        if (ps.type == atype) return;
        if (ps.type != ALIGNPATH::NONE) path.push_back(ps);
        ps.type = atype;
        ps.length = 0;
    }
};



/// bookkeeping variables used during alignment backtrace
template <typename ScoreType>
struct BackTrace
{
    BackTrace() :
        max(0),
        state(AlignState::MATCH),
        queryBegin(0),
        refBegin(0),
        isInit(false)
    {}

    ScoreType max;
    AlignState::index_t state;
    unsigned queryBegin,refBegin;
    bool isInit;
};



template <typename ScoreType>
void
updateBacktrace(
    const ScoreType thisMax,
    const unsigned refIndex,
    const unsigned queryIndex,
    BackTrace<ScoreType>& btrace)
{
    if ( (! btrace.isInit) || (thisMax>btrace.max))
    {
        btrace.max=thisMax;
        btrace.refBegin=refIndex;
        btrace.queryBegin=queryIndex;
        btrace.isInit=true;
    }
}





