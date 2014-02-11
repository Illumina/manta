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

#include "manta/SVScoreInfo.hh"
#include "manta/SVScoreInfoDiploid.hh"
#include "manta/SVScoreInfoSomatic.hh"


/// all scoring info for one sv candidate, including data related to specific scoring models
///
///
struct SVModelScoreInfo
{
    void
    clear()
    {
        base.clear();
        diploid.clear();
        somatic.clear();
    }

    SVScoreInfo base;
    SVScoreInfoDiploid diploid;
    SVScoreInfoSomatic somatic;
};
