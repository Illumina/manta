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


#pragma once

#include "blt_util/id_map.hh"
#include "common/ReadPairOrient.hh"

#include <iosfwd>
#include <vector>
#include <string>
#include <map>

#include "boost/optional.hpp"


struct PairStatSet
{
    PairStatSet()
    {
        clear();
    }

    ///
    /// const interface used by variant callers:
    ///

    // return value for which we observe value or less with prob p
    // (not sure what the exact way to phrase this is for the discrete case)
    double
    quantile(const double p) const;

    // cdf(x)
    double
    cdf(const double x) const;


    ///
    /// remainder of interface for estimation, store/read from disk
    ///
    void clear()
    {
        median = 0.;
        sd = 0.;
    }

    /// TODO: hide implementation details, better estimate:
    double median;
    double sd;
};

std::ostream&
operator<<(std::ostream& os, const PairStatSet& pss);



// Read pair insert stats can be computed for each sample or read group, this
// class represents the statistics for one group:
//
struct ReadGroupStats
{

    ReadGroupStats() {}
    ReadGroupStats(const std::string& bamFile);
    ReadGroupStats(const std::vector<std::string>& data);

    void
    write(std::ostream& os) const;

private:
    // These data are used temporarily during ReadPairStats estimation
    struct PairStatsData
    {
        std::vector<int32_t> fragmentLengths;
    };

    /// If PairStats has converged (or if isForcedConvergence is true)
    /// 1. All stats are computed
    /// 2. return true
    ///
    /// Otherwise:
    /// 1. only insert stats are computed
    /// 2. return false
    ///
    bool computePairStats(PairStatsData& psd, const bool isForcedConvergence = false);

public:
    //////////////// data:
    ReadPairOrient relOrients;

    PairStatSet fragSize;
};

