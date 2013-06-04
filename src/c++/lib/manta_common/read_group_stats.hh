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
#include "common/pair_orient.hh"

#include <iosfwd>
#include <vector>
#include <string>
#include <map>

#include "boost/optional.hpp"


struct PairStatSet {
    PairStatSet() {
        clear();
    }

    void clear() {
        median = 0.;
        mean = 0.;
        sd = 0.;
    }

    double median;
    double mean;
    double sd;
};

std::ostream&
operator<<(std::ostream& os, const PairStatSet& pss);



// Read pair insert stats can be computed for each sample or read group, this
// class represents the statistics for one group:
//
struct read_group_stats {

    read_group_stats() {}
    read_group_stats(const std::string& bamFile);
    read_group_stats(const std::vector<std::string>& data);

    unsigned getReadLen(const unsigned readNum) const;

    void
    store(std::ostream& os) const;

private:
    // These data are used temporarily during ReadPairStats estimation
    struct PairStatsData {
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
    pair_orient relOrients;

    PairStatSet InsSize;

private:
    std::vector<unsigned> readLens;
};

