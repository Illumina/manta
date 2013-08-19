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
#include <ext/hash_map>

//#include "boost/unordered_map.hpp"
#include "boost/optional.hpp"
#include "boost/archive/text_oarchive.hpp"
#include "boost/archive/text_iarchive.hpp"
#include "boost/serialization/vector.hpp"
#include "boost/serialization/hash_map.hpp"


struct PairStatSet
{
    PairStatSet()
    {
        totalCount = 0;
        numOfFragSize = 0;
    }

    int totalCount;
    int numOfFragSize;
    std::vector<int> fragmentSizes;
    static const int quantileNum = 1000;
    float quantiles[quantileNum];

    //typedef boost::unordered_map<int, std::pair<int, float> > hash_map_fragment;
    typedef __gnu_cxx::hash_map <int, std::pair<int, float> > hash_map_fragment;
    hash_map_fragment fragmentSizeHash;

    bool
    calcStats();
    ///
    /// const interface used by variant callers:
    ///
    // return value for which we observe value or less with prob p
    // (not sure what the exact way to phrase this is for the discrete case)
    int
    quantile(const float p) const;

    // cdf(x)
    float
    cdf(const int x) const;

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar& totalCount;
        ar& numOfFragSize;
        ar& fragmentSizes;
        ar& quantiles;
        ar& fragmentSizeHash;
    }

};

std::ostream&
operator<<(std::ostream& os, const PairStatSet& pss);


// Read pair insert stats can be computed for each sample or read group, this
// class represents the statistics for one group:
//
struct ReadGroupStats
{

    ReadGroupStats() {}
    ReadGroupStats(const std::string& statsBamFile);

public:
    PairStatSet fragSize;
    ReadPairOrient relOrients;

    void
    write(std::ostream& os) const;

private:
    /// If PairStats has converged (or if isForcedConvergence is true)
    /// 1. All stats are computed
    /// 2. return true
    ///
    /// Otherwise:
    /// 1. only insert stats are computed
    /// 2. return false
    ///

    bool computePairStats(std::string& statsBamFile,
                          const bool isForcedConvergence = false);

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar& fragSize;
        ar& relOrients;
    }
};

