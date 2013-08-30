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

#include "common/ReadPairOrient.hh"
#include "manta/SizeDistribution.hh"

#include <string>


/// Read pair insert stats can be computed for each sample or read group, this
/// class represents the statistics for one group:
///
struct ReadGroupStats
{

    ReadGroupStats() {}
    ReadGroupStats(const std::string& statsBamFile);

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
    void serialize(Archive& ar, const unsigned /*version*/)
    {
        ar& boost::serialization::make_nvp("fragmentSizeDistribution", fragStats);
        ar& boost::serialization::make_nvp("pairOrientation", relOrients);
    }

    ///////////////////////////// data:
public:
    SizeDistribution fragStats;
    ReadPairOrient relOrients;

};

BOOST_CLASS_IMPLEMENTATION(ReadGroupStats, boost::serialization::object_serializable)

