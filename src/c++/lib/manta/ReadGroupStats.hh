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
/// \author Bret Barnes, Xiaoyu Chen
///

#pragma once

#include "blt_util/SizeDistribution.hh"
#include "common/ReadPairOrient.hh"


/// Read pair insert stats can be computed for each sample or read group, this
/// class represents the statistics for one group:
///
struct ReadGroupStats
{
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(
        Archive& ar,
        const unsigned /*version*/)
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
