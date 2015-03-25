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
/// \author Chris Saunders
///

#pragma once

#include "boost/serialization/nvp.hpp"

#include <cstdint>

#include <iosfwd>


struct SVFinderStats
{
    SVFinderStats() {}

    void
    merge(
        const SVFinderStats& rhs)
    {
        edgeFilter += rhs.edgeFilter;
        semiMappedFilter += rhs.semiMappedFilter;
        ComplexLowCountFilter += rhs.ComplexLowCountFilter;
        ComplexLowSignalFilter += rhs.ComplexLowSignalFilter;
    }

    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& BOOST_SERIALIZATION_NVP(edgeFilter)
        & BOOST_SERIALIZATION_NVP(semiMappedFilter)
        & BOOST_SERIALIZATION_NVP(ComplexLowCountFilter)
        & BOOST_SERIALIZATION_NVP(ComplexLowSignalFilter)
        ;
    }

    void
    report(std::ostream& os) const;


    uint64_t edgeFilter = 0;
    uint64_t semiMappedFilter = 0;
    uint64_t ComplexLowCountFilter = 0;
    uint64_t ComplexLowSignalFilter = 0;
};

BOOST_CLASS_IMPLEMENTATION(SVFinderStats, boost::serialization::object_serializable)

