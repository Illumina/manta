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

#pragma once

#include "blt_types.hh"
#include "known_pos_range2.hh"

#include "boost/icl/interval_map.hpp"


/// accumulates region specific sum(T) give a set of (region,T) associations
///
template <typename T>
struct RegionSum
{
    void
    clear()
    {
        _map.clear();
    }

    void
    add(
        const known_pos_range2 pr,
        const T val)
    {
        _map.add(std::make_pair(boost::icl::interval<pos_t>::right_open(pr.begin_pos(),pr.end_pos()),val));
    }

    /// return peak value for all regions:
    T
    maxVal() const
    {
        T max(0);
        bool isFirst(true);
        for (const auto& val : _map)
        {
            if (isFirst || val.second > max)
            {
                max = val.second;
                isFirst = false;
            }
        }
        return max;
    }

private:
    using map_t = boost::icl::interval_map<pos_t,T>;
    map_t _map;
};
