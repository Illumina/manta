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

///
/// \author Chris Saunders
///

#pragma once

#include "blt_util/known_pos_range2.hh"

#include <iosfwd>


/// single chromosome range
///
/// all internal locations use a chromosome index number
struct GenomeInterval
{
    GenomeInterval(
        const int32_t initTid = 0,
        const pos_t beginPos = 0,
        const pos_t endPos = 0) :
        tid(initTid),
        range(beginPos,endPos)
    {}

    /// does this intersect a second GenomeInterval?
    bool
    isIntersect(const GenomeInterval& gi) const
    {
        if (tid != gi.tid) return false;
        return range.is_range_intersect(gi.range);
    }

    bool
    operator<(const GenomeInterval& rhs) const
    {
        if (tid<rhs.tid) return true;
        if (tid == rhs.tid)
        {
            return (range<rhs.range);
        }
        return false;
    }

    bool
    operator==(const GenomeInterval& rhs) const
    {
        return ((tid==rhs.tid) && (range==rhs.range));
    }

    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& tid& range;
    }

    int32_t tid;
    known_pos_range2 range;
};


std::ostream&
operator<<(std::ostream& os, const GenomeInterval& gi);

BOOST_CLASS_IMPLEMENTATION(GenomeInterval, boost::serialization::object_serializable)
