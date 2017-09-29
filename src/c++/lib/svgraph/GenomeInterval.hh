//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2017 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

#pragma once

#include "blt_util/known_pos_range2.hh"

#include <iosfwd>

///
/// \author Chris Saunders
///
/// \brief GenomeInterval identifies a unique contiguous chromosomal region.
///
/// \details GenomeInterval identifies single contiguous chromosome range. All internal locations use a chromosome
/// index number. GenomeInterval uses boost::serialize to save/load the class.
struct GenomeInterval
{
    /// \brief GenomeInverval represents a single chromosome range
    ///
    /// All internal locations use a chromosome index number
    GenomeInterval(
        const int32_t initTid = 0,
        const pos_t beginPos = 0,
        const pos_t endPos = 0) :
        tid(initTid),
        range(beginPos,endPos)
    {}

    /// \brief Identify if the GenomeIntersect Intersects with another GenomeInterval
    ///
    /// 1. The ids must be the same.
    /// 2. The range of the GenomeIntervals must overlap by more than 1 BP.
    bool
    isIntersect(const GenomeInterval& gi) const
    {
        if (tid != gi.tid) return false;
        return range.is_range_intersect(gi.range);
    }

    /// \brief Identify if the GenomeInterval range is less than the input GenomeInterval
    ///
    /// 1. start < start
    /// 2. if start == start, then end < end
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

    /// \brief Identify if the GenomeIntervals are identical
    ///
    /// 1. id == id && range == range
    bool
    operator==(const GenomeInterval& rhs) const
    {
        return ((tid==rhs.tid) && (range==rhs.range));
    }

    /// \brief Set the id to 0 and clear the chromosome range.
    void
    clear()
    {
        tid = 0;
        range.clear();
    }

    /// \brief Serialize the GenomeInterval
    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& tid& range;
    }

    /// \brief Chromsome range identifier
    int32_t tid;
    /// \brief Chromosome interval range
    known_pos_range2 range;
};

std::ostream&
operator<<(std::ostream& os, const GenomeInterval& gi);

BOOST_CLASS_IMPLEMENTATION(GenomeInterval, boost::serialization::object_serializable)
