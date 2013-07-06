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

#include "manta/GenomeInterval.hh"

#include <iosfwd>


namespace SVBreakendState
{
    enum index_t {
        UNKNOWN,    // Everything else not covered below
        RIGHT_OPEN, // 5' side of region is mapped
        LEFT_OPEN,  // 3' side of region is mapped
        COMPLEX     // A typical small scale assembly locus -- something is happening in a small region,
                    // the event might be local to that region but we don't know
    };

    inline
    const char*
    label(const index_t idx)
    {
        switch(idx)
        {
        case UNKNOWN:    return "UNKNOWN";
        case RIGHT_OPEN: return "RIGHT_OPEN";
        case LEFT_OPEN:  return "LEFT_OPEN";
        case COMPLEX:    return "COMPLEX";
        default:         return "UNKNOWN";
        }
    }
}


struct SVBreakend
{
    SVBreakend() :
        state(SVBreakendState::UNKNOWN),
        readCount(0),
        pairCount(0)
    {}

    bool
    isIntersect(const SVBreakend& rhs) const
    {
        if(state != rhs.state) return false;
        return interval.isIntersect(rhs.interval);
    }

    bool
    operator<(const SVBreakend& rhs) const
    {
        if(state < rhs.state) return true;
        if(state == rhs.state)
        {
            return (interval < rhs.interval);
        }
        return false;
    }

    bool
    merge(const SVBreakend& rhs)
    {
        if(! isIntersect(rhs)) return false;
        interval.range.merge_range(rhs.interval.range);
        readCount += rhs.readCount;
        pairCount += rhs.pairCount;
        return true;
    }

    void
    clear()
    {
        interval.clear();
        state=SVBreakendState::UNKNOWN;
        readCount=0;
        pairCount=0;
    }

    // the interval is the X% confidence interval of the SV breakend, the interface allows for
    // various probability distributions to back this interval, but these must be accessed via
    // SVCandidate:
    GenomeInterval interval;
    SVBreakendState::index_t state;

    // reads which support this SV (inclusive of read pairs)
    unsigned short readCount;

    // number of cases where both ends of a read pair pass QC and have been found to both support this SV in a consistent way:
    unsigned short pairCount;
};


std::ostream&
operator<<(std::ostream& os, const SVBreakend& svb);


struct SVCandidate
{
#if 0
    double
    breakpointProb(pos_t x, pos_t y) const;
#endif

    bool
    isIntersect(const SVCandidate& rhs) const
    {
        return ((bp1.isIntersect(rhs.bp1) && bp2.isIntersect(rhs.bp2)) ||
                (bp1.isIntersect(rhs.bp2) && bp2.isIntersect(rhs.bp1)));
    }

    bool
    merge(const SVCandidate& rhs)
    {
        if(! isIntersect(rhs)) return false;

        if(bp1.isIntersect(rhs.bp1))
        {
            bp1.merge(rhs.bp1);
            bp2.merge(rhs.bp2);
        }
        else
        {
            bp1.merge(rhs.bp2);
            bp2.merge(rhs.bp1);
        }
        return true;
    }

    void
    clear()
    {
        bp1.clear();
        bp2.clear();
    }

#if 0
    std::pair<const SVBreakend&, const SVBreakend&>
    getOrderedBreakends() const
    {
        if(bp2 < bp1) { return std::make_pair(bp2,bp1); }
        else          { return std::make_pair(bp1,bp2); }
    }
#endif

    SVBreakend bp1;
    SVBreakend bp2;
};


std::ostream&
operator<<(std::ostream& os, const SVCandidate& svc);
