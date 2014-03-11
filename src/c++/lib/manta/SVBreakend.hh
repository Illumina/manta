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

#include "svgraph/GenomeInterval.hh"

#include <cassert>

#include <iosfwd>

namespace SVEvidenceType
{
enum index_t
{
    PAIR,        /// a pair observation based on both read BAM records
    LOCAL_PAIR,  /// a pair observation based on one read, htere
    CIGAR,
    SOFTCLIP,
    SEMIALIGN,
    SHADOW,
    SPLIT_ALIGN,
    UNKNOWN,    /// temporary state
    SIZE
};

inline
const char*
label(const int i)
{
    switch (i)
    {
    case PAIR:
        return "pair";
    case LOCAL_PAIR:
        return "local_pair";
    case CIGAR:
        return "cigar";
    case SOFTCLIP:
        return "softclip";
    case SEMIALIGN:
        return "semialign";
    case SHADOW:
        return "shadow";
    case SPLIT_ALIGN:
        return "split_align";
    case UNKNOWN:
        return "unknown";
    default:
        assert(false && "Unknown SVCandidate evidence type");
        return NULL;
    }
}

inline
bool
isPairType(const int i)
{
    switch (i)
    {
    case PAIR:
    case LOCAL_PAIR:
        return true;
    default:
        return false;
    }
}
}



/// enumerate candidate evidence counts:
struct SVBreakendLowResEvidence
{
    typedef SVBreakendLowResEvidence self_t;

    SVBreakendLowResEvidence()
    {
        clear();
    }

    SVBreakendLowResEvidence(
        const self_t& rhs)
    {
        clear();
        merge(rhs);
    }

    const self_t&
    operator=(const self_t& rhs)
    {
        if (this == &rhs) return *this;

        clear();
        merge(rhs);
        return *this;
    }

    unsigned
    getVal(const int i) const
    {
        assert((i>=0) && (i<SVEvidenceType::SIZE));
        return _evidence[i];
    }

    unsigned
    getTotal() const
    {
        unsigned sum(0);
        for (int i(0); i<SVEvidenceType::SIZE; ++i)
        {
            sum += _evidence[i];
        }
        return sum;
    }

    void
    clear()
    {
        for (int i(0); i<SVEvidenceType::SIZE; ++i) _evidence[i] = 0;
    }

    void
    add(const int i,
        const unsigned count = 1)
    {
        assert((i>=0) && (i<SVEvidenceType::SIZE));
        _evidence[i] += count;
    }

    void
    merge(const self_t& rhs)
    {
        for (int i(0); i<SVEvidenceType::SIZE; ++i)
        {
            _evidence[i] += rhs._evidence[i];
        }
    }

private:
    unsigned short _evidence[SVEvidenceType::SIZE];
};

std::ostream&
operator<<(
    std::ostream& os,
    const SVBreakendLowResEvidence& sce);





namespace SVBreakendState
{
enum index_t
{
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
    switch (idx)
    {
    case UNKNOWN:
        return "UNKNOWN";
    case RIGHT_OPEN:
        return "RIGHT_OPEN";
    case LEFT_OPEN:
        return "LEFT_OPEN";
    case COMPLEX:
        return "COMPLEX";
    default:
        return "UNKNOWN";
    }
}

inline
bool
isSimpleBreakend(
    const index_t idx)
{
    return ((idx==RIGHT_OPEN) || (idx==LEFT_OPEN));
}

inline
bool
isSameOrientation(
    const index_t idx1,
    const index_t idx2)
{
    if (! isSimpleBreakend(idx1)) return false;
    if (! isSimpleBreakend(idx2)) return false;
    return (idx1==idx2);
}

inline
bool
isOppositeOrientation(
    const index_t idx1,
    const index_t idx2)
{
    if (! isSimpleBreakend(idx1)) return false;
    if (! isSimpleBreakend(idx2)) return false;
    return (idx1!=idx2);
}

inline
bool
isInnies(
    const bool isIdx1First,
    const index_t idx1,
    const index_t idx2)
{
    if (isIdx1First)
    {
        return ((idx1==RIGHT_OPEN) && (idx2==LEFT_OPEN));
    }
    else
    {
        return ((idx2==RIGHT_OPEN) && (idx1==LEFT_OPEN));
    }
}

inline
bool
isOutties(
    const bool isIdx1First,
    const index_t idx1,
    const index_t idx2)
{
    return isInnies((! isIdx1First),idx1,idx2);
}

}


struct SVBreakend
{
    SVBreakend() :
        state(SVBreakendState::UNKNOWN)
    {}

    bool
    isIntersect(const SVBreakend& rhs) const
    {
        if (state != rhs.state) return false;
        if (SVBreakendState::UNKNOWN == state) return true;
        return interval.isIntersect(rhs.interval);
    }

    bool
    operator<(const SVBreakend& rhs) const
    {
        if (state < rhs.state) return true;
        if (state == rhs.state)
        {
            return (interval < rhs.interval);
        }
        return false;
    }

    bool
    merge(const SVBreakend& rhs)
    {
        if (! isIntersect(rhs)) return false;
        interval.range.merge_range(rhs.interval.range);
        lowresEvidence.merge(rhs.lowresEvidence);
        return true;
    }

    void
    clear()
    {
        interval.clear();
        state=SVBreakendState::UNKNOWN;
        lowresEvidence.clear();
    }

    unsigned
    getPairCount() const
    {
        return lowresEvidence.getVal(SVEvidenceType::PAIR);
    }

    unsigned
    getLocalPairCount() const
    {
        return lowresEvidence.getVal(SVEvidenceType::LOCAL_PAIR);
    }

    /// return true if there is pair evidence, but it is only local
    bool
    isLocalPairOnly() const
    {
        return ((getLocalPairCount() > 0) && (getPairCount() == 0));
    }

    /// return true if all evidence for this breakend is local
    bool
    isLocalOnly() const
    {
        return (getLocalPairCount() == lowresEvidence.getTotal());
    }

    unsigned
    getAnyNonPairCount() const
    {
        using namespace SVEvidenceType;

        unsigned sum(0);
        for (int i(0); i< SVEvidenceType::SIZE; ++i)
        {
            if (i == PAIR) continue;
            if (i == LOCAL_PAIR) continue;
            if (i == UNKNOWN) continue;
            sum += lowresEvidence.getVal(i);
        }
        return sum;
    }

    // include any evidence type which defines a two-region hypothesis
    unsigned
    getSpanningCount() const
    {
        using namespace SVEvidenceType;

        return (lowresEvidence.getVal(PAIR) +
                lowresEvidence.getVal(CIGAR) +
                lowresEvidence.getVal(SPLIT_ALIGN));
    }

public:

    // if ! isPrecise() the interval is the X% confidence interval of the SV breakend, the interface allows for
    // various probability distributions to back this interval, but these must be accessed via
    // SVCandidate:
    GenomeInterval interval;
    SVBreakendState::index_t state;

    SVBreakendLowResEvidence lowresEvidence;
};

std::ostream&
operator<<(std::ostream& os, const SVBreakend& svb);
