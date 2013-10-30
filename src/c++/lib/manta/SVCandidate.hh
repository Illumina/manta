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
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#pragma once

#include "blt_util/align_path.hh"
#include "manta/SVBreakend.hh"

#include <iosfwd>
#include <string>


struct SVCandidate
{
    SVCandidate() :
        _isImprecise(true),
        candidateIndex(0),
        assemblyAlignIndex(0),
        assemblySegmentIndex(0)
    {}

#if 0
    double
    breakpointProb(pos_t x, pos_t y) const;
#endif


    /// if false, the breakend interval is at base-pair resolution
    ///
    /// false does not mean that the interval size is zero, a precise breakend interval range represents microhomology at the breakend site
    bool
    isImprecise() const
    {
        return _isImprecise;
    }

    bool
    isIntersect(const SVCandidate& rhs) const
    {
        return ((bp1.isIntersect(rhs.bp1) && bp2.isIntersect(rhs.bp2)) ||
                (bp1.isIntersect(rhs.bp2) && bp2.isIntersect(rhs.bp1)));
    }

    bool
    merge(const SVCandidate& rhs)
    {
        if (! isIntersect(rhs)) return false;

        if (bp1.isIntersect(rhs.bp1))
        {
            bp1.merge(rhs.bp1);
            bp2.merge(rhs.bp2);
        }
        else
        {
            bp1.merge(rhs.bp2);
            bp2.merge(rhs.bp1);
        }

        _isImprecise = (isImprecise() || rhs.isImprecise());

        return true;
    }

    void
    clear()
    {
        _isImprecise = true;
        bp1.clear();
        bp2.clear();
        insertSeq.clear();
        candidateIndex=0;
        assemblyAlignIndex=0;
        assemblySegmentIndex=0;
    }

    void
    setPrecise()
    {
        _isImprecise = false;
    }

    /// if 1 is added to the position of one breakend (within the homologous breakend range), then is 1 also added to the other breakend?
    ///
    /// if false then breakends move in opposite directions;
    bool
    isBreakendRangeSameShift() const
    {
        return (bp1.state != bp2.state);
    }

#if 0
    std::pair<const SVBreakend&, const SVBreakend&>
    getOrderedBreakends() const
    {
        if (bp2 < bp1)
        {
            return std::make_pair(bp2,bp1);
        }
        else
        {
            return std::make_pair(bp1,bp2);
        }
    }
#endif

private:
    bool _isImprecise;

public:
    SVBreakend bp1;
    SVBreakend bp2;

    // this is either a micro-insertion in a large-scale SV, or the inserted sequence of an actual insertion
    // in case bp1 and bp2 are on opposite strands (ie. an inversion) the insertSeq is oriented to the fwd strand for bp1
    std::string insertSeq;

    // for some small indels, the alignment becomes sufficiently complex that a CIGAR string provides better detail
    // (this is provided for any small SV which is more complicated than a simple insert or deletion)
    ALIGNPATH::path_t insertAlignment;

    unsigned candidateIndex; // low-res candidate index number, used to generate unique SV id
    unsigned assemblyAlignIndex; // high-res assembly index number of alignment, used to generate unique SV id
    unsigned assemblySegmentIndex; // high-res assembly index number of alignment segment, used to generate unique SV id
};


std::ostream&
operator<<(std::ostream& os, const SVCandidate& svc);
