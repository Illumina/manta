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

#include "blt_util/align_path.hh"
#include "manta/SVBreakend.hh"

#include <iosfwd>
#include <string>


struct SVCandidate
{
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
            fwReads += rhs.fwReads;
            rvReads += rhs.rvReads;
        }
        else
        {
            bp1.merge(rhs.bp2);
            bp2.merge(rhs.bp1);
            fwReads += rhs.rvReads;
            rvReads += rhs.fwReads;
        }

        _isImprecise = (isImprecise() || rhs.isImprecise());

        return true;
    }

    bool
    evidenceMerge(const SVCandidate& rhs)
    {
        if (! isIntersect(rhs)) return false;

        if (bp1.isIntersect(rhs.bp1))
        {
            bp1.evidenceMerge(rhs.bp1);
            bp2.evidenceMerge(rhs.bp2);
        }
        else
        {
            bp1.evidenceMerge(rhs.bp2);
            bp2.evidenceMerge(rhs.bp1);
        }

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
        isUnknownSizeInsertion = false;
        unknownSizeInsertionLeftSeq.clear();
        unknownSizeInsertionRightSeq.clear();
    }

    void
    setPrecise()
    {
        _isImprecise = false;
    }

    bool
    isForward() const
    {
        return (fwReads > rvReads);
    }

    /// if 1 is added to the position of one breakend (within the homologous breakend range), then is 1 also added to the other breakend?
    ///
    /// if false then breakends move in opposite directions;
    bool
    isBreakendRangeSameShift() const
    {
        return (bp1.state != bp2.state);
    }

    int
    centerSize() const
    {
        return std::abs(bp2.interval.range.center_pos() - bp1.interval.range.center_pos());
    }

    /// for precise SV report the full spanning count
    /// for imprecise SV report spanning pairs only
    ///
    unsigned
    getPostAssemblySpanningCount() const
    {
        if (isImprecise())
        {
            return bp1.getPairCount();
        }
        else
        {
            return bp1.getSpanningCount();
        }
    }

private:
    bool _isImprecise = true;

public:
    SVBreakend bp1;
    SVBreakend bp2;

    // this is either a micro-insertion in a large-scale SV, or the inserted sequence of an actual insertion
    // in case bp1 and bp2 are on opposite strands (ie. an inversion) the insertSeq is oriented to the fwd strand for bp1
    std::string insertSeq;

    // for some small indels, the alignment becomes sufficiently complex that a CIGAR string provides better detail
    // (this is provided for any small SV which is more complicated than a simple insert or deletion)
    ALIGNPATH::path_t insertAlignment;

    unsigned candidateIndex = 0; ///< low-res candidate index number, used to generate unique SV id
    unsigned assemblyAlignIndex = 0; ///< high-res assembly index number of alignment, used to generate unique SV id
    unsigned assemblySegmentIndex = 0; ///< high-res assembly index number of alignment segment, used to generate unique SV id

    bool isUnknownSizeInsertion = false; ///< these insertions haven't been assembled all the way through

    std::string unknownSizeInsertionLeftSeq; ///< for an incomplete insertion, this is the known left side of the insert sequence
    std::string unknownSizeInsertionRightSeq; ///< for an incomplete insertion, this is the known right side of the insert sequence

    unsigned fwReads = 0; ///< Number of reads (pairs) supporting a direction from bp1 to bp2 (used for stranded RNA data)
    unsigned rvReads = 0; ///< Number of reads (pairs) directed from bp2 to bp1
};

std::ostream&
operator<<(std::ostream& os, const SVCandidate& svc);


namespace FRAGSOURCE
{
enum index_t
{
    UNKNOWN,
    READ1,
    READ2,
    PAIR
};
}

/// when we extract an SV candidate from a single piece of evidence, it can be treated as a special 'observation' class:
///
struct SVObservation : public SVCandidate
{
    SVObservation() :
        SVCandidate(),
        evtype(SVEvidenceType::UNKNOWN),
        fragSource(FRAGSOURCE::UNKNOWN)
    {}

    void
    clear()
    {
        evtype = SVEvidenceType::UNKNOWN;
        fragSource = FRAGSOURCE::UNKNOWN;
        SVCandidate::clear();
    }

    bool
    isSingleReadSource() const
    {
        using namespace FRAGSOURCE;
        return ((fragSource == READ1) || (fragSource == READ2));
    }

    bool
    isRead1Source() const
    {
        using namespace FRAGSOURCE;
        return (fragSource == READ1);
    }

    SVEvidenceType::index_t evtype;
    FRAGSOURCE::index_t fragSource;
};
