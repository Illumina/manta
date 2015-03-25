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

#include "manta/SVBreakend.hh"
#include "manta/SVLocusEvidenceCount.hh"

#include <algorithm>
#include <iosfwd>


/// enumerate evidence type estimated on input for each sample
struct SampleReadInputCounts
{
    void
    clear()
    {
        minMapq = 0;
        evidenceCount.clear();
    }

    double
    total() const
    {
        return (minMapq+evidenceCount.total);
    }

    void
    merge(
        const SampleReadInputCounts& rhs)
    {
        minMapq += rhs.minMapq;
        evidenceCount.merge(rhs.evidenceCount);
    }

    void
    write(
        std::ostream& os) const;


    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& minMapq& evidenceCount;
    }

    // using doubles for integral counts here because (1) counts are potentially very high and (2) exact counts don't matter

    ///< total number of reads filtered for mapq before any classification step
    double minMapq = 0;

    SVLocusEvidenceCount evidenceCount;
};



/// enumerate detailed evidence type counts for each sample
struct SampleEvidenceCounts
{
    void
    clear()
    {
        std::fill(eType.begin(),eType.end(),0);
        closeCount = 0;
    }

    void
    merge(
        const SampleEvidenceCounts& srs)
    {
        for (unsigned i(0); i< SVEvidenceType::SIZE; ++i)
        {
            eType[i] += srs.eType[i];
        }
        closeCount += srs.closeCount;
    }

    void
    write(
        std::ostream& os) const;


    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& eType& closeCount;
    }

    // (don't want to bother with std::array even though size is known at compile-time:
    std::vector<unsigned long> eType = std::vector<unsigned long>(SVEvidenceType::SIZE,0);

    /// these are anomalous pairs which still are close to the proper pair threshold, thus downweighted
    unsigned long closeCount = 0;
};



/// total statistics for each sample
struct SampleCounts
{
    void
    clear()
    {
        input.clear();
        evidence.clear();
    }

    void
    merge(
        const SampleCounts& srs)
    {
        input.merge(srs.input);
        evidence.merge(srs.evidence);
    }

    void
    write(
        std::ostream& os,
        const char* label) const;

    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& input& evidence;
    }

    SampleReadInputCounts input;
    SampleEvidenceCounts evidence;
};



struct AllCounts
{
    void
    clear()
    {
        tumor.clear();
        normal.clear();
    }

    SampleCounts&
    getSample(
        const bool isTumor)
    {
        return ( isTumor ? tumor : normal );
    }

    const SampleCounts&
    getSample(
        const bool isTumor) const
    {
        return ( isTumor ? tumor : normal );
    }

    SampleCounts normal;
    SampleCounts tumor;
};
