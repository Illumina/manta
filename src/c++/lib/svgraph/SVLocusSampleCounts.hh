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

#include "manta/SVBreakend.hh"

#include <algorithm>
#include <iosfwd>


/// enumerate evidence type estimated on input for each sample
struct SampleReadInputCounts
{
    void
    clear()
    {
        anom = 0;
        assm = 0;
        nonAnom = 0;

        remoteRecoveryCandidates = 0;
    }

    unsigned long
    total() const
    {
        return (anom+assm+nonAnom);
    }

    void
    merge(
        const SampleReadInputCounts& srs)
    {
        anom += srs.anom;
        assm += srs.assm;
        nonAnom += srs.nonAnom;
        remoteRecoveryCandidates += srs.remoteRecoveryCandidates;
    }

    void
    write(
        std::ostream& os,
        const char* label) const;


    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& anom& assm& nonAnom & remoteRecoveryCandidates;
    }

    ///< total number of non-filtered anomalous reads scanned
    unsigned long anom = 0;

    ///< total number of non-filtered non-anomolous assembly reads scanned
    unsigned long assm = 0;

    ///< total number of non-filtered non-anomalous reads scanned
    unsigned long nonAnom = 0;

    ///< subset of anom. these are reads which qualify as candidates for remote recovery
    unsigned long remoteRecoveryCandidates = 0;
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
        std::ostream& os,
        const char* label) const;


    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& eType& closeCount;
    }

    // (don't want to bother with std::array even thought size is known at compile-time:
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
        const char* label) const
    {
        input.write(os, label);
        evidence.write(os, label);
    }


    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& input& evidence;
    }

    SampleReadInputCounts input;
    SampleEvidenceCounts evidence;
};
