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


/// enumerate evidence type estimated on input for each sample
struct SVLocusEvidenceCount
{
    void
    clear()
    {
        total = 0;
        ignored = 0;
        anom = 0;
        split = 0;
        indel = 0;
        assm = 0;
        remoteRecoveryCandidates = 0;
    }

    void
    merge(
        const SVLocusEvidenceCount& rhs)
    {
        total += rhs.total;
        ignored += rhs.ignored;
        anom += rhs.anom;
        split += rhs.split;
        indel += rhs.indel;
        assm += rhs.assm;
        remoteRecoveryCandidates += rhs.remoteRecoveryCandidates;
    }

    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& total& ignored& anom& split& indel& assm& remoteRecoveryCandidates;
    }

    // using doubles for integral counts here because (1) counts are potentially very high and (2) exact counts don't matter

    ///< total number of non-filtered anomalous reads scanned
    double total = 0;

    ///< total number of non-filtered reads ignored for SV purposes
    double ignored = 0;

    ///< total number of non-filtered anomalous reads scanned
    double anom = 0;

    ///< total number of non-filtered split (SA-tag) reads scanned
    double split = 0;

    ///< total number of non-filtered CIGAR large indel reads scanned
    double indel = 0;

    ///< total number of non-filtered semi-aligned reads scanned
    double assm = 0;

    ///< subset of anom. these are reads which qualify as candidates for remote recovery
    double remoteRecoveryCandidates = 0;
};
