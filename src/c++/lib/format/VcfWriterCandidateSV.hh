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

#include "format/VcfWriterSV.hh"


struct VcfWriterCandidateSV : public VcfWriterSV
{
    VcfWriterCandidateSV(
        const std::string& referenceFilename,
        const bool isRNA,
        const SVLocusSet& set,
        std::ostream& os) :
        VcfWriterSV(referenceFilename, isRNA, set, os)
    {}

    void
    addHeaderInfo() const override;

    void
    modifyTranslocInfo(
        const SVCandidate& sv,
        const bool isFirstOfPair,
        InfoTag_t& infoTags) const override;

    void
    modifyInvdelInfo(
        const SVCandidate& sv,
        const bool isBp1First,
        InfoTag_t& infoTags) const override;

    void
    writeSV(
        const SVCandidateSetData& svData,
        const SVCandidateAssemblyData& adata,
        const SVCandidate& sv,
        const SVId& svId);
};

