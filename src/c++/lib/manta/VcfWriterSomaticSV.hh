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

#include "manta/SomaticSVScoreInfo.hh"
#include "manta/VcfWriterSV.hh"
#include "options/SomaticCallOptions.hh"


struct VcfWriterSomaticSV : public VcfWriterSV
{
    VcfWriterSomaticSV(
        const SomaticCallOptions& somaticOpt,
        const std::string& referenceFilename,
        const SVLocusSet& set,
        std::ostream& os) :
        VcfWriterSV(referenceFilename,set,os),
        _somaticOpt(somaticOpt)
    {}

    void
    writeSV(
        const EdgeInfo& edge,
        const SVCandidateData& svData,
        const unsigned svIndex,
        const SVCandidate& sv,
        const SomaticSVScoreInfo& ssInfo);

private:
    const SomaticCallOptions& _somaticOpt;
};

