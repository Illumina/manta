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

#include "manta/SVCandidate.hh"


/// SVComplexCandidate represents an associated grouping of multiple breakend pairs
///
/// Examples: The two breakend pairs of a simple inversion would form a complex candidates composed of the two SVCandidates for the
/// forward and reverse junctions.
///
/// for a complex candidates we want to test the concept that the full set of breakends occurred together as part of the same
/// event, and thus be able to call the full event with perhaps less evidence per each single breakend than would be acceptable
/// during regular calling.
///
struct SVMultiJunctionCandidate
{
    SVMultiJunctionCandidate()
    {}

    std::vector<SVCandidate> junction;

    /// TODO: need to design a quick data structure to iterate through complex event breakend regions and pull out the associated candidate's
    /// breakends -- is this just a mapping from MJ breakend 'groups' to literal junction breakends??
};

std::ostream&
operator<<(std::ostream& os, const SVMultiJunctionCandidate& svc);
