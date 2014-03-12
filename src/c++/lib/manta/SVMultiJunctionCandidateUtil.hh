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

#include "manta/SVMultiJunctionCandidate.hh"
#include "manta/SVCandidateUtil.hh"

/// complex in this case means that we have no specific hypothesis for the SV --
/// it is just a single genomic region for which we schedule local assembly
///
inline
bool
isComplexSV(const SVMultiJunctionCandidate& mjSV)
{
    if (mjSV.junction.size() != 1) return false;

    return isComplexSV(mjSV.junction[0]);
}
