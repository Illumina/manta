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


/// options shared by multiple scoring schemes:
///
/// Note that in theory these could be offered once for each scoring scheme, but
/// it would be difficult to do this efficiently because these options have an impact
/// on early scoring likelihoods.
///
struct CallOptionsShared
{
    /// This influences alignments to the ref allele when comparing ref vs alt align quality
    float snpPrior = 1e-3f;
};
