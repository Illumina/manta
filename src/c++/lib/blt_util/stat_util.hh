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

bool
is_chi_sqr_reject(const double xsq,
                  const unsigned df,
                  const double alpha);

bool
is_lrt_reject_null(const double null_loghood,
                   const double alt_loghood,
                   const unsigned df,
                   const double alpha);
