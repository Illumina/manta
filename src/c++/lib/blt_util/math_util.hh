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

/// \file

/// \author Chris Saunders
///
#ifndef __MATH_UTIL_HH
#define __MATH_UTIL_HH

#include "boost/math/special_functions/log1p.hpp"

#include <cmath>

#include <algorithm>


/// returns log(1+x), switches to special libc function when abs(x) is small
///
template <typename FloatType>
FloatType
log1p_switch(const FloatType x) {

    // better number??
    static const FloatType smallx_thresh(0.01);

    if (std::abs(x)<smallx_thresh) {
        return boost::math::log1p(x);
    } else {
        return std::log(1+x);
    }
}


/// returns equiv of log(exp(x1)+exp(x2))
///
template <typename FloatType>
FloatType
log_sum(FloatType x1, FloatType x2) {
    if (x1<x2) std::swap(x1,x2);
    return x1 + log1p_switch(std::exp(x2-x1));
}

#endif
