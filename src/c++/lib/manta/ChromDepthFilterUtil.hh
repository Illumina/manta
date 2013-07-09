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

#include "blt_util/bam_header_info.hh"
#include "blt_util/chrom_depth_map.hh"

#include <cassert>

#include <vector>


/// hold information about chrom depth cutoffs
///
/// preprocess the chrom depth file so that the filter value can be
/// efficiently looked up by bam tid
///
struct ChromDepthFilterUtil
{
    ChromDepthFilterUtil(
            const std::string& chromDepthFile,
            const double maxDepthFactor,
            const bam_header_info& header);

    bool
    isMaxDepthFilter() const
    {
        return _isMaxDepthFilter;
    }

    double
    maxDepth(const int32_t tid) const
    {
        assert((tid >= 0) && (tid < static_cast<int32_t>(_maxDepthFilter.size())));
        return _maxDepthFilter[tid];
    }

private:
    bool _isMaxDepthFilter;
    std::vector<double> _maxDepthFilter;
};
