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
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#include "common/Exceptions.hh"
#include "manta/ChromDepthFilterUtil.hh"

#include "boost/foreach.hpp"

#include <sstream>



ChromDepthFilterUtil::
ChromDepthFilterUtil(
    const std::string& chromDepthFile,
    const double maxDepthFactor,
    const bam_header_info& header) :
    _isMaxDepthFilter(! chromDepthFile.empty())
{
    using namespace illumina::common;

    // read in chrom depth file if one is specified:
    if (! _isMaxDepthFilter) return;

    cdmap_t chromDepth;
    parse_chrom_depth(chromDepthFile,chromDepth);

    // translate string chrom labels into tid values in lookup vector:
    //
    BOOST_FOREACH(const bam_header_info::chrom_info& cdata, header.chrom_data)
    {
        cdmap_t::const_iterator cdi(chromDepth.find(cdata.label));
        if (cdi == chromDepth.end())
        {
            std::ostringstream oss;
            oss << "ERROR: Can't find chromosome: '" << cdata.label
                << "' in chrom depth file: " << chromDepthFile << "\n";
            BOOST_THROW_EXCEPTION(LogicException(oss.str()));
        }

        _maxDepthFilter.push_back(cdi->second*maxDepthFactor);
        assert(_maxDepthFilter.back()>=0.);
    }
}

