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
/// \author Ole Schulz-Trieglaff
///

#include "ShadowReadFinder.hh"

#include "blt_util/bam_record_util.hh"


#ifdef DEBUG_IS_SHADOW
#include <iostream>
#include "blt_util/log.hh"
#endif


static
bool
isGoodShadow(
    const bam_record& bamRead,
    const std::string& lastQname)
{
#ifdef DEBUG_IS_SHADOW
    static const std::string logtag("isGoodShadow");
#endif

    // sanity check that this is a shadow read:
    if (!bamRead.is_unmapped()) return false;
    if (bamRead.is_mate_unmapped()) return false;

    static const unsigned minAvgQualShadow = 25;
    if (get_avg_quality(bamRead) < minAvgQualShadow)
    {
        return false;
    }

    if (bamRead.qname() != lastQname)
    {
        // something went wrong here, shadows should have their singleton partner
        // preceding them in the BAM file.
#ifdef DEBUG_IS_SHADOW
        log_os << logtag << " ERROR: Shadow without matching singleton : " << bamRead.qname() << " vs " << lastQname << std::endl;
#endif
        return false;
    }

#ifdef DEBUG_IS_SHADOW
    log_os << logtag << " Found shadow!\n";
            << logtag << " this mapq  = " << ((unsigned int)bamRead.map_qual()) << std::endl;
            << logtag << " last qname = " << lastQname << std::endl;
#endif

    return true;
}



bool
ShadowReadFinder::
check(
    const bam_record& bamRead)
{
    if (_isLastSet)
    {
        _isLastSet = false;
        if (isGoodShadow(bamRead, _lastQname)) return true;
    }
    else if ((! bamRead.is_unmapped()) && (bamRead.is_mate_unmapped()))
    {

        if ((! _isLeft) && (! bamRead.is_fwd_strand())) return false;
        if ((! _isRight) && bamRead.is_fwd_strand()) return false;

        _lastMapq  = bamRead.map_qual();
        if (_lastMapq < _minMapq) return false;

        _lastQname = bamRead.qname();
        _isLastSet = true;
    }

    return false;
}
