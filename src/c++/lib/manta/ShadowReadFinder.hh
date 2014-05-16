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

#pragma once

#include "blt_util/bam_record.hh"

#include <string>

#define BAMSURGEON_BUG_WORKAROUND


/// encapsulates the logic of checking for shadow reads assuming that they've been placed
/// consecutively after their mapped mate read
///
/// usage:
///
/// ShadowReadFinder checker(blah)
/// for bam_region in regions :
///     for bam_record in bam_region :
///         checker.check(bam_record)
///     checker.reset()
///
struct ShadowReadFinder
{
    ShadowReadFinder(
        const unsigned minMapq,
        const bool isSearchForLeftOpen = true,
        const bool isSearchForRightOpen = true) :
        _minMapq(minMapq),
        _isLeft(isSearchForLeftOpen),
        _isRight(isSearchForRightOpen),
        _isLastSet(false),
        _lastMapq(0)
    {}

    /// reset mate tracking info
    void
    reset()
    {
        _isLastSet=false;
    }

    bool
    check(const bam_record& bamRead);

    /// only valid after check() is true
    unsigned
    getMateMapq() const
    {
        return _lastMapq;
    }

    bool
    isShadowMate() const
    {
        return _isLastSet;
    }

private:

    const unsigned _minMapq;
    const bool _isLeft;
    const bool _isRight;
    bool _isLastSet;
    uint8_t _lastMapq;
    std::string _lastQname;
};
