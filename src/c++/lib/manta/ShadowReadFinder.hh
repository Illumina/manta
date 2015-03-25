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
/// \author Ole Schulz-Trieglaff
///

#pragma once

#include "htsapi/bam_record.hh"

#include <string>


/// encapsulates the logic of checking for shadow reads assuming that they've been placed
/// consecutively after their mapped mate read
///
/// basic usage:
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
        _isLeftDefault(isSearchForLeftOpen),
        _isRightDefault(isSearchForRightOpen),
        _isLastSet(false),
        _lastMapq(0)
    {}

    /// reset mate tracking info
    void
    reset()
    {
        _isLastSet=false;
    }

    /// all in one single method interface to shadow finder
    ///
    /// if this is called only once for each read it will return
    /// true for any unmapped shadow read, assuming the common convention
    /// that unmapped shadows follow their anchor.
    ///
    bool
    check(const bam_record& bamRead)
    {
        if (isShadow(bamRead)) return true;
        if (isShadowAnchor(bamRead)) setAnchor(bamRead);
        return false;
    }

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

    /// the following methods are subcomponents of the check() system above --
    /// you probably only want to use one or the other

    /// check for shadow anchor status
    ///
    /// uses default left-open, right-open values
    bool
    isShadowAnchor(
        const bam_record& bamRead) const
    {
        return isShadowAnchor(bamRead,_isLeftDefault,_isRightDefault);
    }

    /// check for shadow anchor status
    ///
    bool
    isShadowAnchor(
        const bam_record& bamRead,
        const bool isSearchForLeftOpen,
        const bool isSearchForRightOpen) const;

    void
    setAnchor(
        const bam_record& bamRead);

    bool
    isShadow(
        const bam_record& bamRead);


private:

    const unsigned _minMapq;
    const bool _isLeftDefault;
    const bool _isRightDefault;
    bool _isLastSet;
    uint8_t _lastMapq;
    std::string _lastQname;
};
