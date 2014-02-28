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

#include <cassert>
#include <cstdlib>
#include <cstring>

#include <iosfwd>


struct ReadGroupLabel
{
    /// if isCopyPtrs then the strings are copied and alloced/de-alloced by
    /// the object, if false the client is responsible these pointers over
    /// the lifetime of the label:
    explicit
    ReadGroupLabel(
        const char* bamLabelInit,
        const char* rgLabelInit,
        const bool isCopyPtrsInit = true) :
        isCopyPtrs(isCopyPtrsInit),
        bamLabel((isCopyPtrs && (NULL != bamLabelInit)) ? strdup(bamLabelInit) : bamLabelInit),
        rgLabel((isCopyPtrs && (NULL != rgLabelInit)) ? strdup(rgLabelInit) : rgLabelInit)
    {
        assert(NULL != bamLabel);
        assert(NULL != rgLabel);
    }

    ReadGroupLabel(const ReadGroupLabel& rhs) :
        isCopyPtrs(rhs.isCopyPtrs),
        bamLabel(isCopyPtrs ? strdup(rhs.bamLabel) : rhs.bamLabel),
        rgLabel(isCopyPtrs ? strdup(rhs.rgLabel) : rhs.rgLabel)
    {}

    ReadGroupLabel&
    operator=(const ReadGroupLabel& rhs)
    {
        if (this == &rhs) return *this;
        clear();
        isCopyPtrs = rhs.isCopyPtrs;
        bamLabel = (isCopyPtrs ? strdup(rhs.bamLabel) : rhs.bamLabel);
        rgLabel = (isCopyPtrs ? strdup(rhs.rgLabel) : rhs.rgLabel);
        return *this;
    }

public:

    ~ReadGroupLabel()
    {
        clear();
    }

    /// sort allowing for NULL string pointers in primary and secondary key:
    bool
    operator<(
        const ReadGroupLabel& rhs) const
    {
        const int scval(strcmp(bamLabel,rhs.bamLabel));
        if (scval < 0) return true;
        if (scval == 0)
        {
            return (strcmp(rgLabel,rhs.rgLabel) < 0);
        }

        return false;
    }

private:
    void
    clear()
    {
        if (isCopyPtrs)
        {
            if (NULL != bamLabel) free(const_cast<char*>(bamLabel));
            if (NULL != rgLabel) free(const_cast<char*>(rgLabel));
        }
    }

    bool isCopyPtrs;

public:
    const char* bamLabel;
    const char* rgLabel;
};


std::ostream&
operator<<(
    std::ostream& os,
    const ReadGroupLabel& rgl);
