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

#pragma once

#include "blt_util/blt_types.hh"

#include <string>


/// Manages a partial reference sequence segment
///
/// This object holds the reference sequence specified by the current
/// runs begin and end range, plus some padding on each side. To get
/// this integrated into the current code as quickly as possible it
/// currently exposes the internal string object holding the sequence
/// data. When time allows this will be restricted so that a compressed
/// internal object can be used.
///
struct reference_contig_segment
{

    reference_contig_segment()
        : _offset(0)
    {}

    char
    get_base(const pos_t pos) const
    {
        if (pos<_offset || pos>=end()) return 'N';
        return _seq[pos-_offset];
    }

    void
    get_substring(const pos_t pos,
                  const pos_t length,
                  std::string& substr) const
    {

        if (pos<_offset || (pos+length)>end())
        {
            //slow path (minority of calls):
            substr.clear();
            for (int i(0); i<length; ++i)
            {
                substr.push_back(get_base(pos+i));
            }
        }
        else
        {
            //fast path
            substr.assign(_seq,pos-_offset,length);
        }
    }

    std::string& seq()
    {
        return _seq;
    }
    const std::string& seq() const
    {
        return _seq;
    }

    pos_t
    get_offset() const
    {
        return _offset;
    }

    void
    set_offset(const pos_t offset)
    {
        _offset=offset;
    }

    pos_t
    end() const
    {
        return _offset+_seq.size();
    }

private:

    pos_t _offset;
    std::string _seq;
};


