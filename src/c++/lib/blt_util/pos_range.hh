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

#include "blt_util/blt_types.hh"

#include <algorithm>
#include <iosfwd>


/// \brief integer ranges which are potentially unbounded
///
/// Object handles representation, including intersection with positions
/// and other ranges.
///
/// note coding convention for all ranges '_pos fields' is:
/// XXX_begin_pos is zero-indexed position at the beginning of the range
/// XXX_end_pos is zero-index position 1 step after the end of the range
///
/// any non-range pos value is assumed to be zero-indexed
///
struct pos_range {

    pos_range() : is_begin_pos(false), is_end_pos(false), begin_pos(0), end_pos(0) {}

    pos_range(const pos_t bp,const pos_t ep)
        :  is_begin_pos(true), is_end_pos(true), begin_pos(bp), end_pos(ep) {}

    void
    clear() {
        is_begin_pos=false;
        is_end_pos=false;
        begin_pos=0;
        end_pos=0;
    }

    void
    set_begin_pos(const pos_t pos) {
        begin_pos=pos;
        is_begin_pos=true;
    }

    void
    set_end_pos(const pos_t pos) {
        end_pos=pos;
        is_end_pos=true;
    }

    void
    set_range(const pos_t begin,
              const pos_t end) {
        set_begin_pos(begin);
        set_end_pos(end);
    }

    bool
    is_empty() const {
        return ! (is_begin_pos || is_end_pos);
    }

    bool
    is_complete() const {
        return (is_begin_pos && is_end_pos);
    }

    inline
    bool
    is_pos_intersect(const pos_t pos) const {

        return (((! is_begin_pos) || (pos >= begin_pos)) &&
                ((! is_end_pos) || (pos < end_pos)));
    }

    bool
    is_range_intersect(const pos_range& pr) const {
        return (((! pr.is_end_pos) || (! is_begin_pos) || (pr.end_pos > begin_pos)) &&
                ((! pr.is_begin_pos) || (! is_end_pos) || (pr.begin_pos < end_pos)));
    }

    /// does this range completely overlap pr?
    bool
    is_superset_of(const pos_range& pr) const {
        return
            (((! is_end_pos) ||
              ( pr.is_end_pos && (pr.end_pos <= end_pos) )) &&
             ((! is_begin_pos) ||
              ( pr.is_begin_pos && (pr.begin_pos >= begin_pos) )));
    }

    unsigned
    size() const {
        if (! is_complete()) return 0;
        return std::max(0,end_pos-begin_pos);
    }

    bool
    operator<(const pos_range& rhs) const
    {
        if     ((!is_begin_pos) && rhs.is_begin_pos) return true;
        else if ((is_begin_pos) && (!rhs.is_begin_pos)) return false;
        else if (is_begin_pos && rhs.is_begin_pos)
        {
            if (begin_pos < rhs.begin_pos) return true;
            if (begin_pos > rhs.begin_pos) return false;
        }

        if     ((!is_end_pos) && rhs.is_end_pos) return true;
        else if (is_end_pos && rhs.is_end_pos)
        {
            if (end_pos < rhs.end_pos) return true;
        }

        return false;
    }

    bool is_begin_pos;
    bool is_end_pos;
    pos_t begin_pos;
    pos_t end_pos;
};



/// \brief pos_range for bounded intervals only
///
struct known_pos_range : public pos_range {

    known_pos_range(const pos_t bp,const pos_t ep) : pos_range(bp,ep) {}

    bool
    operator<(const pos_range& rhs) const
    {
        if (begin_pos < rhs.begin_pos) return true;
        if (begin_pos == rhs.begin_pos)
        {
            if (end_pos < rhs.end_pos) return true;
        }
        return false;
    }

    bool
    operator==(const pos_range& rhs) const
    {
        return ((begin_pos==rhs.begin_pos) && (end_pos==rhs.end_pos));
    }

    // expand range to extend of a second range:
    void
    merge_range(const known_pos_range& kpr)
    {
        if (kpr.begin_pos<begin_pos) begin_pos=kpr.begin_pos;
        if (kpr.end_pos>end_pos) end_pos=kpr.end_pos;
    }

private:
    void clear();
};


std::ostream& operator<<(std::ostream& os, const pos_range& pr);
