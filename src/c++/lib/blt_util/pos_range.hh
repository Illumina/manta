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

    bool is_begin_pos;
    bool is_end_pos;
    pos_t begin_pos;
    pos_t end_pos;
};


/// \brief pos_range for bounded intervals only
///
struct known_pos_range : public pos_range {

    known_pos_range(const pos_t bp,const pos_t ep) : pos_range(bp,ep) {}

private:
    void clear();
};


std::ostream& operator<<(std::ostream& os, const pos_range& pr);
