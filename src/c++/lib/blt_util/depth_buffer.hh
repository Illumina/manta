
///
/// \author Chris Saunders
///

#pragma once

#include "blt_util/blt_types.hh"

#include <cassert>

#include <map>


/// simple map of position to depth
///
struct depth_buffer
{

    unsigned
    val(const pos_t pos) const
    {
        const citer i(_data.find(pos));
        if (i == _data.end()) return 0;
        else                  return i->second;
    }

    void
    inc(const pos_t pos)
    {
        const iter i(_data.find(pos));
        if (i == _data.end()) _data[pos] = 1;
        else                  i->second += 1;
    }

    void
    clear_pos(const pos_t pos)
    {
        _data.erase(pos);
    }

    /// return true if buffered depth exceeds depth in [begin,end]
    bool
    is_range_ge_than(const pos_t begin,
                     const pos_t end,
                     const unsigned depth) const
    {
        assert(begin <= end);
        citer i(_data.lower_bound(begin));
        const citer i_end(_data.upper_bound(end));
        for (; i!=i_end; ++i) if (i->second >= depth) return true;
        return false;
    }

private:
    typedef std::map<pos_t,unsigned> count_t;
    typedef count_t::iterator iter;
    typedef count_t::const_iterator citer;

    count_t _data;
};

