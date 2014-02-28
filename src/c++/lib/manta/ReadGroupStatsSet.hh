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

#pragma once

#include "blt_util/id_map.hh"
#include "manta/ReadGroupLabel.hh"
#include "manta/ReadGroupStats.hh"

#include "boost/optional.hpp"

#include <iosfwd>
#include <string>


/// \brief manages multiple read_group_stats
///
struct ReadGroupStatsSet
{
    typedef ReadGroupLabel KeyType;

    bool
    empty() const
    {
        return _group.empty();
    }

    unsigned
    size() const
    {
        return _group.size();
    }

    /// \brief get the index of a read group
    ///
    /// the index can be used for fast lookup of the
    /// stats for that group
    ///
    /// if the group does not exist, the returned value
    /// evaluates to false per boost::optional
    ///
    /// Each read group is identified as a combination of a bam filename and
    /// an RG tag label. An empty label refers to the "default" read group
    /// for the file (all records that had no RG tag).
    boost::optional<unsigned>
    getGroupIndex(
            const ReadGroupLabel& rgLabel) const
    {
        return _group.get_optional_id(rgLabel);
    }

    /// get stats associated with index
    const ReadGroupStats&
    getStats(
        const unsigned groupIndex) const
    {
        return _group.get_value(groupIndex);
    }

    const KeyType&
    getKey(
        const unsigned groupIndex) const
    {
        return _group.get_key(groupIndex);
    }

    /// set stats for index
    void
    setStats(
        const ReadGroupLabel& rgLabel,
        const ReadGroupStats& rps)
    {
        _group.insert(rgLabel,rps);
    }

    /// serialize
    void
    save(
        const char* filename) const;

    /// deserialize
    void
    load(
        const char* filename);

private:
    void
    clear()
    {
        _group.clear();
    }

    id_map<KeyType, ReadGroupStats> _group;
};

