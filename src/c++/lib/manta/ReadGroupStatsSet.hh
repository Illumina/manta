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

#pragma once

#include "blt_util/id_map.hh"
#include "manta/ReadGroupStats.hh"

#include "boost/optional.hpp"

#include <iosfwd>
#include <string>


/// \brief manages multiple read_group_stats
///
struct ReadGroupStatsSet
{

    /// \brief get the index of a read group
    ///
    /// the index can be used for fast lookup of the
    /// stats for that group
    ///
    /// if the group does not exist, the returned value
    /// evaluates to false per boost::optional
    ///
    /// for now, a "read group" is fixed to the name of
    /// each bam file
    boost::optional<unsigned>
    getGroupIndex(const std::string& bam_file) const
    {
        return _group.get_optional_id(bam_file);
    }

    /// get stats associated with index
    const ReadGroupStats&
    getStats(const unsigned group_index) const
    {
        return _group.get_value(group_index);
    }

    /// set stats for index
    void
    setStats(const std::string& bam_file,
             const ReadGroupStats& rps)
    {
        _group.insert(bam_file,rps);
    }


    /// serialize
    void
    save(const char* filename) const;

    /// restore from serialization
    void
    load(const char* filename);

#if 0
    // write out brief info of the stats class
    // and some debugging info if under debug mode
    void
    write(std::ostream& os) const;
#endif

private:
    void
    clear()
    {
        _group.clear();
    }

    id_map<std::string, ReadGroupStats> _group;
};

