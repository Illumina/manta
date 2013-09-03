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

#include "ReadGroupStatsSet.hh"

#include "blt_util/log.hh"
#include "blt_util/io_util.hh"
#include "blt_util/parse_util.hh"
#include "blt_util/string_util.hh"

#include "boost/archive/xml_oarchive.hpp"
#include "boost/archive/xml_iarchive.hpp"
#include "boost/serialization/map.hpp"
#include "boost/serialization/string.hpp"
#include "boost/serialization/vector.hpp"

#include <fstream>
#include <iostream>
#include <sstream>


/// this struct exists for the sole purpose of xml output:
struct ReadGroupStatsExporter
{

    template<class Archive>
    void serialize(Archive& ar, const unsigned /*version*/)
    {
        ar& boost::serialization::make_nvp("groupLabel", groupLabel);
        ar& boost::serialization::make_nvp("groupStats", groupStats);
    }

    std::string groupLabel;
    ReadGroupStats groupStats;
};

BOOST_CLASS_IMPLEMENTATION(ReadGroupStatsExporter, boost::serialization::object_serializable)




// serialization
void
ReadGroupStatsSet::
save(const char* filename) const
{
    assert(NULL != filename);
    std::ofstream ofs(filename);
    boost::archive::xml_oarchive oa(ofs);

    const unsigned numGroups(_group.size());
    oa << boost::serialization::make_nvp("numGroups", numGroups);
    ReadGroupStatsExporter se;
    for (unsigned i(0); i<numGroups; ++i)
    {
        se.groupLabel = _group.get_key(i);
        se.groupStats = getStats(i);

        std::ostringstream oss;
        oss << "groupStats_" << i;
        oa << boost::serialization::make_nvp(oss.str().c_str(), se);
    }
}

// restore from serialization
void
ReadGroupStatsSet::
load(const char* filename)
{
    clear();

    assert(NULL != filename);
    std::ifstream ifs(filename);
    boost::archive::xml_iarchive ia(ifs);

    int numGroups;
    ia >> boost::serialization::make_nvp("numGroups", numGroups);
    ReadGroupStatsExporter se;
    for (int i=0; i<numGroups; i++)
    {
        ia >> boost::serialization::make_nvp("bogus", se);

        setStats(se.groupLabel, se.groupStats);
    }
}


