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

#include "boost/archive/text_oarchive.hpp"
#include "boost/archive/text_iarchive.hpp"
#include "boost/serialization/string.hpp"

#include <fstream>
#include <iostream>
#include <sstream>


// serialization
void
ReadGroupStatsSet::
save(const char* filename) const
{
    assert(NULL != filename);
    std::ofstream ofs(filename, std::ios::binary);
    boost::archive::text_oarchive oa(ofs);

    const unsigned numGroups(_group.size());
    oa << numGroups;
    for (unsigned i(0); i<numGroups; ++i)
    {
        oa << _group.get_key(i);
        oa << getStats(i);
    }
}

// restore from serialization
void
ReadGroupStatsSet::
load(const char* filename)
{
    clear();

    assert(NULL != filename);
    std::ifstream ifs(filename, std::ios::binary);
    boost::archive::text_iarchive ia(ifs);

    int numGroups;
    ia >> numGroups;
    for (int i=0; i<numGroups; i++)
    {
        std::string bamFile;
        ReadGroupStats rgs;
        ia >> bamFile;
        ia >> rgs;

        setStats(bamFile, rgs);

        std::cerr<<"i= "<<i<<"\n";
        std::cerr<<"bamFile: "<<bamFile<<"\n";
    }
}


void
ReadGroupStatsSet::
write(std::ostream& os) const
{
    const unsigned n_groups(_group.size());
    for (unsigned i(0); i<n_groups; ++i)
    {
        os << "# Bam_Path\t" << i << "\t" << _group.get_key(i) << '\n';
    }
    // write column header for better readability
    os << "# index"
       << "\treadOrientation"
       << "\tsample-count\tnumber-of-fragment-sizes"
       << '\n';

    for (unsigned i(0); i<n_groups; ++i)
    {
        os << i << '\t';
        getStats(i).write(os);
        os << '\n';
    }
}

