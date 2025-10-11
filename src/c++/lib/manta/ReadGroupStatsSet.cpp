//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2025 Illumina, Inc.
//
// This program is licensed under the terms of the Polyform strict license
//
// ***As far as the law allows, the software comes as is, without
// any warranty or condition, and the licensor will not be liable
// to you for any damages arising out of these terms or the use
// or nature of the software, under any kind of legal claim.***
//
// You should have received a copy of the PolyForm Strict License 1.0.0
// along with this program.  If not, see <https://polyformproject.org/licenses/strict/1.0.0>.
//
//

#include "ReadGroupStatsSet.hpp"

#include "blt_util/io_util.hpp"
#include "blt_util/log.hpp"
#include "blt_util/parse_util.hpp"
#include "blt_util/string_util.hpp"

// workaround intel compiler boost warnings:
#include "boost/config.hpp"
#ifdef BOOST_INTEL_CXX_VERSION
#pragma warning push
#pragma warning(disable : 1944)
#endif

#include "boost/archive/xml_iarchive.hpp"
#include "boost/archive/xml_oarchive.hpp"

#ifdef BOOST_INTEL_CXX_VERSION
#pragma warning pop
#endif

#include "boost/serialization/map.hpp"
#include "boost/serialization/string.hpp"
#include "boost/serialization/vector.hpp"

#include <fstream>
#include <iostream>
#include <sstream>

/// this struct exists for the sole purpose of xml output:
struct ReadGroupStatsExporter {
  template <class Archive>
  void serialize(Archive& ar, const unsigned /*version*/)
  {
#ifdef READ_GROUPS
    ar& boost::serialization::make_nvp("bamFile", bamFile);
    ar& boost::serialization::make_nvp("readGroup", readGroup);
#else
    ar& boost::serialization::make_nvp("groupLabel", bamFile);
#endif
    ar& boost::serialization::make_nvp("groupStats", groupStats);
  }

  std::string    bamFile;
  std::string    readGroup;
  ReadGroupStats groupStats;
};

BOOST_CLASS_IMPLEMENTATION(ReadGroupStatsExporter, boost::serialization::object_serializable)

void ReadGroupStatsSet::merge(const ReadGroupStatsSet& rhs)
{
  const unsigned numGroups(rhs.size());
  for (unsigned i(0); i < numGroups; ++i) {
    const ReadGroupLabel& mkey(rhs.getKey(i));
    if (_group.test_key(mkey)) {
      log_os << "Can't merge stats set objects with repeated key: '" << mkey << "'";
      exit(EXIT_FAILURE);
    }

    setStats(mkey, rhs.getStats(i));
  }
}

void ReadGroupStatsSet::save(const char* filename) const
{
  assert(nullptr != filename);
  std::ofstream                ofs(filename);
  boost::archive::xml_oarchive oa(ofs);

  const unsigned numGroups(size());
  oa << boost::serialization::make_nvp("numGroups", numGroups);
  ReadGroupStatsExporter se;
  for (unsigned i(0); i < numGroups; ++i) {
    const KeyType& key(getKey(i));
    se.bamFile    = key.bamLabel;
    se.readGroup  = key.rgLabel;
    se.groupStats = getStats(i);

    std::ostringstream oss;
    oss << "groupStats_" << i;
    oa << boost::serialization::make_nvp(oss.str().c_str(), se);
  }
}

void ReadGroupStatsSet::load(const char* filename)
{
  clear();

  assert(nullptr != filename);
  std::ifstream                ifs(filename);
  boost::archive::xml_iarchive ia(ifs);

  int numGroups;
  ia >> boost::serialization::make_nvp("numGroups", numGroups);
  ReadGroupStatsExporter se;
  for (int i = 0; i < numGroups; i++) {
    ia >> boost::serialization::make_nvp("bogus", se);

    setStats(KeyType(se.bamFile.c_str(), se.readGroup.c_str()), se.groupStats);
  }
}
