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

#pragma once

#include <iosfwd>
#include <string>

#include "boost/optional.hpp"

#include "blt_util/id_map.hpp"
#include "manta/ReadGroupLabel.hpp"
#include "manta/ReadGroupStats.hpp"

/// \brief manages multiple read_group_stats
///
struct ReadGroupStatsSet {
  typedef ReadGroupLabel KeyType;

  bool empty() const { return _group.empty(); }

  unsigned size() const { return _group.size(); }

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
  boost::optional<unsigned> getGroupIndex(const ReadGroupLabel& rgLabel) const
  {
    return _group.get_optional_id(rgLabel);
  }

  /// get stats associated with index
  const ReadGroupStats& getStats(const unsigned groupIndex) const { return _group.get_value(groupIndex); }

  const KeyType& getKey(const unsigned groupIndex) const { return _group.get_key(groupIndex); }

  /// set stats for index
  void setStats(const ReadGroupLabel& rgLabel, const ReadGroupStats& rps) { _group.insert(rgLabel, rps); }

  /// merge in the contents of another stats set object:
  void merge(const ReadGroupStatsSet& rhs);

  void save(const char* filename) const;

  void load(const char* filename);

  bool isEmpty() { return _group.empty(); }

private:
  void clear() { _group.clear(); }

  id_map<KeyType, ReadGroupStats> _group;
};
