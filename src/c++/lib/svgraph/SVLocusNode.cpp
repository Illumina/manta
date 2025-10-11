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

/// \file
/// \author Chris Saunders
///

#include "svgraph/SVLocusNode.hpp"
#include "common/Exceptions.hpp"

#include <iostream>
#include <sstream>

std::ostream& operator<<(std::ostream& os, const SVLocusEdge& edge)
{
  os << "Edgecount: " << edge.getCount() << " isCountExact?: " << edge.isCountExact();
  return os;
}

const SVLocusEdgesType SVLocusEdgeManager::staticMap;

void SVLocusNode::getEdgeException(const NodeIndexType toIndex, const char* label) const
{
  using namespace illumina::common;

  std::ostringstream oss;
  oss << "SVLocusNode::" << label << "() no edge exists\n";
  oss << "\tfrom node: " << (*this) << "\n";
  oss << "\tto_node index: " << toIndex;
  BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
}

std::ostream& operator<<(std::ostream& os, const SVLocusNode& node)
{
  os << "LocusNode: " << node.getInterval() << " n_edges: " << node.size()
     << " out_count: " << node.outCount() << " evidence: " << node.getEvidenceRange() << "\n";

  const SVLocusEdgeManager edgeMap(node.getEdgeManager());
  for (const SVLocusEdgesType::value_type& edgeIter : edgeMap.getMap()) {
    os << "\tEdgeTo: " << edgeIter.first << " out_count: " << edgeIter.second.getCount() << "\n";
  }
  return os;
}
