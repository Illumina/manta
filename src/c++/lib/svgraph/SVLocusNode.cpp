//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2019 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
