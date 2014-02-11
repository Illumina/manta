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

///
/// \author Chris Saunders
///

#include "common/Exceptions.hh"
#include "svgraph/SVLocusNode.hh"

#include <iostream>
#include <sstream>


std::ostream&
operator<<(std::ostream& os, const SVLocusEdge& edge)
{
    os << "Edgecount: " << edge.getCount() << " isCountExact?: " << edge.isCountExact();
    return os;
}



const SVLocusEdgesType SVLocusEdgeManager::staticMap;



void
SVLocusNode::
getEdgeException(
    const NodeIndexType toIndex,
    const char* label) const
{
    using namespace illumina::common;

    std::ostringstream oss;
    oss << "ERROR: SVLocusNode::" << label << "() no edge exists\n";
    oss << "\tfrom node: " << (*this) << "\n";
    oss << "\tto_node index: " << toIndex << "\n";
    BOOST_THROW_EXCEPTION(LogicException(oss.str()));
}



std::ostream&
operator<<(std::ostream& os, const SVLocusNode& node)
{
    os << "LocusNode: " << node.getInterval()
       << " n_edges: " << node.size()
       << " out_count: " << node.outCount()
       << " evidence: " << node.getEvidenceRange()
       << "\n";

    const SVLocusEdgeManager edgeMap(node.getEdgeManager());
    BOOST_FOREACH(const SVLocusEdgesType::value_type& edgeIter, edgeMap.getMap())
    {
        os << "\tEdgeTo: " << edgeIter.first
           << " out_count: " << edgeIter.second.getCount()
           << "\n";
    }
    return os;
}

