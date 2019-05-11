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

#include "boost/test/unit_test.hpp"

#include "svgraph/SVLocus.hpp"
#include "test/testSVLocusUtil.hpp"

BOOST_AUTO_TEST_SUITE(SVLocusNode_test_suite)

BOOST_AUTO_TEST_CASE(test_SVLocusNode_EdgeManager)
{
  // test the new edge manager for SVLocusNode
  SVLocus       locus1;
  NodeIndexType nodePtr1     = locus1.addNode(GenomeInterval(1, 10, 20));
  NodeIndexType nodePtr1copy = locus1.addNode(GenomeInterval(1, 10, 20));
  locus1.linkNodes(nodePtr1, nodePtr1copy, 1, 1);
  locus1.mergeSelfOverlap();

  BOOST_REQUIRE_EQUAL(locus1.size(), 1u);

  const SVLocusNode& node(static_cast<const SVLocus&>(locus1).getNode(0));

  const SVLocusEdgeManager em = node.getEdgeManager();

  BOOST_REQUIRE_EQUAL(em.getMap().size(), 1u);

  BOOST_REQUIRE_EQUAL(em.getMap().begin()->first, 0u);
  BOOST_REQUIRE_EQUAL(em.getMap().begin()->second.getCount(), 1u);
}

BOOST_AUTO_TEST_SUITE_END()
