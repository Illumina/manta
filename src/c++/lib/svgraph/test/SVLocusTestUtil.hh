//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2018 Illumina, Inc.
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

#pragma once

#include "svgraph/SVLocus.hh"
#include "boost/test/unit_test.hpp"
#include "svgraph/SVLocusSet.hh"

#include <stdio.h>
#include <fstream>
#include <string>

/// \brief Add a pair of nodes to a SVLocus object.
inline
void
locusAddPair(
    SVLocus& locus,
    const int32_t tid1,
    const int32_t beginPos1,
    const int32_t endPos1,
    const int32_t tid2,
    const int32_t beginPos2,
    const int32_t endPos2,
    const bool bothLocal = false,
    const int count = 1)
{
    const NodeIndexType nodePtr1 = locus.addNode(GenomeInterval(tid1,beginPos1,endPos1));
    const NodeIndexType nodePtr2 = locus.addNode(GenomeInterval(tid2,beginPos2,endPos2));
    const int remoteCount(bothLocal ? count : 0);
    locus.linkNodes(nodePtr1,nodePtr2,count,remoteCount);
}


/// \brief Test the size and count of the properties of the SVLocusSet.
///
/// \param set1 The SV locus set to be tested
/// \param setSize The expected number of SVLocus objects in the set
/// \param setNonEmptySize The expected number of SVLocus objects that are not empty in the set
/// \param setNodeCount The expected total number of nodes in the SVLocus objects in the set
/// \param setEdgeCount The expected total number of edges in the SVLocus objects in the set
inline
void
TestSVLocusSetProperties(
    const SVLocusSet& set1,
    int setSize,
    int setNonEmptySize,
    int setNodeCount,
    int setEdgeCount)
{
    // Test the base SVLocusSet Contents
    BOOST_REQUIRE_EQUAL(set1.size(), setSize);
    BOOST_REQUIRE_EQUAL(set1.nonEmptySize(), setNonEmptySize);
    BOOST_REQUIRE_EQUAL(set1.totalNodeCount(), setNodeCount);
    BOOST_REQUIRE_EQUAL(set1.totalEdgeCount(), setEdgeCount);
}


/// \brief Test the distribution of the total edges and node observations for all nodes in the SVLocusSet.
///
/// \param set1 The SV locus set to be tested
/// \param edgeCountDistro The expected number of edges in each SVLocus object in the set
/// \param obsCountDistro The expected number of node observations in each SVLocus object in the set
inline
void
TestSVLocusSetDistro(
    const SVLocusSet& set1,
    std::vector<unsigned> edgeCountDistro,
    std::vector<unsigned> obsCountDistro)
{

    // review merged nodes
    static const unsigned maxCount(set1.totalEdgeCount());
    std::vector<unsigned> edgeCount(maxCount);
    std::vector<unsigned> obsCount(maxCount);

    set1.getNodeEdgeCountDistro(edgeCount);
    set1.getNodeObsCountDistro(obsCount);

    for ( unsigned count = 0; count < maxCount; count++)
    {
        BOOST_REQUIRE_EQUAL(edgeCountDistro[count], edgeCount[count]);
    }

    for ( unsigned count = 0; count < maxCount; count++)
    {
        BOOST_REQUIRE_EQUAL(obsCountDistro[count], obsCount[count]);
    }
}


/// \brief Test the distribution of out edges for a SV locus
///
/// \param set1 The SV locus set to be tested
/// \param locusIndex The locus index to be tested
/// \param expectedEdgeOutCount The expected number of out edges for each node of the tested locus
inline
void
TestSVLocusSetOutEdges(
    const SVLocusSet& set1,
    unsigned locusIndex,
    std::vector<unsigned> expectedEdgeOutCount)
{
    assert(set1.getLocus(locusIndex).size() == expectedEdgeOutCount.size());

    for ( unsigned count = 0; count < expectedEdgeOutCount.size(); count++)
    {
        BOOST_REQUIRE_EQUAL(set1.getLocus(locusIndex).getNode(count).outCount(), expectedEdgeOutCount[count]);
    }
}


/// \brief Check if a query string exists in a file name.
inline
bool
FindStringInFile(std::string fileName, std::string searchString)
{
    std::ifstream iss(fileName);

    bool found = false;
    std::string line;
    while ( std::getline(iss, line))
    {
        if ( line.find(searchString) != std::string::npos)
        {
            found = true;
        }
    }
    return found;
}

