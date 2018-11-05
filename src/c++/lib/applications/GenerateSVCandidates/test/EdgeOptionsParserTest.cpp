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

#include "boost/test/unit_test.hpp"

#include "EdgeOptionsParser.hh"

using namespace std;

BOOST_AUTO_TEST_SUITE( EdgeOptionsParser_test_suite )

// Test the description of the parser for edge options
BOOST_AUTO_TEST_CASE( test_OptionsDescription )
{
    EdgeOptions opt;
    boost::program_options::options_description optdesc = getOptionsDescription(opt);

    // check description of bin count option
    BOOST_REQUIRE_EQUAL(optdesc.options().at(0).get()->format_name(), "--bin-count");
    BOOST_REQUIRE_EQUAL(optdesc.options().at(0).get()->format_parameter(), "arg (=1)");
    BOOST_REQUIRE_EQUAL(optdesc.options().at(0).get()->description(), "Specify how many bins "
                                                                      "the SV candidate problem "
                                                                      "should be divided into, "
                                                                      "where bin-index can be "
                                                                      "used to specify which bin "
                                                                      "to solve");

    // check description bin index option
    BOOST_REQUIRE_EQUAL(optdesc.options().at(1).get()->format_name(), "--bin-index");
    BOOST_REQUIRE_EQUAL(optdesc.options().at(1).get()->format_parameter(), "arg (=0)");
    BOOST_REQUIRE_EQUAL(optdesc.options().at(1).get()->description(), "Specify which bin to solve "
                                                                      "when the SV candidate "
                                                                      "problem is subdivided into "
                                                                      "bins. Value must bin in "
                                                                      "[0,bin-count)");
    // check description of max edge count option
    BOOST_REQUIRE_EQUAL(optdesc.options().at(2).get()->format_name(), "--max-edge-count");
    BOOST_REQUIRE_EQUAL(optdesc.options().at(2).get()->format_parameter(), "arg (=10)");
    BOOST_REQUIRE_EQUAL(optdesc.options().at(2).get()->description(), "Specify the maximum number of edge count. "
                                                                      "If both nodes of an edge have an edge count "
                                                                      "higher than this, the edge is skipped for "
                                                                      "evaluation.");

    // check description of locus index option
    BOOST_REQUIRE_EQUAL(optdesc.options().at(3).get()->format_name(), "--locus-index");
    BOOST_REQUIRE_EQUAL(optdesc.options().at(3).get()->format_parameter(), "arg");
    BOOST_REQUIRE_EQUAL(optdesc.options().at(3).get()->description(), "Instead of solving for all SV candidates in a "
                                                                      "bin, solve for candidates of a particular locus "
                                                                      "or edge. If this argument is specified then bin-index "
                                                                      "is ignored. Argument can be one of { locusIndex , "
                                                                      "locusIndex:nodeIndex , locusIndex:nodeIndex:nodeIndex }, "
                                                                      "which will run an entire locus, all edges connected to "
                                                                      "one node in a locus or a single edge, respectively.");
}

// Test the parsing of edge options.
BOOST_AUTO_TEST_CASE( test_ParseEdgeOptions )
{
    EdgeOptions edgeoption;
    namespace po = boost::program_options;
    const char locusIndexKey[] = "locus-index";
    std::string errorMsg;
    boost::program_options::variables_map vm;
    string locusIndex = "1:2:3:4";
    vm.insert(std::make_pair(locusIndexKey, po::variable_value(locusIndex, true)));

    // locus index format is not correct. It should not have more than 3 colon separated segments.
    // locus index format is either of the following two things:
    // 1. locusIndex:nodeIndex
    // 2. locusIndex:nodeIndex:nodeIndex
    parseOptions(vm, edgeoption, errorMsg);
    BOOST_REQUIRE_EQUAL(errorMsg, "locus-index argument can have no more than 3 colon separated segments");

    locusIndex = "1:2:3"; // correct locus-index format
    vm.clear();
    vm.insert(std::make_pair(locusIndexKey, po::variable_value(locusIndex, true)));
    edgeoption.binCount = 0;

    // Divide all edges in the graph into binCount bins of approx equal complexity.
    // So that at least 1 bin count is required.
    parseOptions(vm, edgeoption, errorMsg);
    BOOST_REQUIRE_EQUAL(errorMsg, "bin-count must be 1 or greater");

    edgeoption.binCount = 1;
    edgeoption.binIndex = 1;
    // Out of binCount bins, iterate through the edges in this bin only.
    // So the bin index should be less than total bin count.
    parseOptions(vm, edgeoption, errorMsg);
    BOOST_REQUIRE_EQUAL(errorMsg, "bin-index must be in range [0,bin-count)");

    // Successful case when all parameters are correct
    edgeoption.binCount = 1;
    edgeoption.binIndex = 0;
    parseOptions(vm, edgeoption, errorMsg);
    BOOST_REQUIRE_EQUAL(errorMsg, "");
}

BOOST_AUTO_TEST_SUITE_END()