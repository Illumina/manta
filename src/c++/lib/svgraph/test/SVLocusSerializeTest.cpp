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

#include "boost/archive/tmpdir.hpp"
#include "boost/archive/text_iarchive.hpp"
#include "boost/archive/text_oarchive.hpp"
#include "boost/archive/binary_iarchive.hpp"
#include "boost/archive/binary_oarchive.hpp"
#include "boost/test/unit_test.hpp"

#include "svgraph/SVLocus.hh"

#include "SVLocusTestUtil.hh"

#include <fstream>

using namespace boost::archive;


BOOST_AUTO_TEST_SUITE( test_SVLocusSerialize )

// test serialization with a very simple class first:
//
template <typename InputArchiver, typename OutputArchiver>
void
GenomeIntervalSerializeTest(const char* extension)
{
    // construct a simple two-node locus
    GenomeInterval gi(1,10,20);

    std::string filename(boost::archive::tmpdir());
    filename += "/testfile";
    filename += extension;

    // serialize
    {
        std::ofstream ofs(filename.c_str(), std::ios::binary);
        OutputArchiver oa(ofs);
        oa << gi;
    }

    GenomeInterval gi_copy;

    // deserialize
    {
        std::ifstream ifs(filename.c_str(), std::ios::binary);
        InputArchiver ia(ifs);
        ia >> gi_copy;
    }

    BOOST_REQUIRE_EQUAL(gi,gi_copy);
}


BOOST_AUTO_TEST_CASE( test_GenomeIntervalSerializeText )
{
    GenomeIntervalSerializeTest<text_iarchive,text_oarchive>(".txt");
}


BOOST_AUTO_TEST_CASE( test_GenomeIntervalSerializeBinary )
{
    GenomeIntervalSerializeTest<binary_iarchive,binary_oarchive>(".bin");
}


BOOST_AUTO_TEST_CASE( test_SVLocusNodeSerialze )
{

    // construct a simple two-node locus
    SVLocus locus1;
    locusAddPair(locus1,1,10,20,1,30,40);

    std::string filename(boost::archive::tmpdir());
    filename += "/testfile.bin";

    const SVLocusNode& node1(static_cast<const SVLocus&>(locus1).getNode(0));

    // serialize
    {
        std::ofstream ofs(filename.c_str(), std::ios::binary);
        binary_oarchive oa(ofs);
        oa << node1;
    }

    SVLocusNode node_copy1;

    // deserialize
    {
        std::ifstream ifs(filename.c_str(), std::ios::binary);
        binary_iarchive ia(ifs);
        ia >> node_copy1;
    }

    BOOST_REQUIRE_EQUAL(node1.outCount(), node_copy1.outCount());
    BOOST_REQUIRE_EQUAL(node1.getInterval(), node_copy1.getInterval());
    BOOST_REQUIRE_EQUAL(node1.size() ,node_copy1.size());

    const SVLocusEdgeManager node1Manager(node1.getEdgeManager());
    const SVLocusEdgeManager node1CopyManager(node_copy1.getEdgeManager());

    SVLocusEdgesType::const_iterator ibegin(node1Manager.getMap().begin());
    SVLocusNode::const_iterator copy_ibegin(node1CopyManager.getMap().begin());

    BOOST_REQUIRE_EQUAL(ibegin->second.getCount(), copy_ibegin->second.getCount());

    BOOST_REQUIRE_EQUAL(ibegin->first, copy_ibegin->first);
}



BOOST_AUTO_TEST_CASE( test_SVLocusSerialze )
{

    // construct a simple two-node locus
    SVLocus locus1;
    locusAddPair(locus1,1,10,20,1,30,40);

    std::string filename(boost::archive::tmpdir());
    filename += "/testfile.bin";

    // serialize
    {
        std::ofstream ofs(filename.c_str(), std::ios::binary);
        binary_oarchive oa(ofs);
        oa << locus1;
    }

    SVLocus locus1_copy;

    // deserialize
    {
        std::ifstream ifs(filename.c_str(), std::ios::binary);
        binary_iarchive ia(ifs);
        ia >> locus1_copy;
    }

    BOOST_REQUIRE_EQUAL(locus1.size(),locus1_copy.size());
    const SVLocus& clocus1(locus1);
    const SVLocus& clocus1_copy(locus1_copy);

    bool isMatchFound(false);
    for (const SVLocusNode& node : clocus1_copy)
    {
        if (node.getInterval() == (clocus1.begin())->getInterval()) isMatchFound=true;
    }
    BOOST_REQUIRE(isMatchFound);
}


BOOST_AUTO_TEST_SUITE_END()

