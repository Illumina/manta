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

///
/// \author Chris Saunders
///

#include "boost/archive/tmpdir.hpp"
#include "boost/archive/text_iarchive.hpp"
#include "boost/archive/text_oarchive.hpp"
#include "boost/archive/binary_iarchive.hpp"
#include "boost/archive/binary_oarchive.hpp"
#include "boost/foreach.hpp"
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

    BOOST_REQUIRE_EQUAL(node1.count,node_copy1.count);
    BOOST_REQUIRE_EQUAL(node1.interval,node_copy1.interval);
    BOOST_REQUIRE_EQUAL(node1.size(),node_copy1.size());

    SVLocusNode::const_iterator ibegin(node1.begin());
    SVLocusNode::const_iterator copy_ibegin(node_copy1.begin());

    BOOST_REQUIRE_EQUAL(ibegin->second.count, copy_ibegin->second.count);

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
    BOOST_FOREACH(const SVLocusNode& node, clocus1_copy)
    {
        if (node.interval == (clocus1.begin())->interval) isMatchFound=true;
    }
    BOOST_REQUIRE(isMatchFound);
}


BOOST_AUTO_TEST_SUITE_END()

