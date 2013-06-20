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
#include "boost/test/unit_test.hpp"

#include "manta/SVLocusSet.hh"

using namespace boost::archive;


BOOST_AUTO_TEST_SUITE( test_SVLocusSetSerialize )


BOOST_AUTO_TEST_CASE( test_SVLocusSetSerialze )
{
    // construct a simple two-node locus
    SVLocus locus1;
    SVLocusNode* nodePtr1 = locus1.addNode(1,10,20);
    locus1.addNode(2,30,40,nodePtr1);

    SVLocus locus2;
    SVLocusNode* nodePtr2 = locus2.addNode(3,10,20);
    locus2.addNode(4,30,40,nodePtr2);

    SVLocusSet set1;
    set1.merge(locus1);
    set1.merge(locus2);

    std::string filename(tmpdir());
    filename += "/testfile.bin";

    // serialize
    set1.save(filename.c_str());

    SVLocusSet set1_copy;

    // deserialize
    set1_copy.load(filename.c_str());

    BOOST_REQUIRE_EQUAL(set1.size(),set1_copy.size());

    typedef SVLocusSet::const_iterator citer;

    citer i(set1.begin());
    citer i_copy(set1_copy.begin());

    const SVLocus& set1_locus1(*i);
    const SVLocus& set1_copy_locus1(*i_copy);
    BOOST_REQUIRE_EQUAL(set1_locus1.size(),set1_copy_locus1.size());
}


BOOST_AUTO_TEST_SUITE_END()

