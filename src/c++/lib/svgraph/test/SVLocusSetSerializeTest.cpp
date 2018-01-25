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

#include "boost/archive/tmpdir.hpp"
#include "boost/test/unit_test.hpp"

#include "svgraph/SVLocusSet.hh"

#include "SVLocusTestUtil.hh"

using namespace boost::archive;


BOOST_AUTO_TEST_SUITE( test_SVLocusSetSerialize )


BOOST_AUTO_TEST_CASE( test_SVLocusSetSerialze )
{
    // construct a simple two-node locus
    SVLocus locus1;
    locusAddPair(locus1,1,10,20,2,30,40);

    SVLocus locus2;
    locusAddPair(locus2,3,10,20,4,30,40);

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


BOOST_AUTO_TEST_CASE( test_SVLocusSetSerialze2 )
{
    SVLocusSet set1;
    {
        SVLocus locus1;
        locusAddPair(locus1,1,10,20,2,30,40);

        SVLocus locus2;
        locusAddPair(locus2,3,10,20,4,30,40);

        set1.merge(locus1);
        set1.merge(locus2);
    }

    SVLocusSet set2;
    {
        SVLocus locus1;
        locusAddPair(locus1,1,15,25,4,30,40);

        SVLocus locus2;
        locusAddPair(locus2,3,30,40,2,30,40);

        set2.merge(locus1);
        set2.merge(locus2);
    }

    SVLocusSet set1_copy;
    {
        std::string filename(tmpdir());
        filename += "/testfile.bin";

        // serialize
        set1.save(filename.c_str());

        // deserialize
        set1_copy.load(filename.c_str());
    }

    SVLocusSet set2_copy;
    {
        std::string filename(tmpdir());
        filename += "/testfile.bin";

        // serialize
        set2.save(filename.c_str());

        // deserialize
        set2_copy.load(filename.c_str());
    }

    set1.merge(set2);
    set1_copy.merge(set2_copy);

    BOOST_REQUIRE_EQUAL(set1.size(),set1_copy.size());

    typedef SVLocusSet::const_iterator citer;

    citer i(set1.begin());
    citer i_copy(set1_copy.begin());

    const SVLocus& set1_locus1(*i);
    const SVLocus& set1_copy_locus1(*i_copy);
    BOOST_REQUIRE_EQUAL(set1_locus1.size(),set1_copy_locus1.size());
}



BOOST_AUTO_TEST_SUITE_END()

