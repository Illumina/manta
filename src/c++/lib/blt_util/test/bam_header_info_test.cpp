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

#include "boost/test/unit_test.hpp"

#include "bam_header_info.hh"

#include "boost/archive/binary_iarchive.hpp"
#include "boost/archive/binary_oarchive.hpp"
#include "boost/archive/tmpdir.hpp"

#include <fstream>

BOOST_AUTO_TEST_SUITE( test_bam_header_info )


BOOST_AUTO_TEST_CASE( test_bam_header_info_serialize)
{
    using namespace boost::archive;

    bam_header_info header;
    header.chrom_data.push_back(bam_header_info::chrom_info("chr1",1000));
    header.chrom_data.push_back(bam_header_info::chrom_info("chr2",2000));

    std::string filename(boost::archive::tmpdir());
    filename += "/testfile.bin";

    // serialize
    {
        std::ofstream ofs(filename.c_str(), std::ios::binary);
        binary_oarchive oa(ofs);
        oa << header;
    }

    bam_header_info header2;

    // deserialize
    {
        std::ifstream ifs(filename.c_str(), std::ios::binary);
        binary_iarchive ia(ifs);
        ia >> header2;
    }

    BOOST_REQUIRE_EQUAL(header,header2);
}


BOOST_AUTO_TEST_SUITE_END()
