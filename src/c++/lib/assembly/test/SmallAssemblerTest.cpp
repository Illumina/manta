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

#include "boost/test/unit_test.hpp"

#include "SmallAssembler.hh"


BOOST_AUTO_TEST_SUITE( test_SmallAssembler )


BOOST_AUTO_TEST_CASE( test_SmallAssembler1 )
{
    SmallAssemblerOptions assembleOpt;

    assembleOpt.wordLength = 5;
    assembleOpt.maxWordLength = 5;

    // ********************  still setting this up....

    //runSmallAssembler(assembleOpt,,);
    
    BOOST_REQUIRE(true);
}


BOOST_AUTO_TEST_SUITE_END()

