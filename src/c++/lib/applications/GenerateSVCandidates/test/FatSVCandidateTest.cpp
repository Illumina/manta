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

#include "FatSVCandidate.hh"

BOOST_AUTO_TEST_SUITE( FatSVCandidate_test_suite )

// Test the output stream of FatSVCandidate
BOOST_AUTO_TEST_CASE( test_Print_FatSVCandidate)
{
    SVCandidate svCandidate;
    svCandidate.bp1.interval = GenomeInterval(0,10, 100);
    svCandidate.bp2.interval = GenomeInterval(0, 1000, 1100);
    svCandidate.bp1.state  = SVBreakendState::RIGHT_OPEN;
    svCandidate.bp2.state  = SVBreakendState::LEFT_OPEN;
    FatSVCandidate fatSVCandidate(svCandidate, 1u);
    fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3443);
    fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3452);
    fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3440);
    fatSVCandidate.bp1EvidenceIndex[0][0].push_back(3489);

    fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1403);
    fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1428);
    fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1480);
    fatSVCandidate.bp2EvidenceIndex[0][0].push_back(1507);
    std::cout << fatSVCandidate << std::endl;


}


/// Two SVCandidates intersect if their breakend regions overlap in the same direction.
/// In the schematic below, the intersecting candidate pairs are (1,2) and (2,3)
/// Candidate1: >bp2>----------------------------------->>bp1>>
/// Candidate2:    >>>>>>>bp1>>>>>>>--------------->>bp2>>
/// Candidate3:               >>>bp2>>>--------------->>>>>>>bp1>>>>>>
/// Candidate4:               <<<bp2<<<---------------<<<<<<<bp1<<<<<<
BOOST_AUTO_TEST_CASE( test_merge)
{
    SVCandidate svCandidate1;
    svCandidate1.bp1.interval = GenomeInterval(0,10, 100);
    svCandidate1.bp2.interval = GenomeInterval(0, 1000, 1100);
    svCandidate1.bp1.state  = SVBreakendState::RIGHT_OPEN;
    svCandidate1.bp2.state  = SVBreakendState::LEFT_OPEN;

    SVCandidate svCandidate2;
    svCandidate2.bp1.interval = GenomeInterval(0,50, 120);
    svCandidate2.bp2.interval = GenomeInterval(0, 1050, 1150);
    svCandidate2.bp1.state  = SVBreakendState::RIGHT_OPEN;
    svCandidate2.bp2.state  = SVBreakendState::LEFT_OPEN;

    SVCandidate svCandidate3;
    svCandidate3.bp1.interval = GenomeInterval(0, 200, 300);
    svCandidate3.bp2.interval = GenomeInterval(0, 1400, 1500);
    svCandidate3.bp1.state  = SVBreakendState::RIGHT_OPEN;
    svCandidate3.bp2.state  = SVBreakendState::LEFT_OPEN;

    SVCandidate svCandidate4;
    svCandidate3.bp1.interval = GenomeInterval(0, 50, 120);
    svCandidate3.bp2.interval = GenomeInterval(0, 990, 1050);
    svCandidate3.bp1.state  = SVBreakendState::LEFT_OPEN;
    svCandidate3.bp2.state  = SVBreakendState::LEFT_OPEN;

    SVCandidate svCandidate5;
    svCandidate5.bp1.interval = GenomeInterval(0,1050, 1150);
    svCandidate5.bp2.interval = GenomeInterval(0, 50, 120);
    svCandidate5.bp1.state  = SVBreakendState::LEFT_OPEN;
    svCandidate5.bp2.state  = SVBreakendState::RIGHT_OPEN;

    SVCandidate svCandidate6;
    svCandidate6.bp1.interval = GenomeInterval(0,10, 100);
    svCandidate6.bp2.interval = GenomeInterval(0, 1000, 1100);
    svCandidate6.bp1.state  = SVBreakendState::RIGHT_OPEN;
    svCandidate6.bp2.state  = SVBreakendState::LEFT_OPEN;

    FatSVCandidate fatSVCandidate1(svCandidate1, 1u);
    FatSVCandidate fatSVCandidate2(svCandidate2, 1u);
    FatSVCandidate fatSVCandidate3(svCandidate3, 1u);
    FatSVCandidate fatSVCandidate4(svCandidate4, 1u);
    FatSVCandidate fatSVCandidate5(svCandidate5, 1u);
    FatSVCandidate fatSVCandidate6(svCandidate6, 1u);

    // Evidence read indices
    std::vector<int > readIndicesBP1 = {3443,3452,3440,3489,3445,3450,3420,4000};
    std::vector<int > readIndicesBP2 = {1403,1428,1480,1507,1400,142,140,150};

    // adding evidence indices to all fatsvcandidates
    fatSVCandidate1.bp1EvidenceIndex[0][0].push_back(readIndicesBP1[0]);
    fatSVCandidate1.bp1EvidenceIndex[0][0].push_back(readIndicesBP1[1]);
    fatSVCandidate1.bp1EvidenceIndex[0][0].push_back(readIndicesBP1[2]);
    fatSVCandidate1.bp1EvidenceIndex[0][0].push_back(readIndicesBP1[3]);

    fatSVCandidate1.bp2EvidenceIndex[0][0].push_back(readIndicesBP2[0]);
    fatSVCandidate1.bp2EvidenceIndex[0][0].push_back(readIndicesBP2[1]);
    fatSVCandidate1.bp2EvidenceIndex[0][0].push_back(readIndicesBP2[2]);
    fatSVCandidate1.bp2EvidenceIndex[0][0].push_back(readIndicesBP2[3]);

    fatSVCandidate2.bp1EvidenceIndex[0][0].push_back(readIndicesBP1[4]);
    fatSVCandidate2.bp1EvidenceIndex[0][0].push_back(readIndicesBP1[5]);
    fatSVCandidate2.bp1EvidenceIndex[0][0].push_back(readIndicesBP1[6]);
    fatSVCandidate2.bp1EvidenceIndex[0][0].push_back(readIndicesBP1[7]);

    fatSVCandidate2.bp2EvidenceIndex[0][0].push_back(readIndicesBP2[4]);
    fatSVCandidate2.bp2EvidenceIndex[0][0].push_back(readIndicesBP2[5]);
    fatSVCandidate2.bp2EvidenceIndex[0][0].push_back(readIndicesBP2[6]);
    fatSVCandidate2.bp2EvidenceIndex[0][0].push_back(readIndicesBP2[7]);

    fatSVCandidate4.bp1EvidenceIndex[0][0].push_back(3423);
    fatSVCandidate4.bp1EvidenceIndex[0][0].push_back(3412);
    fatSVCandidate4.bp1EvidenceIndex[0][0].push_back(3440);
    fatSVCandidate4.bp1EvidenceIndex[0][0].push_back(3489);

    fatSVCandidate4.bp2EvidenceIndex[0][0].push_back(14023);
    fatSVCandidate4.bp2EvidenceIndex[0][0].push_back(14228);
    fatSVCandidate4.bp2EvidenceIndex[0][0].push_back(14820);
    fatSVCandidate4.bp2EvidenceIndex[0][0].push_back(15037);

    // Although the breakend direction of fatSVCandidate1 and fatSVCandidate3 are same,
    // but their breakend coordinates are not intersecting. They cannot be merged.
    BOOST_REQUIRE(!fatSVCandidate1.merge(fatSVCandidate3));

    // Although the breakpoint-1 of fatSVCandidate1 intersects with the the breakpoint-1 of
    // fatSVCandidate2, but breakpoint-2 direction of fatSVCandidate1 is different from
    // breakpoint-2 direction of fatSVCandidate4. They cannot be merged.
    BOOST_REQUIRE(!fatSVCandidate1.merge(fatSVCandidate4));

    // Checking the size of both breakpoint before merging
    BOOST_REQUIRE_EQUAL(fatSVCandidate1.bp1EvidenceIndex[0][0].size(), 4);
    BOOST_REQUIRE_EQUAL(fatSVCandidate1.bp2EvidenceIndex[0][0].size(), 4);

    // Breakend coordinates are intersecting as well as breakend directions are also
    // same for both fatSVCandidate1 and fatSVCandidate2. They can be merged.
    BOOST_REQUIRE(fatSVCandidate1.merge(fatSVCandidate2));

    // Both fatSVCandidate1 and fatSVCandidate2 are merged and stored in fatSVCandidate1.
    // So the total evidence count of fatSVCandidate1 is sum of evidence count of fatSVCandidate1 and
    // fatSVCandidate2.
    BOOST_REQUIRE_EQUAL(fatSVCandidate1.bp1EvidenceIndex[0][0].size(), readIndicesBP1.size());
    BOOST_REQUIRE_EQUAL(fatSVCandidate1.bp2EvidenceIndex[0][0].size(), readIndicesBP2.size());

    // checking the merged values. It should be sequentially merged.
    for (unsigned i(0); i < readIndicesBP1.size(); i++)
        BOOST_REQUIRE_EQUAL(fatSVCandidate1.bp1EvidenceIndex[0][0][i], readIndicesBP1[i]);
    for (unsigned i(0); i < readIndicesBP2.size(); i++)
        BOOST_REQUIRE_EQUAL(fatSVCandidate1.bp2EvidenceIndex[0][0][i], readIndicesBP2[i]);

    // Breakpoint-1 of fatSVCandidate5 intersects breakpoint-2 of fatSVCandidate6 as well as the direction of
    // breakpoint-1 of fatSVCandidate5 matches with direction of breakpoint-2 of fatSVCandidate6. Similar case
    // for Breakpoint-2 of fatSVCandidate5 and breakpoint-1 of fatSVCandidate6. So they can be merged.
    BOOST_REQUIRE(fatSVCandidate5.merge(fatSVCandidate6));

    // default constructor with everything is empty
    FatSVCandidate fatSVCandidate7;
    for (unsigned evidenceTypeIndex(0); evidenceTypeIndex<SVEvidenceType::SIZE; ++evidenceTypeIndex)
    {
        BOOST_REQUIRE(fatSVCandidate7.bp1EvidenceIndex[evidenceTypeIndex].size() == 0);
        BOOST_REQUIRE(fatSVCandidate7.bp2EvidenceIndex[evidenceTypeIndex].size() == 0);
    }

}

BOOST_AUTO_TEST_SUITE_END()