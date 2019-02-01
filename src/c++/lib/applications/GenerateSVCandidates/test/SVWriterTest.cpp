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

#include "SVWriter.cpp" // cpp required to pull in checkJunctionsToFilter


BOOST_AUTO_TEST_SUITE( SVWriter_test_suite )


// A candidate sv needs to be checked after assembly whether it is good candidate sv or not.
// A candidate sv is not good if either of the following three cases is  satisfied:
// 1. Spanning evidence count of a candidate is less than min candidate spanning count threshold (default 3).
// 2. Non-spanning low-res candidate went into assembly but did not produce a successful contig alignment
//    that means SV is imprecise.
// 3. Variant size is smaller than minCandidateVariantSize (default is 10)
// Successful cases are designed in test_tumorOnly, test_RNA, test_Diploid and test_Somatic.
BOOST_AUTO_TEST_CASE( test_JunctionFilter )
{
    GSCOptions options;
    // Creating SV candidate
    SVCandidate candidate1;
    candidate1.setPrecise();
    candidate1.bp1.state = SVBreakendState::RIGHT_OPEN;
    candidate1.bp1.interval = GenomeInterval(0 , 40, 50);
    candidate1.bp1.lowresEvidence.add(0, 1);
    candidate1.bp2.state = SVBreakendState::LEFT_OPEN;
    candidate1.bp2.interval = GenomeInterval(1 , 65, 75);
    // Adding a pair evidence to this candidate
    candidate1.bp2.lowresEvidence.add(0, 1);
    SVMultiJunctionCandidate junctionCandidate;
    junctionCandidate.junction = {candidate1};
    std::vector<bool> isInputJunctionFiltered1;
    isInputJunctionFiltered1.resize(1);
    SVCandidateAssemblyData candidateAssemblyData;
    candidateAssemblyData.isCandidateSpanning = true;
    std::vector <SVCandidateAssemblyData> assemblyData;
    assemblyData.push_back(candidateAssemblyData);
    // Case-1 is designed where number of spanning evidence observations required is 3.
    // But here spanning observation count is 1.
    checkJunctionsToFilter(junctionCandidate, assemblyData, isInputJunctionFiltered1, options);
    BOOST_REQUIRE_EQUAL(isInputJunctionFiltered1[0], true);

    // Case-2 is designed
    // In this case a non-spanning low-res candidate went into assembly but
    // did not produce a successful contig alignment.
    SVCandidate candidate2;
    candidate2.bp1.state = SVBreakendState::COMPLEX;
    candidate2.bp1.interval = GenomeInterval(0 , 40, 50);
    candidate2.bp1.lowresEvidence.add(0, 1);
    candidate2.bp2.state = SVBreakendState::UNKNOWN;
    candidate2.bp2.interval = GenomeInterval(0 , 65, 75);
    candidate2.bp2.lowresEvidence.add(0, 1);
    junctionCandidate.junction.clear();
    junctionCandidate.junction.push_back(candidate2);
    candidateAssemblyData.isCandidateSpanning = false;
    assemblyData.clear();
    assemblyData.push_back(candidateAssemblyData);
    std::vector<bool> isInputJunctionFiltered2;
    isInputJunctionFiltered2.resize(1);
    checkJunctionsToFilter(junctionCandidate, assemblyData, isInputJunctionFiltered2, options);
    BOOST_REQUIRE_EQUAL(isInputJunctionFiltered2[0], true);

    // Case-3 is designed.
    // variant size is less than min variant size(40)
    options.minCandidateSpanningCount = 1;
    SVCandidate candidate3;
    candidate3.bp1.state = SVBreakendState::RIGHT_OPEN;
    candidate3.bp1.interval = GenomeInterval(0 , 40, 50);
    candidate3.bp1.lowresEvidence.add(0, 1);
    candidate3.bp2.state = SVBreakendState::LEFT_OPEN;
    candidate3.bp2.interval = GenomeInterval(0 , 65, 75);
    candidate3.bp2.lowresEvidence.add(0, 1);
    junctionCandidate.junction.clear();
    junctionCandidate.junction.push_back(candidate2);
    candidateAssemblyData.isCandidateSpanning = true;
    assemblyData.clear();
    assemblyData.push_back(candidateAssemblyData);
    options.scanOpt.minCandidateVariantSize = 40;
    std::vector<bool> isInputJunctionFiltered3;
    isInputJunctionFiltered3.resize(1);
    // Variant size = 65 - 40 -1 = 24 which is less than 40.
    checkJunctionsToFilter(junctionCandidate, assemblyData, isInputJunctionFiltered3, options);
    BOOST_REQUIRE_EQUAL(isInputJunctionFiltered3[0], true);
}

BOOST_AUTO_TEST_SUITE_END()
