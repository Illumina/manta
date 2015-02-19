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

#include "GSCEdgeStats.hh"

#include "boost/archive/xml_iarchive.hpp"
#include "boost/archive/xml_oarchive.hpp"

#include <fstream>
#include <iostream>



template <typename A,typename B>
double
safeFrac(
    const A num,
    const B den)
{
    if (den == 0) return 0;
    return (static_cast<double>(num)/den);
}



void
SimpleHist::
report(std::ostream& os) const
{
    for (unsigned i(0);i<histdata.size();++i)
    {
        os << i;
        if (i+1 == histdata.size()) os << "+";
        os << "\t" << histdata[i] << "\n";
    }
}



void
reportTime(
    const char* label,
    const double ltime,
    const uint64_t edgeCount,
    const uint64_t candCount,
    std::ostream& os)
{
    static const double secPerHour(3600);
    os << label << "Hours_SecsPerEdge_SecsPerCand\t" << ltime/secPerHour
        << "\t" << safeFrac(ltime,edgeCount)
        << "\t" << safeFrac(ltime,candCount)
        << "\n";
}



void
GSCEdgeGroupStats::
report(std::ostream& os) const
{
    const double nocatTime(totalTime-(candTime+assemblyTime+scoringTime));

    os << "InputEdgeCount\t" << totalInputEdgeCount << "\n";
    os << "InputEdgeCandidatesPerEdge:\n";
    candidatesPerEdge.report(os);
    os << "CandidateCount\t" << totalCandidateCount << "\n";
    os << "ComplexCandidateCount\t" << totalComplexCandidate << "\n";
    os << "SpanningComplexCandidateFiltered\t" << totalSpanningCandidateFilter << "\n";
    os << "JunctionCount\t" << totalJunctionCount << "\n";
    os << "ComplexJunctionCount\t" << totalComplexJunctionCount << "\n";
    os << "BreaksPerJunction:\n";
    breaksPerJunction.report(os);
    os << "TotalAssemblyCandidates\t" << totalAssemblyCandidates << "\n";
    os << "TotalSpanningAssemblyCandidates\t" << totalSpanningAssemblyCandidates << "\n";
    os << "AssemblyCandidatesPerJunction:\n";
    assemblyCandidatesPerJunction.report(os);
    reportTime("total",totalTime,totalInputEdgeCount,totalCandidateCount, os);
    reportTime("candi",candTime,totalInputEdgeCount,totalCandidateCount, os);
    reportTime("assem",assemblyTime,totalInputEdgeCount,totalCandidateCount, os);
    reportTime("score",scoringTime,totalInputEdgeCount,totalCandidateCount, os);
    reportTime("nocat",nocatTime,totalInputEdgeCount,totalCandidateCount, os);
}



void
GSCEdgeStatsData::
report(std::ostream& os) const
{
    GSCEdgeGroupStats all(remoteEdges);
    all.merge(selfEdges);
    os << "SVGenTotalHours\t" << lifeTime/3600 << "\n";
    os << "NonEdgeHours\t" << (lifeTime-all.totalTime)/3600 << "\n";
    os << "\n[AllEdges]\n";
    all.report(os);
    os << "\n[RemoteEdges]\n";
    remoteEdges.report(os);
    os << "\n[SelfEdges]\n";
    selfEdges.report(os);
}



void
GSCEdgeStats::
load(const char* filename)
{
    assert(nullptr != filename);
    std::ifstream ifs(filename);
    boost::archive::xml_iarchive ia(ifs);
    ia >> BOOST_SERIALIZATION_NVP(edgeData);
}



void
GSCEdgeStats::
save(std::ostream& os) const
{
    boost::archive::xml_oarchive oa(os);
    oa << BOOST_SERIALIZATION_NVP(edgeData);
}



void
GSCEdgeStats::
save(const char* filename) const
{
    assert(nullptr != filename);
    std::ofstream ofs(filename);
    save(ofs);
}



void
GSCEdgeStats::
report(const char* filename) const
{
    assert(nullptr != filename);
    std::ofstream ofs(filename);
    ofs << "EdgeStatsReport\n";
    edgeData.report(ofs);
}
