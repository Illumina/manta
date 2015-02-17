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
GSCEdgeGroupStats::
report(std::ostream& os) const
{
    const double secPerHour(3600);
    const double totalHours(totalTime/secPerHour);
    const double candTime(totalTime-assemblyTime-scoringTime);
    const double candHours(candTime/secPerHour);
    const double assemblyHours(assemblyTime/secPerHour);
    const double scoringHours(scoringTime/secPerHour);

    os << "totalEdgeCount:\t" << totalEdgeCount << "\n";
    os << "totalHours_SecsPerEdge:\t" << totalHours << "\t" << safeFrac(totalTime,totalEdgeCount) << "\n";
    os << "candHours_SecsPerEdge:\t" << candHours << "\t" << safeFrac(candTime,totalEdgeCount) << "\n";
    os << "assmHours_SecsPerEdge:\t" << assemblyHours << "\t" << safeFrac(assemblyTime,totalEdgeCount) << "\n";
    os << "scoreHours_SecsPerEdge:\t" << scoringHours << "\t" << safeFrac(scoringTime,totalEdgeCount) << "\n";
}



void
GSCEdgeStatsData::
report(std::ostream& os) const
{
    GSCEdgeGroupStats all(remoteEdges);
    all.merge(selfEdges);
    os << "AllEdges:\n";
    all.report(os);
    os << "RemoteEdges:\n";
    remoteEdges.report(os);
    os << "SelfEdges:\n";
    selfEdges.report(os);
}



void
GSCEdgeStats::
load(const char* filename)
{
    assert(nullptr != filename);
    std::ifstream ifs(filename);
    boost::archive::xml_iarchive ia(ifs);
    ia >> BOOST_SERIALIZATION_NVP(data);
}



void
GSCEdgeStats::
save(std::ostream& os) const
{
    boost::archive::xml_oarchive oa(os);
    oa << BOOST_SERIALIZATION_NVP(data);
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
    ofs << "EdgeStatsReport:";
    data.report(ofs);
}
