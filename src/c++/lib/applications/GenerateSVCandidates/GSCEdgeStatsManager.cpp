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

#include "GSCEdgeStatsManager.hh"
#include "common/Exceptions.hh"

#include "boost/archive/xml_iarchive.hpp"
#include "boost/archive/xml_oarchive.hpp"

#include <fstream>
#include <iostream>
#include <sstream>



GSCEdgeStatsManager::
GSCEdgeStatsManager(
    const std::string& outputFile)
    : _osPtr(nullptr)
{
    if (outputFile.empty()) return;
    _osPtr = new std::ofstream(outputFile.c_str());
    if (! *_osPtr)
    {
        std::ostringstream oss;
        oss << "ERROR: Can't open output file: " << outputFile << '\n';
        BOOST_THROW_EXCEPTION(illumina::common::LogicException(oss.str()));
    }
}



GSCEdgeStatsManager::
~GSCEdgeStatsManager()
{
    if (_osPtr != nullptr)
    {
        save(*_osPtr);
    }
}



void
GSCEdgeStatsManager::
load(std::istream& is)
{
    boost::archive::xml_iarchive ia(is);
    ia >> BOOST_SERIALIZATION_NVP(edgeStats);
}



void
GSCEdgeStatsManager::
save(std::ostream& os) const
{
    boost::archive::xml_oarchive oa(os);
    oa << BOOST_SERIALIZATION_NVP(edgeStats);
}
