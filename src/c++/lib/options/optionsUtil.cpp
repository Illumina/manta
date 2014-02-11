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

#include "optionsUtil.hh"

#include "boost/filesystem.hpp"

#include <sstream>



bool
checkStandardizeInputFile(
    std::string& filename,
    const char* fileLabel,
    std::string& errorMsg)
{
    errorMsg.clear();

    if (filename.empty())
    {
        std::ostringstream oss;
        oss << "Must specify " << fileLabel << " file";
        errorMsg = oss.str();
    }
    else if (! boost::filesystem::exists(filename))
    {
        std::ostringstream oss;
        oss << "Can't find " << fileLabel << " file '" << filename << "'";
        errorMsg = oss.str();
    }
    else
    {
        filename = boost::filesystem::canonical(filename).string();
    }

    return (! errorMsg.empty());
}
