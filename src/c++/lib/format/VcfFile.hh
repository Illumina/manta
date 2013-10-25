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
// <https://github.com/sequencing/licenses/>
//

///
/// \author Richard Shaw
///

/*****************************************************************************/

#pragma once

#include <string>
#include <map>
#include <fstream>

#include "VcfHeader.hh"
#include "VcfLine.hh"
#include "Variant.hh"

/*****************************************************************************/

class VcfFile
{
public:
    VcfFile(const std::string pathStr, std::map<std::string, int32_t>);
    bool getVariant(Variant& variant, bool& wasLast);
    bool getVariantVec(VariantVec& variantVec);

private:
    bool loadHeader();

    typedef std::map<std::string, Variant::Type> VcfTypeMap;
    typedef VcfTypeMap::const_iterator VcfTypeMapCIter;
    VcfTypeMap myVcfTypeMap;

    typedef std::map<std::string, bool> MateSeenMap;
    typedef MateSeenMap::iterator MateSeenMapIter;
    MateSeenMap myMateSeenMap;

    std::string myPathStr;
    std::map<std::string, int32_t> myChromNameTidMap;
    std::ifstream myStrm;
    bool myVcfHeaderLoadedFlag;
    VcfHeader myVcfHeader;

    size_t mySvTypeFieldId;
    size_t mySvLenFieldId;
    size_t myMateIdFieldId;
};

/*****************************************************************************/

