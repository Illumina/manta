// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2015 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

/// \brief provide access to cmake project version numbers

#pragma once

#include "common/config.h"

namespace manta
{

inline
const char*
getVersion()
{
    return WORKFLOW_VERSION;
}

inline
const char*
getBuildTime()
{
    return BUILD_TIME;
}

inline
const char*
cxxCompilerName()
{
    return CXX_COMPILER_NAME;
}

inline
const char*
compilerVersion()
{
    return COMPILER_VERSION;
}

}
