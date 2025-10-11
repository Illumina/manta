//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2025 Illumina, Inc.
//
// This program is licensed under the terms of the Polyform strict license
//
// ***As far as the law allows, the software comes as is, without
// any warranty or condition, and the licensor will not be liable
// to you for any damages arising out of these terms or the use
// or nature of the software, under any kind of legal claim.***
//
// You should have received a copy of the PolyForm Strict License 1.0.0
// along with this program.  If not, see <https://polyformproject.org/licenses/strict/1.0.0>.
//
//

/// \brief provide access to cmake project version numbers

#pragma once

#include "common/config.h"

namespace illumina {

inline const char* getVersion()
{
  return WORKFLOW_VERSION;
}

inline const char* getBuildTime()
{
  return BUILD_TIME;
}

inline const char* cxxCompilerName()
{
  return CXX_COMPILER_NAME;
}

inline const char* compilerVersion()
{
  return COMPILER_VERSION;
}

}  // namespace illumina
