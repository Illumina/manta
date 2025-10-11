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

/// \file
/// \brief Take care of some (mostly C99) functions not available in VS C++
///

#pragma once

#include <string>

#ifdef _MSC_VER
#define snprintf _snprintf
#define strdup _strdup
#endif

#if ((defined(_MSC_VER)) && (_MSC_VER <= 1800))
#undef noexcept
#define noexcept
#endif

double compat_round(const double x);

const char* compat_basename(const char* s);
