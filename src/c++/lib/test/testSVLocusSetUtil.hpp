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
/// \brief Unit testing utility functions for SVLocusSet elements
///

#pragma once

#include "svgraph/SVLocusSet.hpp"

#include <memory>

/// Serialize, then deserialize object and pass back a pointer to the resulting object copy
std::unique_ptr<SVLocusSet> getSerializedSVLocusSetCopy(const SVLocusSet& set);
