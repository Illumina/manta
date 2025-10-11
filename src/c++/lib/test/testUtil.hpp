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
/// \brief Utility functions for unit testing
/// \author Trevor Ramsay
///

#pragma once

#include <string>

/// \brief Return the name of a new temp file which is unused at the time the function is called
std::string getNewTempFile();

/// \brief Return the value for the given key from a key/value tsv file
///
/// \param tsvFile Input file to search
/// \param key Key value to search for in first column of \p tsvFile
/// \return Second-column value from first line in \p tsvFile matching \p key, or empty string if key not
/// found.
std::string getValueFromTSVKeyValFile(const std::string& tsvFile, const std::string& key);
