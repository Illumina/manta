//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2019 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
