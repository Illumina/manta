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

#pragma once

#include "options/AlignmentFileOptions.hpp"

#include "boost/program_options.hpp"

#include <string>

boost::program_options::options_description getOptionsDescription(AlignmentFileOptions& opt);

/// \brief Parse AlignmentFileOptions
///
/// \param[out] errorMsg If an error occurs this is set to an end-user targeted error message. Any string
/// content on input is cleared
///
/// \return True if an error occurs while parsing options
bool parseOptions(
    const boost::program_options::variables_map& vm, AlignmentFileOptions& opt, std::string& errorMsg);
