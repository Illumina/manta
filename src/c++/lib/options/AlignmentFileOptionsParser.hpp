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
