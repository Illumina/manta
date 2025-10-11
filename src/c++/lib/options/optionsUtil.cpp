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

#include "optionsUtil.hpp"

#include "boost/filesystem.hpp"

#include <sstream>

bool checkAndStandardizeRequiredInputFilePath(
    std::string& filename, const char* fileLabel, std::string& errorMsg)
{
  errorMsg.clear();

  if (filename.empty()) {
    std::ostringstream oss;
    oss << "Must specify " << fileLabel << " file";
    errorMsg = oss.str();
  } else if (!boost::filesystem::exists(filename)) {
    std::ostringstream oss;
    oss << "Can't find " << fileLabel << " file '" << filename << "'";
    errorMsg = oss.str();
  } else {
    filename = boost::filesystem::absolute(filename).string();
  }

  return (!errorMsg.empty());
}
