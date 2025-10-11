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

/**
 ** \file
 ** \brief Implementation of the common exception mechanism.
 **
 ** \author Come Raczy
 **/

#include "Exceptions.hpp"

#include <boost/date_time.hpp>
#include <cerrno>
#include <cstring>

namespace illumina {
namespace common {

std::string ExceptionData::getContext() const
{
  const std::string now = boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time());
  std::string       errorInfo;
  if (_errorNumber != 0) errorInfo = " '" + std::string(strerror(_errorNumber)) + "'";
  return now + errorInfo + " " + boost::diagnostic_information(*this);
}

}  // namespace common
}  // namespace illumina
