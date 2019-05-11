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
