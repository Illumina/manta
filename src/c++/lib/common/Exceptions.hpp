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
 ** \brief Declaration of the common exception mechanism.
 **
 ** All exceptions must carry the same data (independently of the
 ** exception type) to homogenize the reporting and processing of
 ** errors.
 **
 ** \author Come Raczy
 **/

#pragma once

#include "blt_util/thirdparty_push.h"

#include "boost/cerrno.hpp"
#include "boost/exception/all.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/throw_exception.hpp"

#include "blt_util/thirdparty_pop.h"

#include <ios>
#include <stdexcept>
#include <string>

namespace illumina {
namespace common {

/// \brief Virtual base class to all the exception classes
///
/// Use BOOST_THROW_EXCEPTION to get the context info (file, function, line)
/// at the throw site.
///
class ExceptionData : public boost::exception {
public:
  ExceptionData(const std::string& message, const int errorNumber = 0)
    : boost::exception(), _message(message), _errorNumber(errorNumber)
  {
  }

  ExceptionData(const ExceptionData&) = default;
  ExceptionData& operator=(const ExceptionData&) = delete;

  std::string getContext() const;

private:
  const std::string _message;
  const int         _errorNumber;
};

/// A general purpose exception type
///
/// Use BOOST_THROW_EXCEPTION to get the context info (file, function, line)
/// at the throw site as follows:
///
///     BOOST_THROW_EXCEPTION(GeneralException("Error message"));
///
class GeneralException : public std::logic_error, public ExceptionData {
public:
  explicit GeneralException(const std::string& message, const int errorNumber = 0)
    : std::logic_error(message), ExceptionData(message, errorNumber)
  {
  }
};

}  // namespace common
}  // namespace illumina
