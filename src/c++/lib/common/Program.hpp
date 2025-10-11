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
/// \author Chris Saunders
///

#pragma once

#include <iosfwd>

namespace illumina {

/// base-class for all command-line programs
///
/// this is used to standardize bottom-level exception handling
struct Program {
  virtual ~Program() = default;

  int run(int argc, char* argv[]) const;

  virtual const char* name() const = 0;

  const char* version() const;

  const char* compiler() const;

  const char* buildTime() const;

protected:
  virtual void runInternal(int argc, char* argv[]) const = 0;

private:
  void post_catch(int argc, char* argv[], std::ostream& os) const;
};

}  // namespace illumina
