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

#include <cstdlib>
#include <iostream>

#include "blt_util/log.hpp"
#include "common/OutStream.hpp"
#include "common/Program.hpp"

#include "TestAssemblerOptions.hpp"
#include "assembly/IterativeAssembler.hpp"
#include "assembly/SmallAssembler.hpp"
#include "extractAssemblyReads.hpp"

/// test front-end to run the manta assembler from command-line
///
struct TestAssembler : public illumina::Program {
  const char* name() const { return "TestAssembler"; }

  void runInternal(int argc, char* argv[]) const;
};
