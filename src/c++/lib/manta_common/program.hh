// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/downloads/sequencing/licenses/>.
//

/// \author Chris Saunders
///

#pragma once

#include "manta_common/version.hh"

namespace manta {

/// base-class for all command-line programs
///
/// this is used to standardize bottom-level exception handling
struct program {

    virtual
    ~program() {}

    int
    run(int argc, char* argv[]) const;

    virtual
    const char*
    name() const = 0;

    const char*
    version() const {
        return manta::get_version_full();
    }

protected:
    virtual
    void
    run_internal(int argc, char* argv[]) const = 0;
};

}
