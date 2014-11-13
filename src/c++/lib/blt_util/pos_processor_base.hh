// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#pragma once

#include "blt_util/blt_types.hh"


/// \brief base for objects designed to perform work in a single pass over a position range
///
/// Work progress is communicated via the process_pos() method. This base class is designed to
/// link the worker object with the stage_manager object
///
struct pos_processor_base
{
    pos_processor_base()
        : _is_skip_process_pos(false) {}

    virtual
    ~pos_processor_base() {}

    void
    check_process_pos(const int stage_no,
                      const pos_t pos)
    {
        if (_is_skip_process_pos) return;
        process_pos(stage_no,pos);
    }

    virtual
    void
    process_pos(const int stage_no,
                const pos_t pos) = 0;

protected:
    mutable bool _is_skip_process_pos;
};
