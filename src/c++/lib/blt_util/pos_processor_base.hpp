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

#include "blt_util/blt_types.hpp"

/// \brief Interface for objects designed to perform work in a single pass over a position range
///
/// Work progress is communicated via process_pos(). This object is designed to
/// link its child with the stage_manager.
///
struct pos_processor_base {
  pos_processor_base() : _is_skip_process_pos(false) {}

  virtual ~pos_processor_base() {}

  void check_process_pos(const int stage_no, const pos_t pos)
  {
    if (_is_skip_process_pos) return;
    process_pos(stage_no, pos);
  }

  /// Execute position dependent logic associated with a particular stage
  /// in a positional processing pipeline.
  ///
  /// Stages are each offset by a fixed value from the HEAD position in the
  /// positional pipeline. Conventions on stage numbering are assumed to be
  /// enforced the a separate stage_manager object.
  ///
  virtual void process_pos(const int stage_no, const pos_t pos) = 0;

protected:
  mutable bool _is_skip_process_pos;
};
