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
