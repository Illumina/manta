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

#include <memory>
#include <string>
#include <vector>

#include "blt_util/input_stream_handler.hpp"
#include "htsapi/bam_streamer.hpp"

/// \brief Open all bam files as bam_streamer objects, with no genomic region set
///
/// \param[out] bamStreams Vector of bam_streamers corresponding to the input \p bamFilenames. Existing
/// content will be cleared on input.
void openBamStreams(
    const std::string&                          referenceFilename,
    const std::vector<std::string>&             bamFilenames,
    std::vector<std::shared_ptr<bam_streamer>>& bamStreams);

/// \brief Reset the target genomic region of all \p bamStreams
///
/// Note when \p region is empty this function has no effect
///
void resetBamStreamsRegion(const std::string& region, std::vector<std::shared_ptr<bam_streamer>>& bamStreams);

/// \brief Asserts the headers of all \p bamStreams to see if they are compatible for merged iteration
///
/// Throws an exception if headers are not compatible. Here 'compatible' means each bam streamer uses the same
/// chromosomes in the same order.
///
/// \param[in] bamFilenames Alignment filenames are only used for improving the exception message.
void assertCompatibleBamStreams(
    const std::vector<std::string>&                   bamFilenames,
    const std::vector<std::shared_ptr<bam_streamer>>& bamStreams);

/// \brief Merge bam_streamers together so that they can be iterated as a single merged stream
input_stream_data mergeBamStreams(std::vector<std::shared_ptr<bam_streamer>>& bamStreams);
