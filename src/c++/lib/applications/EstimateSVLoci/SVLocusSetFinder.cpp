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

/// \file
/// \author Chris Saunders
///

#include "SVLocusSetFinder.hpp"

#include "blt_util/log.hpp"
#include "common/Exceptions.hpp"
#include "htsapi/align_path_bam_util.hpp"
#include "manta/ChromDepthFilterUtil.hpp"
#include "manta/ReadFilter.hpp"

#include <iostream>
#include <sstream>

/// The compression level used by this object's private depth-per-position buffer.
static const unsigned depthBufferCompression = 16;

SVLocusSetFinder::SVLocusSetFinder(
    const ESLOptions&                               opt,
    const GenomeInterval&                           scanRegion,
    const std::shared_ptr<reference_contig_segment> refSeqPtr,
    std::shared_ptr<SVLocusSet>                     svLociPtr)
  : _isAlignmentTumor(opt.alignFileOpt.isAlignmentTumor),
    _scanRegion(scanRegion),
    _refSeqPtr(refSeqPtr),
    _svLociPtr(svLociPtr),
    _positionReadDepthEstimatePtr(std::make_shared<depth_buffer_compressible>(depthBufferCompression)),
    _regionManager(_scanRegion, _svLociPtr, _positionReadDepthEstimatePtr),
    _readScanner(opt.scanOpt, opt.statsFilename, opt.alignFileOpt.alignmentFilenames, opt.isRNA),
    _isMaxDepthFilter(false),
    _maxDepth(0)
{
  assert(_svLociPtr);

  //
  // initialize max depth filtration values:
  //
  const ChromDepthFilterUtil dFilter(opt.chromDepthFilename, opt.scanOpt.maxDepthFactor, _bamHeader());
  _isMaxDepthFilter = dFilter.isMaxDepthFilter();
  if (_isMaxDepthFilter) {
    _maxDepth = dFilter.maxDepth(scanRegion.tid);
  }

  assert(opt.alignFileOpt.alignmentFilenames.size() == _svLociPtr->getAllSampleReadCounts().size());
}

void SVLocusSetFinder::addToDepthBuffer(const bam_record& bamRead)
{
  if (!_isMaxDepthFilter) return;

  const pos_t refPos(bamRead.pos() - 1);

  // Estimated read depth uses a very simple approximation that each input read aligns without any indels.
  // This is done to reduce the depth estimation runtime overhead.
  const pos_t readSize(bamRead.read_size());
  _positionReadDepthEstimatePtr->inc(refPos, readSize);
}

void SVLocusSetFinder::update(
    const stream_state_reporter& streamErrorReporter,
    const bam_record&            bamRead,
    const unsigned               defaultReadGroupIndex)
{
  // This is the primary read filtration step for the purpose of graph building
  //
  // Although unmapped reads themselves are filtered out, reads with unmapped mates are still
  // accepted, because these contribute signal for assembly (indel) regions.
  if (isReadUnmappedOrFilteredCore(bamRead)) return;

  const bool isTumor(_isAlignmentTumor[defaultReadGroupIndex]);
  if (!isTumor) addToDepthBuffer(bamRead);

  // Filter out reads from high-depth chromosome regions
  if (_isMaxDepthFilter) {
    if (_positionReadDepthEstimatePtr->val(bamRead.pos() - 1) > _maxDepth) return;
  }

  // Verify the input reads conform to BAM input restrictions before testing if they could be input evidence
  // The reason this is done here is that the restriction against the '=' sign in the sequence could otherwise
  // silently alter downstream evidence checks and assembly steps:
  {
    const bam_seq  readSeq(bamRead.get_bam_read());
    const unsigned readSize(readSeq.size());
    for (unsigned baseIndex(0); baseIndex < readSize; baseIndex++) {
      if (readSeq.get_code(baseIndex) == BAM_BASE::REF) {
        std::ostringstream oss;
        oss << "Unsupported use of the '=' symbol in the BAM/CRAM SEQ field from read:\n";
        streamErrorReporter.report_state(oss);
        BOOST_THROW_EXCEPTION(illumina::common::GeneralException(oss.str()));
      }
    }
  }

  // inputCounts is part of the statistics tracking framework used for
  // methods diagnostics. It does not impact the graph build.
  SampleReadInputCounts& inputCounts(_getLocusSet().getSampleReadInputCounts(defaultReadGroupIndex));

  // Filter out reads below the minimum MAPQ threshold
  if (bamRead.map_qual() < _readScanner.getMinMapQ()) {
    inputCounts.minMapq++;
    return;
  }

  // Filter out reads which are not found to be indicative of an SV
  //
  // For certain types of SV evidence, these tests are designed to be a faster approximation of the
  // full test which will be repeated when the read(-pairs) are converted to SV locus graph elements,
  // therefor a small number of reads which make it through this filter will be effectively filtered
  // later as non-evidence -- they will not be filtered per se, but rather converted into zero SVLocus
  // objects.
  //
  if (!_readScanner.isSVEvidence(bamRead, defaultReadGroupIndex, _refSeq(), &(inputCounts.evidenceCount)))
    return;

#ifdef DEBUG_SFINDER
  log_os << __FUNCTION__ << ": Accepted read. isNonCompressedAnomalous " << isNonCompressedAnomalous
         << " is Local assm evidence: " << isLocalAssemblyEvidence << " read: " << bamRead << "\n";
#endif

  // check that this read starts in our scan region:
  if (!_scanRegion.range.is_pos_intersect(bamRead.pos() - 1)) return;

  // QC check of read length
  SVLocusScanner::checkReadSize(streamErrorReporter, bamRead);

  // update the region manager to move the head pointer forward to the current read's position
  //
  // This may trigger a denoising step or clear buffered information for all positions at a
  // given offset below the new position
  //
  _regionManager.handle_new_pos_value(bamRead.pos() - 1);

  // convert the given read into zero to many SVLocus objects
  //
  // In almost all cases, each read should be converted into one SVLocus object, and that
  // SVLocus will consist of either one or two SVLocusNodes.
  //
  std::vector<SVLocus> loci;

  // evidenceCounts is part of the statistics tracking framework used for
  // methods diagnostics. It does not impact the graph build.
  SampleEvidenceCounts& evidenceCounts(_getLocusSet().getSampleEvidenceCounts(defaultReadGroupIndex));
  _readScanner.getSVLoci(bamRead, defaultReadGroupIndex, _bamHeader(), _refSeq(), loci, evidenceCounts);

  // merge each non-empty SV locus into this genome segment graph:
  for (const SVLocus& locus : loci) {
    if (locus.empty()) continue;
    _getLocusSet().merge(locus);
  }
}
