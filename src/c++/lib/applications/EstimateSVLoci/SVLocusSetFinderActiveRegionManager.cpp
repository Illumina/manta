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
/// \author Trevor Ramsay
///

#include "SVLocusSetFinderActiveRegionManager.hpp"

#include "blt_util/depth_buffer.hpp"
#include "svgraph/SVLocusSet.hpp"

#include <cassert>

//#define DEBUG_SFINDER

namespace STAGE {
enum index_t { HEAD, DENOISE, CLEAR_DEPTH };

static stage_data getStageData(const unsigned denoiseRegionProtectedBorderSize)
{
  // the number of positions below the HEAD position where the
  // depth buffer can be safely erased
  static const unsigned clearDepthBorderSize(10);

  stage_data sd;
  sd.add_stage(HEAD);
  sd.add_stage(DENOISE, HEAD, denoiseRegionProtectedBorderSize);
  sd.add_stage(CLEAR_DEPTH, HEAD, clearDepthBorderSize);

  return sd;
}
}  // namespace STAGE

/// \brief Compute the genomic region in which graph denoising is allowed.
///
/// The allowed denoising region is the same as the current scan region, except that the region is shortened
/// by a fixed-size buffer wherever an adjacent genome segment is potentially being analyzed by another
/// process in parallel. This protected zone near adjacent segment boundaries helps to ensure that SV
/// breakend regions spanning the adjacent process boundary are not pessimistically marked as noise before
/// the evidence from both segments has been merged.
///
/// \param scanRegion The region of the genome scanned by this process for SV locus evidence.
///
/// \param bamHeader Bam header information. Used to extract chromosome length here.
///
/// \param denoiseRegionProtectedBorderSize Length of the protected region where denoising is skipped at
/// adjacent process boundaries.
///
/// \return Region where denoising is allowed in this process
static GenomeInterval computeDenoiseRegion(
    const GenomeInterval&  scanRegion,
    const bam_header_info& bamHeader,
    const int              denoiseRegionProtectedBorderSize)
{
  GenomeInterval denoiseRegion = scanRegion;

  known_pos_range2& range(denoiseRegion.range);
  if (range.begin_pos() > 0) {
    range.set_begin_pos(range.begin_pos() + denoiseRegionProtectedBorderSize);
  }

  bool isEndBorder(true);
  if (static_cast<int32_t>(bamHeader.chrom_data.size()) > denoiseRegion.tid) {
    const pos_t chromEndPos(bamHeader.chrom_data[denoiseRegion.tid].length);
    isEndBorder = (range.end_pos() < chromEndPos);
  }

  if (isEndBorder) {
    range.set_end_pos(range.end_pos() - denoiseRegionProtectedBorderSize);
  }

#ifdef DEBUG_SFINDER
  log_os << __FUNCTION__ << ": " << denoiseRegion << "\n";
#endif

  return denoiseRegion;
}

SVLocusSetFinderActiveRegionManager::SVLocusSetFinderActiveRegionManager(
    const GenomeInterval&                      scanRegion,
    std::shared_ptr<SVLocusSet>                svLociPtr,
    std::shared_ptr<depth_buffer_compressible> positionReadDepthEstimatePtr,
    unsigned                                   denoiseRegionProtectedBorderSize)
  : _scanRegion(scanRegion),
    _svLociPtr(svLociPtr),
    _positionReadDepthEstimatePtr(positionReadDepthEstimatePtr),
    _denoiseRegion(
        computeDenoiseRegion(scanRegion, _svLociPtr->getBamHeader(), denoiseRegionProtectedBorderSize)),
    _stageManager(
        STAGE::getStageData(denoiseRegionProtectedBorderSize),
        pos_range(_scanRegion.range.begin_pos(), _scanRegion.range.end_pos()),
        *this),
    _isInDenoiseRegion(false),
    _denoiseStartPos(0)
{
  assert(_svLociPtr);
}

void SVLocusSetFinderActiveRegionManager::process_pos(const int stage_no, const pos_t pos)
{
#ifdef DEBUG_SFINDER
  log_os << __FUNCTION__ << ": stage_no: " << stage_no << " pos: " << pos << "\n";
#endif

  if (stage_no == STAGE::HEAD) {
    // pass
  } else if (stage_no == STAGE::DENOISE) {
    // denoise the SV locus graph in regions of at least this size
    //
    // this parameter batches the denoising process to balance runtime vs.
    // maintaining the SV locus graph as a reasonably compact data structure
    //
    static const pos_t minDenoiseRegionSize(1000);

    if (_denoiseRegion.range.is_pos_intersect(pos)) {
#ifdef DEBUG_SFINDER
      log_os << __FUNCTION__ << ": pos intersect. pos: " << pos << " dnRegion: " << _denoiseRegion
             << " is in region: " << _isInDenoiseRegion << "\n";
#endif

      if (!_isInDenoiseRegion) {
        _denoiseStartPos   = _denoiseRegion.range.begin_pos();
        _isInDenoiseRegion = true;
      }

      if ((1 + pos - _denoiseStartPos) >= minDenoiseRegionSize) {
        _getLocusSet().cleanRegion(GenomeInterval(_denoiseRegion.tid, _denoiseStartPos, (pos + 1)));
        _denoiseStartPos = (pos + 1);
      }
    } else {
#ifdef DEBUG_SFINDER
      log_os << __FUNCTION__ << ": no pos intersect. pos: " << pos << " dnRegion: " << _denoiseRegion
             << " is in region: " << _isInDenoiseRegion << "\n";
#endif

      if (_isInDenoiseRegion) {
        if ((_denoiseRegion.range.end_pos() - _denoiseStartPos) > 0) {
          _getLocusSet().cleanRegion(
              GenomeInterval(_denoiseRegion.tid, _denoiseStartPos, _denoiseRegion.range.end_pos()));
          _denoiseStartPos = _denoiseRegion.range.end_pos();
        }
        _isInDenoiseRegion = false;
      }
    }
  } else if (stage_no == STAGE::CLEAR_DEPTH) {
    _positionReadDepthEstimatePtr->clear_pos(pos);
  } else {
    assert(false && "Unexpected stage id");
  }
}
