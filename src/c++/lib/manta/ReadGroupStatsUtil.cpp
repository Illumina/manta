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
/// \author Bret Barnes, Xiaoyu Chen
///

#include "manta/ReadGroupStatsUtil.hpp"

#include "blt_util/ReadKey.hpp"
#include "blt_util/log.hpp"
#include "common/Exceptions.hpp"
#include "htsapi/align_path_bam_util.hpp"
#include "htsapi/bam_record_util.hpp"
#include "htsapi/bam_streamer.hpp"
#include "manta/ReadFilter.hpp"
#include "manta/ReadGroupLabel.hpp"

#include <array>
#include <iostream>
#include <set>
#include <sstream>
#include <vector>

//#define DEBUG_RPS

/// compare distributions to determine stats convergence
static bool isStatSetMatch(const SizeDistribution& pss1, const SizeDistribution& pss2)
{
  static const float cdfPrecision(0.001f);

  for (float prob(0.05f); prob < 1; prob += 0.1f) {
    // check if percentile values equal
    if (std::abs(pss1.quantile(prob) - pss2.quantile(prob)) >= 1) {
      return false;
    }

    // check the convergence of fragsize cdf
    const int fragSize(pss2.quantile(prob));
    if (std::abs(pss1.cdf(fragSize) - pss2.cdf(fragSize)) >= cdfPrecision) {
      return false;
    }
  }

  return true;
}

/// This produces a useful result only when both reads align to the same
/// chromosome.
static PAIR_ORIENT::index_t getRelOrient(const bam_record& br)
{
  pos_t pos1           = br.pos();
  bool  is_fwd_strand1 = br.is_fwd_strand();
  pos_t pos2           = br.mate_pos();
  bool  is_fwd_strand2 = br.is_mate_fwd_strand();

  if (!br.is_first()) {
    std::swap(pos1, pos2);
    std::swap(is_fwd_strand1, is_fwd_strand2);
  }

  return PAIR_ORIENT::get_index(pos1, is_fwd_strand1, pos2, is_fwd_strand2);
}

/// given an input integer, return an integer with all but the highest 4 decimal digits set to zero
///
/// this method is not written effeciently, and not intended for general integer truncation.
/// it is used as part of a simple compression scheme for the fragment sizes of the frag size
/// distribution
///
static unsigned getSimplifiedFragSize(unsigned fragmentSize)
{
  unsigned fragSize(fragmentSize);

  // reduce fragsize resolution for very large sizes:
  unsigned steps(0);
  while (fragSize > 1000) {
    fragSize /= 10;
    steps++;
  }
  for (unsigned stepIndex(0); stepIndex < steps; ++stepIndex) fragSize *= 10;
  return fragSize;
}

/// get insert size from bam record removing refskip (e.g. spliced) segments
static int getFragSizeMinusSkip(const bam_record& bamRead)
{
  using namespace ALIGNPATH;
  ALIGNPATH::path_t _apath;
  bam_cigar_to_apath(bamRead.raw_cigar(), bamRead.n_cigar(), _apath);

  int fragSize(std::abs(bamRead.template_size()));
  if (fragSize == 0) return 0;
  for (const path_segment& ps : _apath) {
    if (ps.type == SKIP) fragSize -= ps.length;
  }

  if (fragSize <= 0) {
    using namespace illumina::common;

    std::ostringstream oss;
    oss << "Unexpected fragment size (" << fragSize << ") deduced from bam record: " << bamRead << "\n"
        << "\tPossible invalid template size in bam record.";
    BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
  }

  return fragSize;
}

/// Does this read contain any refskip operation
static bool hasRefSkip(const bam_record& bamRead)
{
  using namespace ALIGNPATH;
  ALIGNPATH::path_t _apath;
  bam_cigar_to_apath(bamRead.raw_cigar(), bamRead.n_cigar(), _apath);

  for (const path_segment& ps : _apath) {
    if (ps.type == SKIP) return true;
  }
  return false;
}

/// track pair orientation so that a consensus can be found for a read group
///
struct ReadGroupOrientTracker {
  ReadGroupOrientTracker(const char* bamLabel, const char* rgLabel)
    : _isFinalized(false), _totalOrientCount(0), _rgLabel(bamLabel, rgLabel)
  {
    std::fill(_orientCount.begin(), _orientCount.end(), 0);
  }

  void addOrient(const PAIR_ORIENT::index_t ori)
  {
    static const unsigned maxOrientCount(100000);
    if (_totalOrientCount >= maxOrientCount) return;
    if (ori == PAIR_ORIENT::UNKNOWN) return;
    addOrientImpl(ori);
  }

  const ReadPairOrient& getConsensusOrient(const ReadCounter& readCounter)
  {
    finalize(readCounter);
    return _finalOrient;
  }

  unsigned getMinCount()
  {
    static const unsigned minCount(100);
    return minCount;
  }

  bool isOrientCountGood() { return (_totalOrientCount >= getMinCount()); }

private:
  void addOrientImpl(const PAIR_ORIENT::index_t ori)
  {
    assert(!_isFinalized);
    assert(ori < PAIR_ORIENT::SIZE);

    _orientCount[ori]++;
    _totalOrientCount++;
  }

  void finalize(const ReadCounter& readCounter)
  {
    if (_isFinalized) return;
    bool     isMaxIndex(false);
    unsigned maxIndex(0);
    for (unsigned i(0); i < _orientCount.size(); ++i) {
      if ((!isMaxIndex) || (_orientCount[i] > _orientCount[maxIndex])) {
        isMaxIndex = true;
        maxIndex   = i;
      }
    }

    assert(isMaxIndex);

    _finalOrient.setVal(maxIndex);

    {
      // make sure there's a dominant consensus orientation and that we have a minimum number of samples:
      static const float minMaxFrac(0.9f);

      using namespace illumina::common;

      if (!isOrientCountGood()) {
        std::ostringstream oss;
        oss << "Too few high-confidence read pairs (" << _totalOrientCount
            << ") to determine pair orientation for " << _rgLabel << "'\n"
            << "\tAt least " << getMinCount()
            << " high-confidence read pairs are required to determine pair orientation.\n"
            << readCounter << "\n";

        BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
      }

      const unsigned minMaxCount(static_cast<unsigned>(minMaxFrac * _totalOrientCount));
      if (_orientCount[maxIndex] < minMaxCount) {
        const unsigned     maxPercent((_orientCount[maxIndex] * 100) / _totalOrientCount);
        std::ostringstream oss;
        oss << "Can't determine consensus pair orientation of " << _rgLabel << ".\n"
            << "\tThe most frequent orientation is '" << _finalOrient << "' (" << maxPercent << "% of "
            << _totalOrientCount << " total used read pairs)\n"
            << "\tThe fraction of '" << _finalOrient
            << "' among total high-confidence read pairs needs to be more than " << minMaxFrac
            << " to determine consensus pair orientation.\n"
            << readCounter << "\n";

        BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
      }
    }

    _isFinalized = true;
  }

  bool                                    _isFinalized;
  unsigned                                _totalOrientCount;
  const ReadGroupLabel                    _rgLabel;
  std::array<unsigned, PAIR_ORIENT::SIZE> _orientCount;
  ReadPairOrient                          _finalOrient;
};

struct SimpleRead {
  SimpleRead(PAIR_ORIENT::index_t ort, unsigned sz) : _orient(ort), _insertSize(sz) {}

  PAIR_ORIENT::index_t _orient;
  unsigned             _insertSize;
};

struct ReadGroupBuffer {
  ReadGroupBuffer() : _abnormalRpCount(0), _observationRpCount(0) {}

  void updateBuffer(PAIR_ORIENT::index_t ort, unsigned sz)
  {
    _readInfo.emplace_back(ort, sz);

    if (ort == PAIR_ORIENT::Rp) {
      _observationRpCount++;
      if (sz >= 5000) _abnormalRpCount++;
    }
  }

  bool isBufferFull() const { return (_observationRpCount >= 1000); }

  bool isBufferNormal() const
  {
    if (_observationRpCount == 0) return false;

    return ((_abnormalRpCount / (float)_observationRpCount) < 0.01);
  }

  unsigned getAbnormalCount() { return _abnormalRpCount; }

  unsigned getObservationCount() { return _observationRpCount; }

  const std::vector<SimpleRead>& getBufferedReads() { return _readInfo; }

  void clearBuffer()
  {
    _abnormalRpCount    = 0;
    _observationRpCount = 0;
    _readInfo.clear();
  }

private:
  unsigned                _abnormalRpCount;
  unsigned                _observationRpCount;
  std::vector<SimpleRead> _readInfo;
};

/// all data required to build ReadGroupStats during estimation from the bam file
///
/// ultimately the only information we want to keep is the ReadGroupStats object itself,
/// which can be exported from this object
///
struct ReadGroupTracker {
  explicit ReadGroupTracker(
      const char*        bamLabel             = nullptr,
      const char*        rgLabel              = nullptr,
      const std::string& defaultStatsFilename = "")
    : _isFinalized(false),
      _rgLabel(bamLabel, rgLabel),
      _orientInfo(bamLabel, rgLabel),
      _isChecked(false),
      _isInsertSizeConverged(false),
      _defaultStatsFilename(defaultStatsFilename)
  {
  }

  unsigned insertSizeObservations() const { return _stats.fragStats.totalObservations(); }

  unsigned getMinObservationCount() const
  {
    static const unsigned minObservations(100);
    return minObservations;
  }

  bool isObservationCountGood() const { return (insertSizeObservations() >= getMinObservationCount()); }

  void checkInsertSizeCount()
  {
    static const unsigned statsCheckCnt(100000);
    const bool            isCheck((insertSizeObservations() % statsCheckCnt) == 0);
    if (isCheck) _isChecked = true;
  }

  bool isInsertSizeChecked() const { return _isChecked; }

  void clearChecked() { _isChecked = false; }

  bool isInsertSizeConverged() const { return _isInsertSizeConverged; }

  bool isCheckedOrConverged() const { return (_isChecked || isInsertSizeConverged()); }

  void updateInsertSizeConvergenceTest()
  {
    // check convergence
    if (_oldInsertSize.totalObservations() > 0) {
      _isInsertSizeConverged = isStatSetMatch(_oldInsertSize, _stats.fragStats);
    }
    _oldInsertSize = _stats.fragStats;
  }

  ReadGroupBuffer& getBuffer() { return _buffer; }

  void addBufferedData()
  {
    for (const SimpleRead& srd : _buffer.getBufferedReads()) {
      // get orientation stats before final filter for innie reads below:
      //
      // we won't use anything but innie reads for insert size stats, but sampling
      // orientation beforehand allows us to detect, ie. a mate-pair library
      // so that we can blow-up with an informative error msg
      //
      const PAIR_ORIENT::index_t ori(srd._orient);
      addOrient(ori);

      // we define "high-confidence" read pairs as those reads passing all filters
      addHighConfidenceReadPairCount();

      // filter mapped innies on the same chrom
      //
      // note we don't rely on the proper pair bit because this already contains an
      // arbitrary length filter and  subjects the method to aligner specific variation
      //
      // TODO: ..note this locks-in standard ilmn orientation -- okay for now but this function needs major
      // re-arrangement for mate-pair support, we could still keep independence from each aligner's proper
      // pair decisions by estimating a fragment distro for each orientation and only keeping the one with
      // the most samples
      //
      if (ori != PAIR_ORIENT::Rp) continue;
      addInsertSize(srd._insertSize);
    }
  }

  /// Add one observation to the buffer
  /// If the buffer is full, AND if the fragment size distribution in the buffer looks normal, add the buffered data;
  /// otherwise, discard the buffer and move to the next region
  bool addObservation(PAIR_ORIENT::index_t ori, unsigned sz)
  {
    bool isNormal(true);

    _buffer.updateBuffer(ori, sz);

    if (_buffer.isBufferFull()) {
      // check abnormal fragment-size distribution in the buffer
      if (_buffer.isBufferNormal()) {
        addBufferedData();
        checkInsertSizeCount();
      } else {
        isNormal = false;
#ifdef DEBUG_RPS
        std::cerr << "The previous region (buffered) contains too many abnormal reads. "
                  << "abnormalCount=" << _buffer.getAbnormalCount()
                  << " observationCount=" << _buffer.getObservationCount() << "\n";
#endif
      }
      _buffer.clearBuffer();
    }

    return isNormal;
  }

  const ReadGroupOrientTracker& getOrientInfo() const { return _orientInfo; }

  /// getting a const ref of the stats forces finalization steps:
  const ReadGroupStats& getStats() const
  {
    assert(_isFinalized);
    return _stats;
  }

  void addReadCount() { _stats.readCounter.addReadCount(); }

  void addPairedReadCount() { _stats.readCounter.addPairedReadCount(); }

  void addUnpairedReadCount() { _stats.readCounter.addUnpairedReadCount(); }

  void addPairedLowMapqReadCount() { _stats.readCounter.addPairedLowMapqReadCount(); }

  void addHighConfidenceReadPairCount() { _stats.readCounter.addHighConfidenceReadPairCount(); }

  void finalize()
  {
    if (_isFinalized) return;

    // add the remaining data in the buffer
    if (_buffer.isBufferNormal()) {
      addBufferedData();
    }
    _buffer.clearBuffer();

    if ((_defaultStatsFilename != "") && (!_orientInfo.isOrientCountGood() || !isObservationCountGood())) {
      log_os << "Can't generate pair statistics for " << _rgLabel << "\n"
             << "\tTotal high-confidence read pairs (FR) used for insert size estimation: "
             << insertSizeObservations() << "\n"
             << "\tAt least " << getMinObservationCount()
             << " high-confidence read pairs (FR) are required to estimate insert size.\n"
             << "\tUsing existing stats as default: " << _defaultStatsFilename << "\n"
             << _stats.readCounter << "\n";
      ReadGroupStatsSet defaultStats;
      defaultStats.load(_defaultStatsFilename.c_str());
      // The current implementation for defaultStats assumes that the
      // provided defaultStats file contains only 1 ReadGroupStatsSet
      // and uses this for every _rgLabel for which stats determination
      // fails. However, in more sophisticated scenarios including
      // variant calling of paired tumor-normal samples, the developer
      // may want to provide a functionality of using a defaultStats
      // file containing multiple ReadGroupStatsSet which will be chosen
      // based on origin of the sample.
      _stats       = defaultStats.getStats(0);
      _isFinalized = true;
      return;
    }

    // finalize pair orientation:
    _stats.relOrients = _orientInfo.getConsensusOrient(_stats.readCounter);

    if (_stats.relOrients.val() != PAIR_ORIENT::Rp) {
      using namespace illumina::common;

      std::ostringstream oss;
      oss << "Unexpected consensus read orientation (" << _stats.relOrients << ") for " << _rgLabel << "\n"
          << "\tManta currently handles paired-end (FR) reads only.\n";
      BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
    }

    // finalize insert size distro:
    if (!isInsertSizeConverged()) {
      if (!isObservationCountGood()) {
        using namespace illumina::common;

        std::ostringstream oss;
        oss << "Can't generate pair statistics for " << _rgLabel << "\n"
            << "\tTotal high-confidence read pairs (FR) used for insert size estimation: "
            << insertSizeObservations() << "\n"
            << "\tAt least " << getMinObservationCount()
            << " high-confidence read pairs (FR) are required to estimate insert size.\n"
            << _stats.readCounter << "\n";
        BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
      } else if (!isInsertSizeChecked()) {
        updateInsertSizeConvergenceTest();
      }

      if (!isInsertSizeConverged()) {
        log_os << "WARNING: read pair statistics did not converge for " << _rgLabel << "\n"
               << "\tTotal high-confidence read pairs (FR) used for insert size estimation: "
               << insertSizeObservations() << "\n"
               << _stats.readCounter << "\n";
      }
    }

    // final step before saving is to cut-off the extreme end of the fragment size distribution, this
    // is similar the some aligner's proper-pair bit definition of (3x the standard mean, etc.)
    static const float filterQuant(0.9995f);
    _stats.fragStats.filterObservationsOverQuantile(filterQuant);

    _isFinalized = true;
  }

private:
  void addOrient(const PAIR_ORIENT::index_t ori)
  {
    assert(!_isFinalized);

    _orientInfo.addOrient(ori);
  }

  void addInsertSize(const int size)
  {
    assert(!_isFinalized);

    _stats.fragStats.addObservation(size);
  }

  bool                   _isFinalized;
  const ReadGroupLabel   _rgLabel;
  ReadGroupOrientTracker _orientInfo;

  bool             _isChecked;
  bool             _isInsertSizeConverged;
  SizeDistribution _oldInsertSize;  // previous fragment distribution is stored to determine convergence

  ReadGroupBuffer   _buffer;
  ReadGroupStats    _stats;
  const std::string _defaultStatsFilename;
};

struct ReadAlignFilter {
  /// use only the most conservative alignments to generate fragment stats --
  /// filter reads containing any cigar types besides MATCH with optional
  /// trailing SOFTCLIP and single REFSKIP
  bool isFilterRead(const bam_record& bamRead)
  {
    using namespace ALIGNPATH;

    bam_cigar_to_apath(bamRead.raw_cigar(), bamRead.n_cigar(), _apath);

    if (!bamRead.is_fwd_strand()) std::reverse(_apath.begin(), _apath.end());

    bool isMatched(false);
    bool isSkip(false);
    bool isClipped(false);
    for (const path_segment& ps : _apath) {
      if (is_segment_align_match(ps.type)) {
        if (isClipped) return true;
        isMatched = true;
      } else if (ps.type == SKIP) {
        if (isSkip) return true;
        isSkip = true;
      } else if (ps.type == SOFT_CLIP) {
        isClipped = true;
      } else {
        return true;
      }
    }
    return (!isMatched);
  }

private:
  ALIGNPATH::path_t _apath;
};

/// this object handles filtration of reads which are:
///
/// 1. not downstream or (if both reads start at same position) not the second in order in the bam file
/// 2. part of a high depth pileup region
///
struct ReadPairDepthFilter {
  bool isFilterRead(const bam_record& bamRead)
  {
    static const unsigned maxPosCount(1);

    if (bamRead.target_id() != _lastTargetId) {
      _goodMates.clear();
      _lastTargetId = bamRead.target_id();
      _posCount     = 0;
      _lastPos      = bamRead.pos();
    } else if (bamRead.pos() != _lastPos) {
      _posCount = 0;
      _lastPos  = bamRead.pos();
    }

    // Assert only two reads per fragment
    const unsigned readNum(bamRead.is_first() ? 1 : 2);
    assert(bamRead.is_second() == (readNum == 2));

    // Filter pairs with templateSize 0 (unknown)
    if (bamRead.template_size() == 0) return true;

    // sample each read pair once by sampling stats from
    // downstream read only, or whichever read is encountered
    // second if the read and its mate start at the same position:
    const bool isDownstream(bamRead.pos() > bamRead.mate_pos());
    const bool isSamePos(bamRead.pos() == bamRead.mate_pos());

    if (isDownstream || isSamePos) {
      const int     mateReadNo(bamRead.is_first() ? 2 : 1);
      const ReadKey mateKey(bamRead.qname(), mateReadNo, false);

      mateMap_t::iterator i(_goodMates.find(mateKey));

      if (i == _goodMates.end()) {
        if (isDownstream) return true;
      } else {
        _goodMates.erase(i);
        return false;
      }
    }

    // to prevent high-depth pileups from overly biasing the
    // read stats, we only take maxPosCount read pairs from each start
    // pos. by not inserting a key in goodMates, we also filter
    // the downstream mate:
    if (_posCount >= maxPosCount) return true;
    ++_posCount;

    // crude mechanism to manage total set memory
    static const unsigned maxMateSetSize(100000);
    if (_goodMates.size() > maxMateSetSize) _goodMates.clear();

    // Ignore pairs where the upstream mate has a refskip, since we cannot
    // compute the correct insert size later when looking at the downstream mate
    // (Or we would have to save the total refskip length here)
    if (hasRefSkip(bamRead)) return true;

    _goodMates.insert(ReadKey(bamRead));
    return true;
  }

private:
  typedef std::set<ReadKey> mateMap_t;

  unsigned  _posCount = 0;
  mateMap_t _goodMates;

  int _lastTargetId = 0;
  int _lastPos      = 0;
};

struct CoreInsertStatsReadFilter {
  bool isFilterRead(const bam_record& bamRead)
  {
    // filter common categories of undesirable reads:
    if (isReadFilteredCore(bamRead)) return true;

    if (bamRead.isNonStrictSupplement()) return true;

    if (!is_mapped_chrom_pair(bamRead)) return true;
    if (bamRead.map_qual() == 0) return true;

    // filter any split reads with an SA tag:
    if (bamRead.isSASplit()) return true;

    // remove alignments other than {X}M({Z}N{X2}M)?({Y}S)? (or reverse for reverse strand)
    if (alignFilter.isFilterRead(bamRead)) return true;

    // filter out upstream reads and high depth regions:
    if (pairFilter.isFilterRead(bamRead)) return true;

    return false;
  }

  ReadAlignFilter     alignFilter;
  ReadPairDepthFilter pairFilter;
};

#if 0
/// samples around various short segments of the genome so
/// that stats generation isn't biased towards a single region
struct GenomeSampler
{
    GenomeSampler(
        const std::string& bamFile) :
        _readStream(bamFile.c_str()),
        _chromCount(0),
        _isInitRegion(false),
        _isInitRecord(false),
        _isActiveChrom(true),
        _currentChrom(0)
    {
        const bam_hdr_t& header(*_readStream.get_header());

        _chromCount = (header.n_targets);

        _isActiveChrom = (_chromCount > 0);

        _chromSize.resize(_chromCount,0);
        _chromHighestPos.resize(_chromCount,-1);

        for (int32_t i(0); i<_chromCount; ++i)
        {
            _chromSize[i] = (header.target_len[i]);
        }

        init();
    }

    /// advance to next region of the genome, this must be called before using nextRecord():
    bool
    nextRegion()
    {
        if (! _isActiveChrom) return false;

        _isInitRegion=true;
        _isInitRecord=false;
        return true;
    }

    /// advance to next record in the current region, this must be called before using getBamRecord():
    bool
    nextRecord()
    {
        assert(_isInitRegion);

        if (! _isActiveChrom) return false;

        _isInitRecord=true;
        return true;
    }

    /// access current bam record
    const bam_record&
    getBamRecord()
    {
        assert(_isInitRegion);
        assert(_isInitRecord);

        return *(_readStream.get_record_ptr());
    }

private:
    void
    init()
    {
        if (! _isActiveChrom) return;

        const int32_t startPos(_chromHighestPos[chromIndex]+1);

        _currentChrom=0;
    }

    bam_streamer _readStream;
    int32_t _chromCount;
    std::vector<int32_t> _chromSize;

    std::vector<int32_t> _chromHighestPos;
    bool _isInitRegion;
    bool _isInitRecord;
    bool _isActiveChrom; ///< this is used to track whether we've reached the end of all chromosomes
    int32_t _currentChrom;
};
#endif

/// manage the info structs for each RG
struct ReadGroupManager {
  typedef std::map<ReadGroupLabel, ReadGroupTracker> RGMapType;

  explicit ReadGroupManager(const std::string& statsBamFile)
    : _isFinalized(false), _statsBamFile(statsBamFile)
  {
  }

  ReadGroupTracker& getTracker(const bam_record& bamRead, const std::string& defaultStatsFilename = "")
  {
    const char* readGroup(getReadGroup(bamRead));

    return getTracker(readGroup, defaultStatsFilename);
  }

  ReadGroupTracker& getTracker(const char* readGroup, const std::string& defaultStatsFilename = "")
  {
    ReadGroupLabel rgKey(_statsBamFile.c_str(), readGroup, false);

    RGMapType::iterator rgIter(_rgTracker.find(rgKey));
    if (rgIter == _rgTracker.end()) {
      std::pair<RGMapType::iterator, bool> retval;
      retval = _rgTracker.insert(std::make_pair(
          ReadGroupLabel(_statsBamFile.c_str(), readGroup),
          ReadGroupTracker(_statsBamFile.c_str(), readGroup, defaultStatsFilename)));

      assert(retval.second);
      rgIter = retval.first;
    }

    return rgIter->second;
  }

  /// check if all read groups have been sufficiently sampled for a slice
  /// for each read group, either 100k samples has been collected, or the insert size distrubution has
  /// converged
  bool isFinishedSlice()
  {
    for (RGMapType::value_type& val : _rgTracker) {
      if (!val.second.isCheckedOrConverged()) return false;
    }

    for (RGMapType::value_type& val : _rgTracker) {
      val.second.clearChecked();
    }

    return true;
  }

  /// test if all read groups have converged or hit other stopping conditions
  bool isStopEstimation()
  {
    static const unsigned maxRecordCount(5000000);
    for (RGMapType::value_type& val : _rgTracker) {
      if (!(val.second.isInsertSizeConverged() || (val.second.insertSizeObservations() > maxRecordCount)))
        return false;
    }
    return true;
  }

  const RGMapType& getMap()
  {
    finalize();
    return _rgTracker;
  }

private:
  void finalize()
  {
    if (_isFinalized) return;
    for (RGMapType::value_type& val : _rgTracker) {
      val.second.finalize();
    }
    _isFinalized = true;
  }

  bool              _isFinalized;
  const std::string _statsBamFile;
  RGMapType         _rgTracker;
};

void extractReadGroupStatsFromAlignmentFile(
    const std::string& referenceFilename,
    const std::string& alignmentFilename,
    const std::string& defaultStatsFilename,
    ReadGroupStatsSet& rstats)
{
  bam_streamer read_stream(alignmentFilename.c_str(), referenceFilename.c_str());

  const bam_hdr_t&     header(read_stream.get_header());
  const int32_t        chromCount(header.n_targets);
  std::vector<int32_t> chromSize(chromCount, 0);
  std::vector<int32_t> chromHighestPos(chromCount, -1);
  for (int32_t i(0); i < chromCount; ++i) {
    chromSize[i] = (header.target_len[i]);
  }

  bool isStopEstimation(false);
  bool isActiveChrom(true);

  CoreInsertStatsReadFilter coreFilter;
  ReadGroupManager          rgManager(alignmentFilename.c_str());

#ifndef READ_GROUPS
  static const char defaultReadGroup[] = "";
  ReadGroupTracker& rgInfo(rgManager.getTracker(defaultReadGroup, defaultStatsFilename));
#endif

  while (isActiveChrom && (!isStopEstimation)) {
    isActiveChrom = false;
    for (int32_t chromIndex(0); chromIndex < chromCount; ++chromIndex) {
      if (isStopEstimation) break;

      // keep sampling until
      // either the chromosome has been exhuasted
      // or the current chunk has been sufficiently sampled
      bool isFinishedSlice(false);
      while (!isFinishedSlice) {
        const int32_t startPos(chromHighestPos[chromIndex] + 1);
        if (startPos >= chromSize[chromIndex]) break;

#ifdef DEBUG_RPS
        std::cerr << "INFO: Stats requesting bam region starting from: chrid: " << chromIndex
                  << " start: " << startPos << "\n";
#endif
        read_stream.resetRegion(chromIndex, startPos, chromSize[chromIndex]);

        while (read_stream.next()) {
          const bam_record& bamRead(*(read_stream.get_record_ptr()));
          if (bamRead.pos() < startPos) continue;

          chromHighestPos[chromIndex] = bamRead.pos();
          isActiveChrom               = true;

          rgInfo.addReadCount();
          if (bamRead.is_paired()) {
            rgInfo.addPairedReadCount();
            if (bamRead.map_qual() == 0) {
              rgInfo.addPairedLowMapqReadCount();
            }
          } else {
            rgInfo.addUnpairedReadCount();
          }

          if (coreFilter.isFilterRead(bamRead)) continue;

#ifdef READ_GROUPS
          ReadGroupTracker& rgInfo(rgManager.getTracker(bamRead, defaultStatsFilename));
#endif
          if (rgInfo.isInsertSizeConverged()) continue;

          const PAIR_ORIENT::index_t ori(getRelOrient(bamRead));
          unsigned                   fragSize(0);
          if (ori == PAIR_ORIENT::Rp) {
            fragSize = getSimplifiedFragSize(getFragSizeMinusSkip(bamRead));
          }
          const bool isNormal = rgInfo.addObservation(ori, fragSize);

          if (!isNormal) {
            chromHighestPos[chromIndex] += std::max(1, chromSize[chromIndex] / 100);
#ifdef DEBUG_RPS
            std::cerr << " Jump to chrid: " << chromIndex << " position: " << chromHighestPos[chromIndex]
                      << "\n";
#endif
            break;
          }

          if (!rgInfo.isInsertSizeChecked()) continue;

          // check convergence
          rgInfo.updateInsertSizeConvergenceTest();

          isFinishedSlice = rgManager.isFinishedSlice();
          if (!isFinishedSlice) continue;

          isStopEstimation = rgManager.isStopEstimation();

          // break from reading the current chromosome
          break;
        }

        // move to next region if no read falling in the current region
        if (chromHighestPos[chromIndex] <= startPos) {
#ifdef DEBUG_RPS
          std::cerr << "No read found in the previous region.\n";
#endif
          chromHighestPos[chromIndex] += std::max(1, chromSize[chromIndex] / 100);
        }
      }
    }
  }

  for (const ReadGroupManager::RGMapType::value_type& val : rgManager.getMap()) {
    rstats.setStats(val.first, val.second.getStats());
  }
}
