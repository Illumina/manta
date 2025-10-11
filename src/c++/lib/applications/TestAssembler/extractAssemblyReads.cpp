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

#include "extractAssemblyReads.hpp"
#include "htsapi/bam_streamer.hpp"
#include "manta/ReadFilter.hpp"
#include "manta/SVLocusScanner.hpp"
#include "manta/ShadowReadFinder.hpp"

void extractAssemblyReadsFromBam(
    const ReadScannerOptions& scanOpt,
    const AssemblerOptions&   asmOpt,
    const std::string&        referenceFilename,
    const std::string&        alignmentFilename,
    AssemblyReadInput&        reads)
{
  bam_streamer bamStream(alignmentFilename.c_str(), referenceFilename.c_str());

  ShadowReadFinder shadow(scanOpt.minSingletonMapqCandidates);

  while (bamStream.next()) {
    const bam_record& bamRead(*(bamStream.get_record_ptr()));

    // filter out reads we ALWAYS filter out from manta
    //
    // don't filter out MAPQ0 because the split reads tend to have reduced mapping scores:
    if (isReadFilteredCore(bamRead)) continue;

    if (bamRead.isNonStrictSupplement()) continue;

    const bool isShadowKeeper(shadow.check(bamRead));

    // only keep unmapped shadows????
    /// TODO --- is this appropriate for this tool?
    if ((not isShadowKeeper) and bamRead.is_unmapped()) continue;

    bool isReversed(false);

    // if shadow read, determine if we need to reverse:
    if (isShadowKeeper) {
      if (bamRead.is_mate_fwd_strand()) {
        isReversed = (!isReversed);
      }
    }

    reads.push_back(bamRead.get_bam_read().get_string());

    // should we recreate manta's fragmentation of reads at low-quality bases?
    /// TODO --- is this appropriate for this tool?
    const uint8_t minQval(asmOpt.minQval);
    {
      std::string& nread(reads.back());

      const unsigned size(nread.size());
      const uint8_t* qual(bamRead.qual());

      for (unsigned i(0); i < size; ++i) {
        if (qual[i] < minQval) nread[i] = 'N';
      }
    }

    if (isReversed) reverseCompStr(reads.back());
  }
}
