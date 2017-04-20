//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2017 Illumina, Inc.
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

#include "extractAssemblyReads.hh"
#include "htsapi/bam_streamer.hh"
#include "manta/ShadowReadFinder.hh"
#include "manta/SVLocusScanner.hh"



void
extractAssemblyReadsFromBam(
    const ReadScannerOptions& scanOpt,
    const AssemblerOptions& asmOpt,
    const char* bamFile,
    AssemblyReadInput& reads)
{
    bam_streamer bamStream(bamFile);

    ShadowReadFinder shadow(scanOpt.minSingletonMapqCandidates);

    while (bamStream.next())
    {
        const bam_record& bamRead(*(bamStream.get_record_ptr()));

        // filter out reads we ALWAYS filter out from manta
        //
        // don't filter out MAPQ0 because the split reads tend to have reduced mapping scores:
        if (SVLocusScanner::isReadFilteredCore(bamRead)) continue;

        if (bamRead.isNonStrictSupplement()) continue;

        const bool isShadowKeeper(shadow.check(bamRead));

        // only keep unmapped shadows????
        ///TODO --- is this appropriate for this tool?
        if ((not isShadowKeeper) and bamRead.is_unmapped()) continue;

        bool isReversed(false);

        // if shadow read, determine if we need to reverse:
        if (isShadowKeeper)
        {
            if (bamRead.is_mate_fwd_strand())
            {
                isReversed = (! isReversed);
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

            for (unsigned i(0); i < size; ++i)
            {
                if (qual[i] < minQval) nread[i] = 'N';
            }
        }

        if (isReversed) reverseCompStr(reads.back());
    }
}
