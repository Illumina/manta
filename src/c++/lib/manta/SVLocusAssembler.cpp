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

///
/// \author Ole Schulz-Trieglaff
///


#include "blt_util/log.hh"
#include "manta/SVLocusAssembler.hh"

#include "boost/foreach.hpp"


// compile with this macro to get verbose output:
//#define DEBUG_ASBL

#ifdef DEBUG_ASBL
#include <iostream>
std::ostream& dbg_os(std::cerr);
#endif



SVLocusAssembler::
SVLocusAssembler(
    const GSCOptions& opt,
    const SmallAssemblerOptions& assembleOpt) :
    _scanOpt(opt.scanOpt),
    _assembleOpt(assembleOpt),
    _readScanner(_scanOpt,opt.statsFilename,opt.alignmentFilename)
{
    // setup regionless bam_streams:
    // setup all data for main analysis loop:
    BOOST_FOREACH(const std::string& afile, opt.alignmentFilename)
    {
        // avoid creating shared_ptr temporaries:
        streamPtr tmp(new bam_streamer(afile.c_str()));
        _bamStreams.push_back(tmp);
    }
}



void
SVLocusAssembler::
getBreakendReads(
    const SVBreakend& bp,
    ReadIndexType& readIndex,
    AssemblyReadInput& reads) const
{
    static const size_t minIntervalSize(300);

    known_pos_range2 searchRange;
    if (bp.interval.range.size() >= minIntervalSize)
    {
        searchRange = bp.interval.range;
    }
    else
    {
        const size_t missing = minIntervalSize - bp.interval.range.size();
        assert(missing > 0);
        const size_t wobble = missing/2;
        // FIXME : not sure what happens if (end_pos + wobble) > chromosome size?
        static const size_t zero(0);
        searchRange.set_range(std::max((bp.interval.range.begin_pos()-wobble),zero),(bp.interval.range.end_pos()+wobble));
    }

    static const unsigned minClipLen(3);

    const unsigned bamCount(_bamStreams.size());
    for (unsigned bamIndex(0); bamIndex < bamCount; ++bamIndex)
    {
        bam_streamer& bamStream(*_bamStreams[bamIndex]);

        // set bam stream to new search interval:
        bamStream.set_new_region(bp.interval.tid, searchRange.begin_pos(), searchRange.end_pos());

        static const unsigned MAX_NUM_READS(1000);

		while (bamStream.next() && (reads.size() < MAX_NUM_READS))
		{
			const bam_record& bamRead(*(bamStream.get_record_ptr()));

            // FIXME: add some criteria to filter for "interesting" reads here, for now we add
            // only clipped reads and reads without N
            if ((bamRead.pos()-1) >= searchRange.end_pos()) break;
            if (!_readScanner.isClipped(bamRead) ) continue;
            if (_readScanner.getClipLength(bamRead)<minClipLen ) continue;
            if (bamRead.get_bam_read().get_string().find('N') != std::string::npos) continue;
            //if ( bamRead.pe_map_qual() == 0 ) continue;
            std::string flag = "1";
            if (bamRead.is_second())
            {
                flag = "2";
            }
            const std::string readKey = std::string(bamRead.qname()) + "_" + flag + "_" + boost::lexical_cast<std::string>(bamIndex);

#ifdef DEBUG_ASBL
            ALIGNPATH::path_t apath;
            bam_cigar_to_apath(bamRead.raw_cigar(), bamRead.n_cigar(), apath);
            dbg_os << "Adding " << readKey << " " << apath << " " << bamRead.pe_map_qual() << " " << bamRead.pos() << endl;
            dbg_os << bamRead.get_bam_read().get_string() << endl;
#endif

            if (readIndex.count(readKey) == 0)
            {
                // the API gives us always the sequence w.r.t to the fwd ref, so no need
                // to reverse complement here
                readIndex.insert(std::make_pair(readKey,reads.size()));
                reads.push_back(bamRead.get_bam_read().get_string());
            }
            else
            {
                log_os << "WARNING : Read name collision : " << readKey << "\n"; 
            }
        }
    }
}



void
SVLocusAssembler::
assembleSingleSVBreakend(const SVBreakend& bp,
                         Assembly& as) const
{
    ReadIndexType readIndex;
    AssemblyReadInput reads;
    getBreakendReads(bp,readIndex, reads);
    AssemblyReadOutput readInfo;
    runSmallAssembler(_assembleOpt, reads, readInfo, as);
}



void
SVLocusAssembler::
assembleSVBreakends(const SVBreakend& bp1,
                    const SVBreakend& bp2,
                    Assembly& as) const
{
    ReadIndexType readIndex;
    AssemblyReadInput reads;
    getBreakendReads(bp1,readIndex, reads);
    getBreakendReads(bp2,readIndex, reads);
    AssemblyReadOutput readInfo;
    runSmallAssembler(_assembleOpt, reads, readInfo, as);
}

