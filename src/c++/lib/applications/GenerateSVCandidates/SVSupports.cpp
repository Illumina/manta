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

/// \file
/// \author Xiaoyu Chen
///

#include "SVSupports.hh"

#include <iostream>

//#define DEBUG_SUPPORT

#ifdef DEBUG_SUPPORT
#include "blt_util/log.hh"
#endif


std::ostream&
operator<<( std::ostream& os, const SupportRead& suppRd)
{
    os << suppRd.tid << ":" << suppRd.pos << "\t";
    for (const auto& sv : suppRd.SVs)
    {
        os << sv.first << ",";
    }

    return os;
}


std::ostream&
operator<<( std::ostream& os, const SupportFragment& suppFrg)
{
    os << suppFrg.read1 << "\n"
       << suppFrg.read2 << "\n";
    return os;
}


std::ostream&
operator<<( std::ostream& os, const SupportFragments& suppFrgs)
{
    for (const auto& frg : suppFrgs.supportFrags)
    {
        os << "qname =" << frg.first << "\n"
           << frg.second;
    }

    return os;
}


std::ostream&
operator<<( std::ostream& os, const SupportSamples& suppSmps)
{
    unsigned size(suppSmps.supportSamples.size());
    for (unsigned i=0; i<size; i++)
    {
        os << "sample index = " << i << "\n"
           << suppSmps.supportSamples[i];
    }

    return os;
}


void
processBamRecords(
    bam_streamer& origBamStream,
    const GenomeInterval& interval,
    const support_fragments_t& supportFrags,
    bam_dumper& bamDumper)
{
#ifdef DEBUG_SUPPORT
    log_os << __FUNCTION__ << "  target interval: "
           << interval << "\n";
#endif

    origBamStream.resetRegion(interval.tid, interval.range.begin_pos(), interval.range.end_pos());
    while (origBamStream.next())
    {
        const bam_record* origBamRec(origBamStream.get_record_ptr());
        bam_record bamRec(*origBamRec);

        const std::string qname(bamRec.qname());
        support_fragments_t::const_iterator suppFragsIter(supportFrags.find(qname));
        if (suppFragsIter != supportFrags.end())
        {
            const SupportFragment& supportFrag(suppFragsIter->second);
            const bool isR1Matched(bamRec.is_first() &&
                                   (bamRec.target_id() == supportFrag.read1.tid) &&
                                   (bamRec.pos() == supportFrag.read1.pos));
            const bool isR2Matched((!bamRec.is_first()) &&
                                   (bamRec.target_id() == supportFrag.read2.tid) &&
                                   (bamRec.pos() == supportFrag.read2.pos));

            if (isR1Matched || isR2Matched)
            {
                const SupportRead& read(isR1Matched ? supportFrag.read1 : supportFrag.read2);
#ifdef DEBUG_SUPPORT
                log_os << __FUNCTION__ << "  matched supporting read: "
                       << read << "\n";
#endif

                bam1_t& br(*(bamRec.get_data()));
                // add new customized field of SV IDs that the read supports
                bool isFirst(true);
                std::string svStr;
                for (const auto& sv : read.SVs)
                {
                    if (! isFirst) svStr.append(",");
                    svStr.append(sv.first);
                    for (const auto& svType : sv.second)
                    {
                        svStr.append('|' + svType);
                    }
                    if (isFirst) isFirst = false;
                }

                static const char svtag[] = {'Z','M'};
                // Do lots of ugly casting on svStr to fit htsapi signature. Usage is actually const in htslib:
                bam_aux_append(&br,svtag,'Z',(svStr.size()+1),
                               reinterpret_cast<uint8_t*>(const_cast<char*>(svStr.c_str())));

                // Update bam record bin value
                bam_update_bin(br);
                // write to bam
                bamDumper.put_record(&br);
            }
        }
    }
}


void
writeSupportBam(const bam_streamer_ptr& origBamStreamPtr,
                const SupportFragments& svSupportFrags,
                const bam_dumper_ptr& supportBamDumperPtr)
{
    std::vector<SupportRead> supportReads;
    const support_fragments_t& supportFrags(svSupportFrags.supportFrags);
    for (const auto& frg : supportFrags)
    {
        supportReads.push_back(frg.second.read1);
        supportReads.push_back(frg.second.read2);
    }
    // sort all the reads w.r.t. genomic positions
    std::sort(supportReads.begin(), supportReads.end());

    // generate a set of intervals containing overlapping reads
    const int readDistance(100);
    int lastTid = -1;
    int lastPos = -1;
    std::vector<GenomeInterval> intervals;
    for (const auto& suppRd : supportReads)
    {
        if  ((lastTid == suppRd.tid) && (lastPos + readDistance >= suppRd.pos))
        {
            GenomeInterval& interval(intervals.back());
            interval.range.set_end_pos(suppRd.pos);
        }
        else
        {
            GenomeInterval interval(suppRd.tid,suppRd.pos-1,suppRd.pos);
            intervals.push_back(interval);
        }

        lastTid=suppRd.tid;
        lastPos=suppRd.pos;
    }

    bam_streamer& origBamStream(*origBamStreamPtr);
    bam_dumper& supportBamDumper(*supportBamDumperPtr);
    for (const auto& interval : intervals)
    {
        processBamRecords(origBamStream, interval,
                          supportFrags, supportBamDumper);
    }

}




