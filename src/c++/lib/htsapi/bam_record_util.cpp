//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2018 Illumina, Inc.
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

#include "bam_record_util.hh"
#include "align_path_bam_util.hh"
#include "htsapi/SimpleAlignment_bam_util.hh"

#include <iostream>


bool
is_mapped_pair(
    const bam_record& bam_read)
{
    if (! bam_read.is_paired()) return false;
    if (bam_read.is_unmapped() || bam_read.is_mate_unmapped()) return false;
    return true;
}



bool
is_mapped_chrom_pair(
    const bam_record& bam_read)
{
    if (! is_mapped_pair(bam_read)) return false;
    if (bam_read.target_id() != bam_read.mate_target_id()) return false;
    return true;
}



bool
is_innie_pair(
    const bam_record& bam_read)
{
    if (! is_mapped_chrom_pair(bam_read)) return false;
    if (bam_read.is_fwd_strand() == bam_read.is_mate_fwd_strand()) return false;

    if     (bam_read.pos() < bam_read.mate_pos())
    {
        if (! bam_read.is_fwd_strand()) return false;
    }
    else if (bam_read.pos() > bam_read.mate_pos())
    {
        if (  bam_read.is_fwd_strand()) return false;
    }

    return true;
}


bool
is_possible_adapter_pair(
    const bam_record& bamRead)
{
    if (! is_mapped_chrom_pair(bamRead)) return false;
    if (bamRead.is_fwd_strand() == bamRead.is_mate_fwd_strand()) return false;
    int posDiff(bamRead.mate_pos() - bamRead.pos());
    if (! bamRead.is_fwd_strand())
    {
        posDiff *= -1;
    }
    // Aggressive check: the pos difference between the reverse read and the forwarding read is [-50, 10]
    return (posDiff < 10) && (posDiff > -50);
}


bool
is_adapter_pair(
    const bam_record& bamRead
)
{
    // if the read (primary) contains a SA tag, keep it for assembly
    if (bamRead.isSASplit()) return false;

    using namespace ALIGNPATH;
    const SimpleAlignment aln(getAlignment(bamRead));
    if (bamRead.hasMateCigar())
    {
        // If we have mate cigar information, check for exact adapter condition:
        // Does this read extend (3') past the start (5') of its mate (including softclipping on either)?
        const SimpleAlignment mate(getKnownOrFakedMateAlignment(bamRead));
        if (aln.is_fwd_strand)
        {
            unsigned const endpos = aln.pos + apath_ref_length(aln.path) + apath_soft_clip_trail_size(aln.path);
            unsigned const mateStartPos = mate.pos + apath_ref_length(mate.path) + apath_soft_clip_lead_size(mate.path);
            return (endpos > mateStartPos);
        }
        else
        {
            unsigned const endpos = aln.pos - apath_soft_clip_trail_size(aln.path);
            unsigned const mateStartPos = mate.pos - apath_soft_clip_lead_size(mate.path);
            return (endpos < mateStartPos);
        }
    }
    else
    {
        // If we do not have mate cigar information use an aggressive heuristic:
        // if the read contains soft clip on the 3' end, it likely runs into adapter.
        unsigned const softClipSize(aln.is_fwd_strand ?
                                    apath_soft_clip_trail_size(aln.path) :
                                    apath_soft_clip_lead_size(aln.path));
        return (softClipSize > 0);
    }
}



bool
is_overlapping_pair(
    const bam_record& bamRead,
    const SimpleAlignment& matchedAlignment)
{
    if (! is_mapped_chrom_pair(bamRead)) return false;
    if (bamRead.is_fwd_strand() == bamRead.is_mate_fwd_strand()) return false;

    // get range of alignment after matching all softclip:
    if (bamRead.is_fwd_strand())
    {
        const pos_t matchedEnd(matchedAlignment.pos + apath_ref_length(matchedAlignment.path));
        return (matchedEnd >= bamRead.mate_pos());
    }
    else
    {
        // NOTE: We do not have the mate alignment information here (missing mate end position)
        // so this check uses a very lenient heuristic (overlap mate start position).
        const pos_t matchedBegin(matchedAlignment.pos);
        return (matchedBegin <= bamRead.mate_pos());
    }
}



unsigned
get_avg_quality(
    const bam_record& bam_read)
{
    const unsigned len(bam_read.read_size());
    if (0 == len) return 0;

    const uint8_t* qual(bam_read.qual());
    unsigned sum(0);
    for (unsigned i(0); i<len; ++i)
    {
        sum+=qual[i];
    }
    // this does not capture the decimal remainder but well...
    return (sum/len);
}



static
std::string
getChromName(
    const bam_header_info& bamHeader,
    const int tid)
{
    if (tid >= 0)
    {
        assert(tid < static_cast<int>(bamHeader.chrom_data.size()));
        return bamHeader.chrom_data[tid].label;
    }
    else
    {
        return "UNKNOWN";
    }
}



void
summarizeAlignmentRecord(
    const bam_header_info& bamHeader,
    const bam_record& bamRead,
    std::ostream& os)
{
    if (bamRead.empty())
    {
        os << "NONE";
        return;
    }

    os << bamRead.qname() << "/" << bamRead.read_no()
       << " chrom:pos:strand " << getChromName(bamHeader, bamRead.target_id()) << ":" << bamRead.pos() << ":" << (bamRead.is_fwd_strand() ? '+' : '-');

    ALIGNPATH::path_t apath;
    bam_cigar_to_apath(bamRead.raw_cigar(),bamRead.n_cigar(),apath);
    os << " cigar: " << apath;


    os << " templateSize: " << bamRead.template_size();

    // print SAtag if present:
    static const char satag[] = {'S','A'};
    const char* saStr(bamRead.get_string_tag(satag));
    if (nullptr != saStr)
    {
        os  << " SA: " << saStr;
    }
    if (bamRead.is_secondary())
    {
        os << " secondary";
    }
    if (bamRead.is_supplementary())
    {
        os << " supplementary";
    }

    if (bamRead.is_paired())
    {
        os  << " mateChrom:pos:strand " << getChromName(bamHeader, bamRead.mate_target_id()) << ":" << bamRead.mate_pos() << ":" << (bamRead.is_mate_fwd_strand() ? '+' : '-');
    }
}
