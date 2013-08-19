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
/// \author Chris Saunders
///

#include "format/VcfWriterSV.hh"

#include "blt_util/samtools_fasta_util.hh"
#include "blt_util/string_util.hh"
#include "blt_util/seq_util.hh"
#include "blt_util/vcf_util.hh"

#include <iostream>


VcfWriterSV::
VcfWriterSV(
    const std::string& referenceFilename,
    const SVLocusSet& set,
    std::ostream& os) :
    _referenceFilename(referenceFilename),
    _minPairCount(set.getMinMergeEdgeCount()),
    _header(set.header),
    _os(os),
    _idFormatter("MantaBND:%i:%i:%i:%i:")
{
}



void
VcfWriterSV::
writeHeaderPrefix(
    const char* progName,
    const char* progVersion)
{
    _os << "##fileformat=VCFv4.1\n";
    _os << "##fileData=" << vcf_fileDate << "\n";
    _os << "##source=" << progName << " " << progVersion << "\n";
    _os << "##reference=file://" << _referenceFilename << "\n";

    BOOST_FOREACH(const bam_header_info::chrom_info& cdata, _header.chrom_data)
    {
        _os << "##contig=<ID=" << cdata.label << ",length=" << cdata.length << ">\n";
    }

    /// vcf 4.1 reserved/suggested INFO tags:
    _os << "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n";
    _os << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n";
    _os << "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">\n";
    _os << "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">\n";
    _os << "##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakend\">\n";
#if 0
    _os << "##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description=\"Length of base pair identical micro-homology at event breakpoints\">\n";
    _os << "##INFO=<ID=HOMSEQ,Number=.,Type=String,Description=\"Sequence of base pair identical micro-homology at event breakpoints\">\n";
#endif

    /// custom INFO tags:
    _os << "##INFO=<ID=SVINSLEN,Number=.,Type=Integer,Description=\"Length of micro-insertion at event breakpoints\">\n";
    _os << "##INFO=<ID=SVINSSEQ,Number=.,Type=String,Description=\"Sequence of micro-insertion at event breakpoints\">\n";
    _os << "##INFO=<ID=PAIR_SUPPORT,Number=1,Type=Integer,Description=\"Read pairs supporting this variant where both reads are confidently mapped\">\n";
    _os << "##INFO=<ID=BND_PAIR_SUPPORT,Number=1,Type=Integer,Description=\"Confidently mapped reads supporting this variant at this breakend (mapping may not be confident at remote breakend)\">\n";
    _os << "##INFO=<ID=UPSTREAM_PAIR_SUPPORT,Number=1,Type=Integer,Description=\"Confidently mapped reads supporting this variant at the upstream breakend (mapping may not be confident at downstream breakend)\">\n";
    _os << "##INFO=<ID=DOWNSTREAM_PAIR_SUPPORT,Number=1,Type=Integer,Description=\"Confidently mapped reads supporting this variant at this downstream breakend (mapping may not be confident at upstream breakend)\">\n";

    addHeaderInfo();

    addHeaderFilters();

    _os << "##ALT=<ID=BND,Description=\"Translocation Breakend\">\n";
    _os << "##ALT=<ID=INV,Description=\"Inversion\">\n";
    _os << "##ALT=<ID=DEL,Description=\"Deletion\">\n";
    _os << "##ALT=<ID=DUP:TANDEM,Description=\"Tandem Duplication\">\n";
    _os << "##ALT=<ID=COMPLEX,Description=\"Unknown Candidate Type\">\n";
}



void
VcfWriterSV::
writeHeaderSuffix()
{
    _os << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
}




static
void
makeInfoField(
    const std::vector<std::string>& info,
    std::ostream& os)
{
    static const char sep(';');
    bool isFirst(true);
    BOOST_FOREACH(const std::string& is, info)
    {
        if (! isFirst) os << sep;
        else          isFirst = false;
        os << is;
    }
}



void
VcfWriterSV::
writeTransloc(
    const SVCandidate& sv,
    const std::string& idPrefix,
    const bool isFirstBreakend,
    const SVCandidateSetData& svData,
    const SVCandidateAssemblyData& adata)
{
    const bool isImprecise(sv.isImprecise());
    const bool isBreakendRangeSameShift(sv.isBreakendRangeSameShift());

    const SVBreakend& bpA( isFirstBreakend ? sv.bp1 : sv.bp2);
    const SVBreakend& bpB( isFirstBreakend ? sv.bp2 : sv.bp1);

    std::vector<std::string> infotags;

    // get CHROM
    const std::string& chrom(_header.chrom_data[bpA.interval.tid].label);
    const std::string& mateChrom(_header.chrom_data[bpB.interval.tid].label);

    const known_pos_range2& bpArange(bpA.interval.range);
    const known_pos_range2& bpBrange(bpB.interval.range);

    if (! isImprecise)
    {
        assert(bpArange.size() == bpBrange.size());
    }

    // get POS
    pos_t pos(bpArange.center_pos()+1);
    pos_t matePos(bpBrange.center_pos()+1);
    if (! isImprecise)
    {
        pos = bpArange.begin_pos()+1;
        if (isBreakendRangeSameShift)
        {
            matePos = bpBrange.begin_pos()+1;
        }
        else
        {
            matePos = bpBrange.end_pos();
        }
    }

    // get ID
    const std::string localId(idPrefix + (isFirstBreakend ? '0' : '1'));
    const std::string mateId(idPrefix + (isFirstBreakend ? '1' : '0'));

    // get REF
    std::string ref;
    get_standardized_region_seq(_referenceFilename,chrom,pos-1,pos-1,ref);

    assert(1 == ref.size());

    // build alt:
    boost::format altFormat("%4%%3%%1%:%2%%3%%5%");
    {
        std::string altPrefix;
        std::string altSuffix;
        if     (bpA.state == SVBreakendState::RIGHT_OPEN)
        {
            altPrefix = ref + sv.insertSeq;
        }
        else if (bpA.state == SVBreakendState::LEFT_OPEN)
        {
            altSuffix = sv.insertSeq + ref;
        }
        else
        {
            assert(! "Unexpected bpA.state");
        }


        char altSep('?');
        if     (bpB.state == SVBreakendState::RIGHT_OPEN)
        {
            altSep=']';
        }
        else if (bpB.state == SVBreakendState::LEFT_OPEN)
        {
            altSep='[';
        }
        else
        {
            assert(! "Unexpected bpB.state");
        }

        altFormat % mateChrom % matePos % altSep % altPrefix % altSuffix;
    }

    // build INFO field
    infotags.push_back("SVTYPE=BND");
    infotags.push_back("MATEID="+mateId);
    infotags.push_back( str(boost::format("BND_PAIR_SUPPORT=%i") % bpA.readCount) );
    infotags.push_back( str(boost::format("PAIR_SUPPORT=%i") % bpA.pairCount) );
    if (isImprecise)
    {
        infotags.push_back("IMPRECISE");
    }

    if (bpArange.size() > 1)
    {
        infotags.push_back( str( boost::format("CIPOS=%i,%i") % ((bpArange.begin_pos()+1) - pos) % (bpArange.end_pos() - pos) ));
#if 0
        infotags.push_back( str( boost::format("HOMLEN=%i") % (bpArange.size()) ));
#endif
    }


    if (! sv.insertSeq.empty())
    {
        infotags.push_back( str( boost::format("SVINSLEN=%i") % (sv.insertSeq.size()) ));
        if (isFirstBreakend || (bpA.state != bpB.state))
        {
            infotags.push_back( str( boost::format("SVINSSEQ=%s") % (sv.insertSeq) ));
        }
        else
        {
            infotags.push_back( str( boost::format("SVINSSEQ=%s") % reverseCompCopyStr(sv.insertSeq) ));
        }
    }

    modifyInfo(bpA,bpB,isFirstBreakend,svData, adata, infotags);

    // write out record:
    _os << chrom
        << '\t' << pos
        << '\t' << localId // ID
        << '\t' << ref // REF
        << '\t' << str( altFormat ) // ALT
        << '\t' << '.' // QUAL
        << '\t' << getFilter() // FILTER
        << '\t';
    makeInfoField(infotags,_os); // INFO
    _os << '\n';
}



void
VcfWriterSV::
writeTranslocPair(
    const EdgeInfo& edge,
    const unsigned svIndex,
    const SVCandidate& sv,
    const SVCandidateSetData& svData,
    const SVCandidateAssemblyData& adata)
{
    const std::string idPrefix( str(_idFormatter % edge.locusIndex % edge.nodeIndex1 % edge.nodeIndex2 % svIndex ) );

    writeTransloc(sv, idPrefix, true, svData, adata);
    writeTransloc(sv, idPrefix, false, svData, adata);
}


void
VcfWriterSV::
writeInvdel(
    const SVCandidate& sv,
    const std::string& label)
{
    const bool isImprecise(sv.isImprecise());
    const bool isBreakendRangeSameShift(sv.isBreakendRangeSameShift());

    const bool isBp1First(sv.bp1.interval.range.end_pos()<sv.bp2.interval.range.begin_pos());

    const SVBreakend& bpA(isBp1First ? sv.bp1 : sv.bp2);
    const SVBreakend& bpB(isBp1First ? sv.bp2 : sv.bp1);

    std::vector<std::string> infotags;

    // get CHROM
    const std::string& chrom(_header.chrom_data[sv.bp1.interval.tid].label);

    const known_pos_range2& bpArange(bpA.interval.range);
    const known_pos_range2& bpBrange(bpB.interval.range);

    if (! isImprecise)
    {
        assert(bpArange.size() == bpBrange.size());
    }

    // get POS and endPos
    pos_t pos(bpArange.center_pos()+1);
    pos_t endPos(bpBrange.center_pos()+1);
    if (! isImprecise)
    {
        pos = bpArange.begin_pos()+1;
        if (isBreakendRangeSameShift)
        {
            endPos = bpBrange.begin_pos()+1;
        }
        else
        {
            endPos = bpBrange.end_pos();
        }
    }

    if (pos<1) return;

    // get ID
    static const std::string localId(".");

    // get REF
    std::string ref;
    get_standardized_region_seq(_referenceFilename,chrom,pos-1,pos-1,ref);

    assert(1 == ref.size());

    // build alt:
    const std::string alt( str( boost::format("<%s>") % label));

    // build INFO field
    std::vector<std::string> words;
    split_string(label,':',words);
    infotags.push_back( str(boost::format("SVTYPE=%s") % words[0]));
    infotags.push_back( str(boost::format("END=%i") % endPos));
    infotags.push_back( str(boost::format("SVLEN=%i") % (endPos-pos+1)));
    infotags.push_back( str(boost::format("UPSTREAM_PAIR_SUPPORT=%i") % bpA.readCount) );
    infotags.push_back( str(boost::format("DOWNSTREAM_PAIR_SUPPORT=%i") % bpB.readCount) );
    infotags.push_back( str(boost::format("PAIR_SUPPORT=%i") % bpA.pairCount) );
    if (isImprecise)
    {
        infotags.push_back("IMPRECISE");
    }

    if (bpArange.size() > 1)
    {
        infotags.push_back( str( boost::format("CIPOS=%i,%i") % ((bpArange.begin_pos()+1) - pos) % (bpArange.end_pos() - pos) ));
    }
    if (bpBrange.size() > 1)
    {
        infotags.push_back( str( boost::format("CIEND=%i,%i") % ((bpBrange.begin_pos()+1) - endPos) % (bpBrange.end_pos() - endPos) ));
    }

#if 0
    if (bpArange.size() > 1)
    {
        infotags.push_back( str( boost::format("HOMLEN=%i") % (bpArange.size()) ));
    }
#endif

    if (! sv.insertSeq.empty())
    {
        infotags.push_back( str( boost::format("SVINSLEN=%i") % (sv.insertSeq.size()) ));
        if (isBp1First || (bpA.state != bpB.state))
        {
            infotags.push_back( str( boost::format("SVINSSEQ=%s") % (sv.insertSeq) ));
        }
        else
        {
            infotags.push_back( str( boost::format("SVINSSEQ=%s") % reverseCompCopyStr(sv.insertSeq) ));
        }
    }

    //modifyInfo(isFirstOfPair, infotags);

    // write out record:
    _os << chrom
        << '\t' << pos
        << '\t' << localId // ID
        << '\t' << ref // REF
        << '\t' << alt // ALT
        << '\t' << '.' // QUAL
        << '\t' << getFilter() // FILTER
        << '\t';
    makeInfoField(infotags,_os); // INFO
    _os << '\n';

}



void
VcfWriterSV::
writeInversion(
    const SVCandidate& sv)
{
    writeInvdel(sv,"INV");
}



void
VcfWriterSV::
writeDeletion(
    const SVCandidate& sv)
{
    writeInvdel(sv,"DEL");
}


void
VcfWriterSV::
writeTanDup(
    const SVCandidate& sv)
{
    writeInvdel(sv,"DUP:TANDEM");
}


void
VcfWriterSV::
writeComplex(
    const SVCandidate& sv)
{
    writeInvdel(sv,"COMPLEX");
}



void
VcfWriterSV::
writeSVCore(
    const EdgeInfo& edge,
    const SVCandidateSetData& svData,
    const SVCandidateAssemblyData& adata,
    const unsigned svIndex,
    const SVCandidate& sv)
{
    const SV_TYPE::index_t svType(getSVType(sv));
    if      (svType == SV_TYPE::INTERTRANSLOC)
    {
        writeTranslocPair(edge, svIndex, sv, svData, adata);
    }
    else if (svType == SV_TYPE::INVERSION)
    {
        writeInversion(sv);
    }
    else if (svType == SV_TYPE::DELETION)
    {
        writeDeletion(sv);
    }
    else if (svType == SV_TYPE::TANDUP)
    {
        writeTanDup(sv);
    }
    else if (svType == SV_TYPE::COMPLEX)
    {
        writeComplex(sv);
    }
    else
    {
        assert(! "Unknown SV type");
    }
}

