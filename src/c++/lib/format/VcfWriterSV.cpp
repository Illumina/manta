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
    _os << "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n";
    _os << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n";
    _os << "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">\n";
    _os << "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">\n";
    _os << "##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakend\">\n";
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
    const SVBreakend& bp1,
    const SVBreakend& bp2,
    const std::string& idPrefix,
    const bool isFirstOfPair,
    const SVCandidateData& svData )
{
    std::vector<std::string> infotags;

    // get CHROM
    const std::string& chrom(_header.chrom_data[bp1.interval.tid].label);
    const std::string& mate_chrom(_header.chrom_data[bp2.interval.tid].label);

    const known_pos_range2& bp1range(bp1.interval.range);
    const known_pos_range2& bp2range(bp2.interval.range);

    // get POS
    const pos_t pos(bp1range.center_pos()+1);
    const pos_t mate_pos(bp2range.center_pos()+1);

    // get ID
    const std::string localId(idPrefix + (isFirstOfPair ? '0' : '1'));
    const std::string mateId(idPrefix + (isFirstOfPair ? '1' : '0'));

    // get REF
    std::string ref;
    get_standardized_region_seq(_referenceFilename,chrom,pos-1,pos-1,ref);

    assert(1 == ref.size());

    // build alt:
    boost::format altFormat("%4%%3%%1%:%2%%3%%5%");
    {
        std::string altPrefix;
        std::string altSuffix;
        if     (bp1.state == SVBreakendState::RIGHT_OPEN)
        {
            altPrefix=ref;
        }
        else if (bp1.state == SVBreakendState::LEFT_OPEN)
        {
            altSuffix=ref;
        }
        else
        {
            assert(! "Unexpected bp2.state");
        }


        char altSep('?');
        if     (bp2.state == SVBreakendState::RIGHT_OPEN)
        {
            altSep=']';
        }
        else if (bp2.state == SVBreakendState::LEFT_OPEN)
        {
            altSep='[';
        }
        else
        {
            assert(! "Unexpected bp2.state");
        }

        altFormat % mate_chrom % mate_pos % altSep % altPrefix % altSuffix;
    }

    // build INFO field
    infotags.push_back("SVTYPE=BND");
    infotags.push_back("MATEID="+mateId);
    infotags.push_back( str(boost::format("BND_PAIR_SUPPORT=%i") % bp1.readCount) );
    infotags.push_back( str(boost::format("PAIR_SUPPORT=%i") % bp1.pairCount) );
    if (! bp1.isPrecise())
    {
        infotags.push_back("IMPRECISE");
        infotags.push_back( str( boost::format("CIPOS=%i,%i") % (bp1range.begin_pos()+1) % bp1range.end_pos() ) );
    }

    modifyInfo(bp1,bp2,isFirstOfPair,infotags,svData);

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
    const SVCandidateData& svData)
{
    const std::string idPrefix( str(_idFormatter % edge.locusIndex % edge.nodeIndex1 % edge.nodeIndex2 % svIndex ) );

    writeTransloc(sv.bp1, sv.bp2, idPrefix, true, svData);
    writeTransloc(sv.bp2, sv.bp1, idPrefix, false, svData);
}


void
VcfWriterSV::
writeInvdel(
    const SVCandidate& sv,
    const std::string& label)
{
    const bool isBp1First(sv.bp1.interval.range.end_pos()<sv.bp2.interval.range.begin_pos());

    const SVBreakend& bpA(isBp1First ? sv.bp1 : sv.bp2);
    const SVBreakend& bpB(isBp1First ? sv.bp2 : sv.bp1);

    std::vector<std::string> infotags;

    // get CHROM
    const std::string& chrom(_header.chrom_data[sv.bp1.interval.tid].label);

    const known_pos_range2& bpArange(bpA.interval.range);
    const known_pos_range2& bpBrange(bpB.interval.range);

    // get POS
    const pos_t pos(bpArange.center_pos()+1);
    if (pos<1) return;

    const pos_t endPos(bpBrange.center_pos()+1);

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
    if ((! bpA.isPrecise()) || (! bpB.isPrecise()))
    {
        infotags.push_back("IMPRECISE");
    }

    if (! bpA.isPrecise())
    {
        infotags.push_back( str( boost::format("CIPOS=%i,%i") % (bpArange.begin_pos()+1) % bpArange.end_pos() ) );
    }
    if (! bpB.isPrecise())
    {
        infotags.push_back( str( boost::format("CIEND=%i,%i") % (bpBrange.begin_pos()+1) % bpBrange.end_pos() ) );
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
    const unsigned svIndex,
    const SVCandidate& sv,
    const SVCandidateData& svData)
{
    const SV_TYPE::index_t svType(getSVType(sv));
    if      (svType == SV_TYPE::INTERTRANSLOC)
    {
        writeTranslocPair(edge, svIndex, sv, svData);
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

