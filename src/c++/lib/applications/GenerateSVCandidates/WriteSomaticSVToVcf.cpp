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

#include "WriteSomaticSVToVcf.hh"

#include "WriteSVToVcf.hh"

#include "blt_util/samtools_fasta_util.hh"

#include "boost/foreach.hpp"
#include "boost/format.hpp"

#include <iostream>



void
writeSomaticSVVcfHeader(
    const char* referenceFilename,
    const char* version,
    std::ostream& os)
{
    writeSVVcfHeaderPrefix(referenceFilename,version,os);
    writeSVVcfHeaderSuffix(os);
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



static
void
writeTransloc(
    const std::string& referenceFilename,
    const bam_header_info& header,
    const SVBreakend& bp1,
    const SVBreakend& bp2,
    const std::string& idPrefix,
    const bool isFirstOfPair,
    std::ostream& os)
{
    std::vector<std::string> infotags;

    // get CHROM
    const std::string& chrom(header.chrom_data[bp1.interval.tid].label);
    const std::string& mate_chrom(header.chrom_data[bp2.interval.tid].label);

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
    get_standardized_region_seq(referenceFilename,chrom,pos-1,pos-1,ref);

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

    // write out record:
    os << chrom
       << '\t' << pos
       << '\t' << localId // ID
       << '\t' << ref // REF
       << '\t' << str( altFormat ) // ALT
       << '\t' << '.' // QUAL
       << '\t' << '.' // FILTER
       << '\t';
    makeInfoField(infotags,os); // INFO
    os << '\n';
}



void
writeSomaticSVToVcf(
    const std::string& referenceFilename,
    const SVLocusSet& set,
    const EdgeInfo& edge,
    const SVCandidateData& ,
    const unsigned svIndex,
    const SVCandidate& sv,
    const SomaticSVScoreInfo& ssInfo,
    std::ostream& os)
{
    const unsigned minPairCount(set.getMinMergeEdgeCount());
    const bam_header_info& header(set.header);

    boost::format idFormatter("MantaBND:%i:%i:%i:%i:");

    {
        if (sv.bp1.pairCount < minPairCount) return;

        const SV_TYPE::index_t svType(getSVType(sv));
        if (svType == SV_TYPE::INTERTRANSLOC)
        {
            const std::string idPrefix( str(idFormatter % edge.locusIndex % edge.nodeIndex1 % edge.nodeIndex2 % svIndex ) );
            writeTransloc(referenceFilename, header, sv.bp1, sv.bp2, idPrefix, true, os);
            writeTransloc(referenceFilename, header, sv.bp2, sv.bp1, idPrefix, false, os);
        }
        else
        {
            assert(! "Unknown SV type");
        }
    }
}


