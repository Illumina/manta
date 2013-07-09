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

#include "manta/VcfWriterSomaticSV.hh"


void
VcfWriterSomaticSV::
modifyInfo(
    const bool isFirstOfPair,
    std::vector<std::string>& infotags)
{
    assert(_ssInfoPtr != NULL);
    const SomaticSVScoreInfo& ssInfo(*_ssInfoPtr);

    infotags.push_back( str(boost::format("SOMATIC_SCORE=%i") % ssInfo.somaticScore) );
    infotags.push_back( str(boost::format("NORMAL_PAIR_SUPPORT=%i") % ssInfo.normal.spanPairs) );
    infotags.push_back( str(boost::format("TUMOR_PAIR_SUPPORT=%i") % ssInfo.tumor.spanPairs) );
    infotags.push_back( str(boost::format("NORMAL_BND_PAIR_SUPPORT=%i") %
                            (isFirstOfPair ? ssInfo.normal.bp1SpanReads : ssInfo.normal.bp2SpanReads) ) );
    infotags.push_back( str(boost::format("TUMOR_BND_PAIR_SUPPORT=%i") %
                            (isFirstOfPair ? ssInfo.tumor.bp1SpanReads : ssInfo.tumor.bp2SpanReads) ) );
}



void
VcfWriterSomaticSV::
writeSV(
    const EdgeInfo& edge,
    const SVCandidateData& ,
    const unsigned svIndex,
    const SVCandidate& sv,
    const SomaticSVScoreInfo& ssInfo)
{

    if (ssInfo.somaticScore < _somaticOpt.minOutputSomaticScore) return;

    //TODO: this is a lame way to customize subclass behavior:
    _ssInfoPtr=&ssInfo;
    writeTranslocPair(edge, svIndex, sv);
    _ssInfoPtr=NULL;
}
