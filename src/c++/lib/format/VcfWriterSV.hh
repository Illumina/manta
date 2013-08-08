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

#pragma once

#include "svgraph/EdgeInfo.hh"
#include "manta/SVCandidate.hh"
#include "manta/SVCandidateData.hh"
#include "svgraph/SVLocusSet.hh"
#include "manta/SomaticSVScoreInfo.hh"
#include "options/SomaticCallOptions.hh"
#include "boost/format.hpp"

#include <iosfwd>


struct VcfWriterSV
{
    VcfWriterSV(
        const std::string& referenceFilename,
        const SVLocusSet& set,
        std::ostream& os);

    virtual
    ~VcfWriterSV() {}

    void
    writeHeader(
        const char* progName,
        const char* progVersion)
    {
        writeHeaderPrefix(progName, progVersion);
        writeHeaderSuffix();
    }

protected:
    void
    writeHeaderPrefix(
        const char* progName,
        const char* progVersion);

    void
    writeHeaderSuffix();

    virtual
    void
    addHeaderInfo() const {}

    virtual
    void
    addHeaderFilters() const {}

    void
    writeSVCore(
        const EdgeInfo& edge,
        const unsigned svIndex,
        const SVCandidate& sv,
        const SVCandidateData& svData);

    virtual
    void
    modifyInfo(
    	const SVBreakend& /*bp1*/,
    	const SVBreakend& /*bp2*/,
        const bool /*isFirstOfPair*/,
        std::vector<std::string>& /*infotags*/,
        const SVCandidateData& /*svData*/) const
    {}

    virtual
    std::string
    getFilter() const
    {
        return ".";
    }

private:

    void
    writeTransloc(
        const SVBreakend& bp1,
        const SVBreakend& bp2,
        const std::string& idPrefix,
        const bool isFirstOfPair,
        const SVCandidateData& svData);

    void
    writeTranslocPair(
        const EdgeInfo& edge,
        const unsigned svIndex,
        const SVCandidate& sv,
        const SVCandidateData& svData);

    void
    writeInvdel(
        const SVCandidate& sv,
        const std::string& label);

    void
    writeInversion(
        const SVCandidate& sv);

    void
    writeDeletion(
        const SVCandidate& sv);

    void
    writeTanDup(
        const SVCandidate& sv);

    void
    writeComplex(
        const SVCandidate& sv);

protected:
    const std::string& _referenceFilename;
    const unsigned _minPairCount;

private:
    const bam_header_info& _header;
protected:
    std::ostream& _os;
private:
    boost::format _idFormatter;
};

