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

#include "manta/EdgeInfo.hh"
#include "manta/SVCandidate.hh"
#include "manta/SVCandidateData.hh"
#include "manta/SVLocusSet.hh"
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

    virtual
    void
    writeHeader(const char* progVersion)
    {
        writeHeaderPrefix(progVersion);
        writeHeaderSuffix();
    }

protected:
    void
    writeHeaderPrefix(const char* progVersion);

    void
    writeHeaderSuffix();

    void
    writeTranslocPair(
            const EdgeInfo& edge,
            const unsigned svIndex,
            const SVCandidate& sv);

    virtual
    void
    modifyInfo(
            const bool /*isFirstOfPair*/,
            std::vector<std::string>& /*infotags*/)
    {}

private:

    void
    writeTransloc(
        const SVBreakend& bp1,
        const SVBreakend& bp2,
        const std::string& idPrefix,
        const bool isFirstOfPair);

protected:
    const std::string& _referenceFilename;
    const unsigned _minPairCount;

private:
    const bam_header_info& _header;
    std::ostream& _os;

    boost::format _idFormatter;
};

