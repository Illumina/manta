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
#include "manta/SVCandidateAssemblyData.hh"
#include "manta/SVCandidateSetData.hh"
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
        const SVCandidateSetData& svData,
        const SVCandidateAssemblyData& adata,
        const unsigned svIndex,
        const SVCandidate& sv);

    virtual
    void
    modifyInfo(
        const SVBreakend& /*bp1*/,
        const SVBreakend& /*bp2*/,
        const bool /*isFirstOfPair*/,
        const SVCandidateSetData& /*svData*/,
        const SVCandidateAssemblyData& /*adata*/,
        std::vector<std::string>& /*infotags*/) const
    {}

    virtual
    std::string
    getFilter() const
    {
        return ".";
    }

private:

    /// \param[in] isBp1 if true report bp1, else report bp2
    void
    writeTransloc(
        const SVCandidate& sv,
        const bool isBp1,
        const std::string& idPrefix,
        const bool isFirstOfPair,
        const SVCandidateSetData& svData,
        const SVCandidateAssemblyData& adata);

    void
    writeTranslocPair(
        const EdgeInfo& edge,
        const unsigned svIndex,
        const SVCandidate& sv,
        const SVCandidateSetData& svData,
        const SVCandidateAssemblyData& adata);

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

