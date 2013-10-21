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
// <https://github.com/sequencing/licenses/>
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
        writeHeaderColumnKey();
    }

    typedef std::vector<std::string> InfoTag_t;
    typedef std::vector<std::pair<std::string,std::vector<std::string> > > SampleTag_t;

protected:
    void
    writeHeaderPrefix(
        const char* progName,
        const char* progVersion);

    void
    writeHeaderColumnKey();

    virtual
    void
    addHeaderInfo() const {}

    virtual
    void
    addHeaderFormat() const {}

    virtual
    void
    addHeaderFilters() const {}

    virtual
    void
    addHeaderFormatSampleKey() const {}

    void
    writeSVCore(
        const EdgeInfo& edge,
        const SVCandidateSetData& svData,
        const SVCandidateAssemblyData& adata,
        const SVCandidate& sv);

    /// add info tags which can be customized by sub-class
    virtual
    void
    modifyInfo(
        InfoTag_t& /*infotags*/) const
    {}

    /// add info tags specific to translocations:
    virtual
    void
    modifyTranslocInfo(
        const bool /*isFirstOfPair*/,
        InfoTag_t& /*infotags*/) const
    {}

    virtual
    void
    writeQual() const
    {
        _os << '.';
    }

    virtual
    void
    writeFilter() const
    {
        _os << '.';
    }

    virtual
    void
    modifySample(
        const SVCandidate& /*sv*/,
        SampleTag_t& /*sampletags*/) const
    {}

    void
    writeFilters(
        const std::set<std::string>& filters) const;

private:

    /// \param[in] idPrefix prefix used for ID/MATEID tags in the vcf ID fields
    /// \param[in] isFirstBreakend if true report bp1, else report bp2
    void
    writeTransloc(
        const SVCandidate& sv,
        const std::string& idPrefix,
        const bool isFirstBreakend,
        const SVCandidateSetData& svData,
        const SVCandidateAssemblyData& adata);

    void
    writeTranslocPair(
        const EdgeInfo& edge,
        const SVCandidate& sv,
        const SVCandidateSetData& svData,
        const SVCandidateAssemblyData& adata);

    /// \param isIndel if true, the variant is a simple right/left breakend insert/delete combination
    void
    writeInvdel(
        const EdgeInfo& edge,
        const SVCandidate& sv,
        const SVCandidateAssemblyData& adata,
        const std::string& label,
        const bool isIndel = false);

    void
    writeInversion(
        const EdgeInfo& edge,
        const SVCandidate& sv,
        const SVCandidateSetData& svData,
        const SVCandidateAssemblyData& adata);

    void
    writeIndel(
        const EdgeInfo& edge,
        const SVCandidate& sv,
        const SVCandidateSetData& svData,
        const SVCandidateAssemblyData& adata);

    void
    writeTanDup(
        const EdgeInfo& edge,
        const SVCandidate& sv,
        const SVCandidateSetData& svData,
        const SVCandidateAssemblyData& adata);

    void
    writeComplex(
        const EdgeInfo& edge,
        const SVCandidate& sv,
        const SVCandidateSetData& svData,
        const SVCandidateAssemblyData& adata);

protected:
    const std::string& _referenceFilename;

private:
    const bam_header_info& _header;
protected:
    std::ostream& _os;
private:
    boost::format _transLocIdFormatter;
    boost::format _otherSVIdFormatter;
};

