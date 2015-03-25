// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2015 Illumina, Inc.
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

#include "manta/EventInfo.hh"
#include "manta/JunctionIdGenerator.hh"
#include "manta/SVCandidate.hh"
#include "manta/SVCandidateAssemblyData.hh"
#include "manta/SVCandidateSetData.hh"
#include "svgraph/SVLocusSet.hh"

#include <iosfwd>


struct VcfWriterSV
{
    VcfWriterSV(
        const std::string& referenceFilename,
        const bool isRNA,
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
        const SVCandidateSetData& svData,
        const SVCandidateAssemblyData& adata,
        const SVCandidate& sv,
        const SVId& svId,
        const EventInfo& event);

    /// add info tags which can be customized by sub-class
    virtual
    void
    modifyInfo(
        const EventInfo& /*event*/,
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

    /// \param[in] isFirstBreakend if true report bp1, else report bp2
    void
    writeTransloc(
        const SVCandidate& sv,
        const SVId& svId,
        const bool isFirstBreakend,
        const SVCandidateSetData& svData,
        const SVCandidateAssemblyData& adata,
        const EventInfo& event);

    void
    writeTranslocPair(
        const SVCandidate& sv,
        const SVId& svId,
        const SVCandidateSetData& svData,
        const SVCandidateAssemblyData& adata,
        const EventInfo& event);

    /// \param isIndel if true, the variant is a simple right/left breakend insert/delete combination
    void
    writeInvdel(
        const SVCandidate& sv,
        const SVId& svId,
        const SVCandidateAssemblyData& adata,
        const bool isIndel,
        const EventInfo& event);

protected:
    const std::string& _referenceFilename;
    const bool _isRNA;

private:
    const bam_header_info& _header;
protected:
    std::ostream& _os;
};

