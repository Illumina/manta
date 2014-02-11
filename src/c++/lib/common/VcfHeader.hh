// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Refactored by Richard Shaw from PUMA by Mauricio Varea
///

/*****************************************************************************/

#pragma once

#include <string>
#include <vector>

#include "VcfFields.hh"
#include "ChromosomeMetadata.hh"
#include "VcfMetaInformation.hh"

/*****************************************************************************/

class VcfHeader
{
    // TODO: check if "const char *" is better than "char *"
    //typedef boost::unordered_map<const char *, size_t, boost::hash<const char *>, Vcf::equal_to<const char *> > Index;
public:
    VcfHeader()
        : fileFormat_(""), fileDate_(""), source_(""), reference_(""), phasing_("")
        , chrList_(0), infoList_(0), filterList_(0), formatList_(0), altList_(0)
        , blockCompressed_(false)
        , strict_(true) {}
    VcfHeader& self()
    {
        return *this;
    }
public:
    const std::string& getFileFormat() const
    {
        return fileFormat_;
    }
    const std::string& getFileDate() const
    {
        return fileDate_;
    }
    const std::string& getSource() const
    {
        return source_;
    }
    const std::string& getReference() const
    {
        return reference_;
    }
    const std::string& getPhasing() const
    {
        return phasing_;
    }

    size_t getContigCount() const
    {
        return chrList_.size();
    }
    size_t getInfoCount() const
    {
        return infoList_.size();
    }
    size_t getFilterCount() const
    {
        return filterList_.size();
    }
    size_t getFormatCount() const
    {
        return formatList_.size();
    }
    size_t getAltCount() const
    {
        return altList_.size();
    }

    size_t getContigIndex(const char* name) const;

    //size_t getInfoIndex(const std::string &name) const {
    //    return getInfoIndex(name.c_str());}

    size_t getInfoIndex(const char* name) const;

    //size_t getFilterIndex(const std::string &name) const {
    //    return getFilterIndex(name.c_str());}

    size_t getFilterIndex(const char* name) const;

    //size_t getFormatIndex(const std::string &name) const {
    //     return getFormatIndex(name.c_str());}

    size_t getFormatIndex(const char* name) const;

    //size_t getAltIndex(const std::string &name) const {
    //   return getAltIndex(name.c_str());}

    size_t getAltIndex(const char* name) const;

    const ChromosomeMetadata& getContig(size_t i) const
    {
        return chrList_.at(i);
    }

    const VcfMetaInformation& getInfo(size_t i) const
    {
        return infoList_.at(i);
    }

    const VcfMetaInformation& getFilter(size_t i) const
    {
        return filterList_.at(i);
    }

    const VcfMetaInformation& getFormat(size_t i) const
    {
        return formatList_.at(i);
    }

    const VcfMetaInformation& getAlt(size_t i) const
    {
        return altList_.at(i);
    }

    const std::vector<ChromosomeMetadata>& getContigList() const
    {
        return chrList_;
    }

    const std::vector<VcfMetaInformation>& getInfoList() const
    {
        return infoList_;
    }

    const std::vector<VcfMetaInformation>& getFilterList() const
    {
        return filterList_;
    }

    const std::vector<VcfMetaInformation>& getFormatList() const
    {
        return formatList_;
    }

    const std::vector<VcfMetaInformation>& getAltList() const
    {
        return altList_;
    }

    unsigned numSamples() const
    {
        return sampleNameList_.size();
    }

    const std::vector<std::string>& getSampleNameList() const
    {
        return sampleNameList_;
    }

    //size_t addContig(const ChromosomeMetadata &chr);

    void clear();
    bool isBlockCompressed() const
    {
        return blockCompressed_;
    }
    bool isStrict() const
    {
        return strict_;
    }
    bool hasContigList() const
    {
        return !chrList_.empty();
    }
private:
    std::string fileFormat_;
    std::string fileDate_;
    std::string source_;
    std::string reference_;
    std::string phasing_;

    std::vector<ChromosomeMetadata> chrList_;
    std::vector<VcfMetaInformation> infoList_;
    std::vector<VcfMetaInformation> filterList_;
    std::vector<VcfMetaInformation> formatList_;
    std::vector<VcfMetaInformation> altList_;

    std::vector<std::string> sampleNameList_;

    bool blockCompressed_;
    bool strict_;

    template< class T>
    size_t getIndex(const char* field, const std::vector<T>& idx,
                    const char* name) const;
    friend std::ostream& operator<<(std::ostream& os,
                                    const VcfHeader& vcfHeader);
    friend std::istream& operator>>(std::istream& is, VcfHeader& vcfHeader);
};

/*****************************************************************************/

std::istream& operator>>(std::istream& is, VcfHeader& vcfHeader);
std::ostream& operator<<(std::ostream& os, const VcfHeader& vcfHeader);

/// specialize FileHeader for VCF use case
// typedef FileHeader<VcfHeader> AnnotatedVcfHeader;

/*****************************************************************************/
