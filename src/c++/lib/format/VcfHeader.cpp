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
/// \author Refactored by Richard Shaw from PUMA by Mauricio Varea
///

/*****************************************************************************/

#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/throw_exception.hpp>

#include "blt_util/log.hh"

#include "common/Exceptions.hh"
using namespace illumina::common;

#include "String.hh"
#include "VcfHeader.hh"

/*****************************************************************************/

void VcfHeader::clear()
{
    this->chrList_.clear();
    this->infoList_.clear();
    this->filterList_.clear();
    this->formatList_.clear();
    this->altList_.clear();
}

/*****************************************************************************/

size_t VcfHeader::getContigIndex(const char* name) const
{
    return getIndex<ChromosomeMetadata>("contig", chrList_, name);
}

/*****************************************************************************/

size_t VcfHeader::getInfoIndex(const char* name) const
{
    return getIndex<VcfMetaInformation>( VcfFields::FIXED[7],
                                         infoList_, name );
}

/*****************************************************************************/

size_t VcfHeader::getFilterIndex(const char* name) const
{
    return getIndex<VcfMetaInformation>( VcfFields::FIXED[6],
                                         filterList_, name );
}

/*****************************************************************************/

size_t VcfHeader::getFormatIndex(const char* name) const
{
    return getIndex<VcfMetaInformation>( VcfFields::GENOTYPE[0],
                                         formatList_, name );
}

/*****************************************************************************/

size_t VcfHeader::getAltIndex(const char* name) const
{
    return getIndex<VcfMetaInformation>( VcfFields::FIXED[4],
                                         altList_, name );
}

/*****************************************************************************/

template< typename T>
size_t VcfHeader::getIndex(const char* field, const std::vector<T>& idx,
                           const char* name) const
{
    typename std::vector<T>::const_iterator it = idx.begin();

    while (it != idx.end())
    {
        if (!strcmp(name, it->getKey()) )
        {
            return (size_t)(it - idx.begin());
        }

        ++it;
    }

    BOOST_THROW_EXCEPTION( illumina::common::OutOfBoundsException(
                               ( boost::format( "Unable to find %s's id '%s' in VCF header")
                                 % std::string(field) % std::string(name)).str() ));
    return 0L;
}


/*****************************************************************************/

std::istream& operator>>(std::istream& is, VcfHeader& vcfHeader)
{
    std::string line("###");
    vcfHeader.clear();
    while ( is.good() && (2 < line.size())
            && ('#' == line[0]) && ('#' == line[1]) )
    {
        do
        {
            std::getline(is,line);
        }
        while (is.good() && line.empty());

        if (!is.good() || 2 > line.size() || '#' != line[0]) break;

        if ('#' != line[1])
        {
            // Should be the #CHROM line.
            line = line.substr(1); // strip the leading hash.
            SplitString<'\t'> labels(line);
            const unsigned int numFixedFields(VcfFields::numFixedFields());
            const unsigned int numLabels(labels.size());

            if (numLabels < numFixedFields)
            {
                BOOST_THROW_EXCEPTION(VcfException((boost::format("Not enough fields in labels line '%s'") % line).str()));
            }

            for (unsigned int fieldInd(0); fieldInd < numFixedFields;
                 ++fieldInd)
            {
                if (labels[fieldInd]
                    != std::string(VcfFields::FIXED[fieldInd]))
                {
                    BOOST_THROW_EXCEPTION(VcfException((boost::format("In labels line found '%s' when expecting '%s'") % (labels[fieldInd],
                                                        VcfFields::FIXED[fieldInd])).str()));
                }
            }

            // If have sample(s), must be preceded by format but format
            // without sample(s) does not make sense.
            if (numLabels > numFixedFields)
            {
                if (labels[numFixedFields]
                    != std::string(VcfFields::GENOTYPE[0]))
                {
                    BOOST_THROW_EXCEPTION(VcfException((boost::format("In labels line found '%s' when expecting '%s'") % (labels[numFixedFields],
                                                        VcfFields::GENOTYPE[0])).str()));
                }

                if (numLabels == (numFixedFields + 1))
                {
                    BOOST_THROW_EXCEPTION(VcfException((boost::format("In labels line found '%s' but no sample names") % VcfFields::GENOTYPE[0]).str()));
                }

                for (unsigned int fieldInd(numFixedFields + 1);
                     fieldInd < numLabels; ++fieldInd)
                {
                    // DEBUG
                    // std::cerr << "Sample name `" << labels[fieldInd] << "'"
                    //           << std::endl;
                    vcfHeader.sampleNameList_.push_back(labels[fieldInd]);
                }
            }

            break;
        }

        size_t equalPos = line.find('=',2);
        if (std::string::npos == equalPos) continue;
        std::string key = line.substr(2,equalPos-2);
        std::string value = line.substr(equalPos+1);
        if      ("fileformat" == key)  vcfHeader.fileFormat_ = value;
        else if ("fileDate" == key)    vcfHeader.fileDate_ = value;
        else if ("source" == key)      vcfHeader.source_ = value;
        else if ("reference" == key)   vcfHeader.reference_ = value;
        else if ("phasing" == key)     vcfHeader.phasing_ = value;
        else if ("contig" == key)
        {
            if ('<' != value[0] || '>' != value[value.length()-1])
            {
                BOOST_THROW_EXCEPTION(VcfException( (boost::format("Unrecognized structure in header meta-information. "
                                                                   "Did not understand '##contig=%s'") % value).str() ));
            }
            ChromosomeMetadata chr;
            SplitString<','> buffer( std::string( value.begin()+1, value.begin()+value.length()-1 ));

            for ( unsigned int i = 0; i < buffer.size(); i++ )
            {
                SplitString<'='> field(buffer[i]);

                if (2 != field.size())
                {
                    BOOST_THROW_EXCEPTION(VcfException( (boost::format("Stray '=' in header meta-information. "
                                                                       "Did not understand '%s' in ##contig")
                                                         % field.toString()).str() ));
                }

                if ( "ID" == field[0] )
                {
                    chr.setId( field[1] );
                }
                else if ( "length" == field[0] )
                {
                    chr.setLength( boost::lexical_cast<size_t>(field[1]) );
                }
                else if ( "assembly" == field[0] )
                {
                    chr.setAssembly( field[1] );
                }
                else if ( "md5" == field[0] )
                {
                    chr.setMd5( field[1] );
                }
                else if ( "species" == field[0] )
                {
                    chr.setSpecies( field[1] );
                }
                else if ( "URL" == field[0] )
                {
                    chr.setUrl( field[1] );
                }
                else if ( "taxonomy" == field[0] )
                {
                    // ##contig=<...taxonomy=""...> not supported
                }
                else
                {
                    BOOST_THROW_EXCEPTION(VcfException( (boost::format("Unrecognized keyword in ##contig. "
                                                                       "Did not understand '%s'") % field[0]).str() ));
                }
            }
            vcfHeader.chrList_.push_back( chr );
        }
        else
        {
            std::istringstream buffer;
            buffer.str(value);
            VcfMetaInformation metaInfo;
            // TODO: make sure that the key is persistent!
            // TODO: check for doublons and other stuff
            if ( !strcmp(VcfFields::FIXED[7],key.c_str()) )
            {
                buffer >> metaInfo;
                vcfHeader.infoList_.push_back(  metaInfo );
                if ( "END" == metaInfo.getId() )
                {
                    // this dirty hack allows us to speed up parsing by not having to check each line for the END tag
                    vcfHeader.blockCompressed_ = true;
                    log_os << "VCF file is block-compressed" << std::endl;
                }
            }
            else if ( !strcmp(VcfFields::FIXED[6],key.c_str()) )
            {
                buffer >> metaInfo;
                vcfHeader.filterList_.push_back(  metaInfo );
            }
            else if ( !strcmp(VcfFields::GENOTYPE[0],key.c_str()) )
            {
                buffer >> metaInfo;
                vcfHeader.formatList_.push_back(  metaInfo );
            }
            else if ( !strcmp(VcfFields::FIXED[4],key.c_str()) )
            {
                buffer >> metaInfo;
                vcfHeader.altList_.push_back(  metaInfo );
            }
            else
            {
                if (vcfHeader.isStrict())
                {
                    log_os << "Unrecognized keyword in header meta-information. "
                           << (boost::format("Did not understand '%s': Ignoring.") % key).str() << std::endl;
                }
            }
        }
    }
    return is;
}

/*****************************************************************************/

std::ostream& operator<<(std::ostream& os, const VcfHeader& vcfHeader)
{
    std::string buffer = vcfHeader.getFileFormat();
    if (!buffer.empty())  os << "##fileformat=" << buffer << std::endl;
    buffer = vcfHeader.getFileDate();
    if (!buffer.empty())  os << "##fileDate=" << buffer << std::endl;
    buffer = vcfHeader.getSource();
    if (!buffer.empty())  os << "##source=" << buffer << std::endl;
    buffer = vcfHeader.getReference();
    if (!buffer.empty())  os << "##reference=" << buffer << std::endl;
    buffer = vcfHeader.getPhasing();
    if (!buffer.empty())  os << "##phasing=" << buffer << std::endl;

    BOOST_FOREACH(const ChromosomeMetadata &chr,
                  vcfHeader.getContigList())
    {
        os << "##contig=<ID=" << chr.getId() << ",length=" << chr.getLength()
           << ">" << std::endl;
    }
    BOOST_FOREACH(const VcfMetaInformation &info, vcfHeader.getInfoList())
    {
        os << "##INFO=" << info << std::endl;
    }
    BOOST_FOREACH(const VcfMetaInformation &filter,
                  vcfHeader.getFilterList())
    {
        os << "##FILTER=" << filter << std::endl;
    }
    BOOST_FOREACH(const VcfMetaInformation &format,
                  vcfHeader.getFormatList())
    {
        os << "##FORMAT=" << format << std::endl;
    }
    BOOST_FOREACH(const VcfMetaInformation &alt, vcfHeader.getAltList())
    {
        os << "##ALT=" << alt << std::endl;
    }
    return os;
}

/*****************************************************************************/
