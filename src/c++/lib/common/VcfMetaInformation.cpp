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
/// \author Refactored by Richard Shaw from PUMA by Mauricio Varea
///

/*****************************************************************************/

#include <boost/format.hpp>
#include <boost/throw_exception.hpp>

#include "common/Exceptions.hh"
using namespace illumina::common;

#include "VcfMetaInformation.hh"

/*****************************************************************************/

struct VcfTransform
{
    static const char* type2string[];
    static const char* XmlFields[];

    static inline void inc(VCF_VALUE_TYPE::index_t& m)
    {
        m = static_cast<VCF_VALUE_TYPE::index_t>(1 + (int)m);
    }

    static inline VCF_VALUE_TYPE::index_t string2type(const std::string& s)
    {
        using namespace VCF_VALUE_TYPE;

        index_t r = NONE;
        while (r != SIZE)
        {
            if ( s == std::string( type2string[r] ) )
            {
                return r;
            }
            inc(r);
        }
        BOOST_THROW_EXCEPTION(PostConditionException(
                                  (boost::format("Inexistent VcfMetainformation type: %s") % s).str() ));
    }
};

/*****************************************************************************/

const char* VcfTransform::type2string[]
    = {"UNDEF","Integer","Float","Flag","Character","String","UNDEF"};

const char* VcfTransform::XmlFields[]
    = {"ID","Number","Type","Description"};


/*****************************************************************************/
// VcfMetaInformation
/*****************************************************************************/

VcfMetaInformation::VcfMetaInformation()
    : type_(VCF_VALUE_TYPE::NONE)
{
    ;
}

/*****************************************************************************/

void VcfMetaInformation::validateFlag() const
{
    if ( VCF_VALUE_TYPE::FLAG == type_ && "0" != number_ )
    {
        BOOST_THROW_EXCEPTION(VcfException( (boost::format("The 'Flag' type indicates that the INFO field does not contain a Value entry, "
                                                           "and hence the Number should be 0 in this case. "
                                                           "Please check the header meta-information (found a 'Flag' with '%s' in the 'Number' field)")
                                             % number_).str() ));
    }
}

/*****************************************************************************/

void VcfMetaInformation::validateNumber(const std::string& num) const
{
    static const std::string validNumbers("0123456789.AG");
    if (VCF_VALUE_TYPE::NONE != type_)
    {
        size_t found = num.find_first_not_of( validNumbers );
        if ( std::string::npos != found )
        {
            BOOST_THROW_EXCEPTION(VcfException( (boost::format("The 'Number' entry is a pseudo-integer that describes "
                                                               "the number of values that can be included with the INFO field. "
                                                               "Valid values are \"%s\", found '%c'")
                                                 % validNumbers % num[found]).str() ));
        }
    }
}

/*****************************************************************************/

void VcfMetaInformation::validateDescription(const std::string& desc) const
{
    if ( '"' != *(desc.begin()) || '"' != *(desc.rbegin()) )
    {
        BOOST_THROW_EXCEPTION(VcfException( (boost::format("The Description value must be surrounded by double-quotes. "
                                                           "Please check the header meta-information (Description=%s)")
                                             % desc).str() ));
    }
}

/*****************************************************************************/

std::istream& operator>>(std::istream& is, VcfMetaInformation& info)
{
    using namespace VCF_VALUE_TYPE;

//    const unsigned int bufferSize = sizeof(VcfTransform::XmlFields) / sizeof(VcfTransform::XmlFields[0]);
    SplitString<','> buffer;
    if ('<' == is.peek())
    {
        is.ignore();
        std::getline( is, buffer.ref(), '>' );

        if (2 > buffer.size())
        {
            BOOST_THROW_EXCEPTION(VcfException( (boost::format("Unrecognized structure in header meta-information. "
                                                               "Did not understand '<%s>'") % buffer.toString()).str() ));
        }
        SplitString<'='> id(buffer[0]);
        info.setId( id[1] );
        if ( 4 > buffer.size())
        {
            buffer.mergeAfter( 2 );
            info.setNumber( "." );
            info.setType( NONE );
            SplitString<'='> desc(buffer[1]);
            desc.mergeAfter( 2 );
            info.setDescription( desc[1] );
        }
        else
        {
            buffer.mergeAfter( 4 );
            SplitString<'='> type( buffer[2] );
            if ( 2 != type.size() || strcmp(VcfTransform::XmlFields[2], type[0].c_str()) )
            {
                info.setNumber( "." );
                info.setType( NONE );
                SplitString<'='> desc(buffer[1]);
                info.setDescription( desc[1] + buffer[2] + buffer[3] );
            }
            else
            {
                info.setType( VcfTransform::string2type( type[1] ));
                SplitString<'='> number(buffer[1]);
                if ( strcmp(VcfTransform::XmlFields[1], number[0].c_str()) )
                {
                    BOOST_THROW_EXCEPTION(VcfException( (boost::format("Unrecognized assignment in '%s'. "
                                                                       "Expecting '%s' as a key.")
                                                         % number.toString()
                                                         % VcfTransform::XmlFields[1]).str() ));
                }
                info.setNumber( number[1] );
                SplitString<'='> desc(buffer[3]);
                desc.mergeAfter( 2 );
                if ( strcmp(VcfTransform::XmlFields[3], desc[0].c_str()) )
                {
                    BOOST_THROW_EXCEPTION(VcfException( (boost::format("Unrecognized assignment in '%s'. "
                                                                       "Expecting '%s' as a key.")
                                                         % desc.toString()
                                                         % VcfTransform::XmlFields[3]).str() ));
                }
                info.setDescription( desc[1] );
            }
        }
    }
    info.validateFlag();
    return is;
}

/*****************************************************************************/

std::ostream& operator<<(std::ostream& os, const VcfMetaInformation& info)
{
    using namespace VCF_VALUE_TYPE;
    os << "<" << VcfTransform::XmlFields[0] << "=" << info.getId();

    if ( NONE != info.getType() )
    {
        os << "," << VcfTransform::XmlFields[1] << "=" << info.getNumber()
           << "," << VcfTransform::XmlFields[2]
           << "=" << VcfTransform::type2string[ info.getType() ];
    }

    os << "," << VcfTransform::XmlFields[3]
       << "=" << info.getDescription() << ">";
    return os;
}


/*****************************************************************************/
