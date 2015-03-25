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

#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "blt_util/log.hh"

#include "VcfFields.hh"

#include "VcfLine.hh"

using namespace illumina::common;


/*****************************************************************************/
// VcfLine helper fns
/*****************************************************************************/

/**
 * \brief Validates bases in VCF file, especially REF and ALT cols
 *
 * We use this conversion table rather than the one in Nucleotide.hh to make sure we conform to the VCF standard
 */
static
std::vector<bool> getValidReferenceBases()
{
    std::vector<bool> validReferenceBases(256, false);
    validReferenceBases['A'] = true;
    validReferenceBases['C'] = true;
    validReferenceBases['G'] = true;
    validReferenceBases['T'] = true;
    validReferenceBases['N'] = true;
    validReferenceBases['a'] = true;
    validReferenceBases['c'] = true;
    validReferenceBases['g'] = true;
    validReferenceBases['t'] = true;
    validReferenceBases['n'] = true;
    return validReferenceBases;
}

/*****************************************************************************/

static
bool validateRef(const char* ref)
{
    static const std::vector<bool> validReferenceBases = getValidReferenceBases();
    while (*ref)
    {
        // if(!genome::Nucleotide::valid( tmp )
        if (!(validReferenceBases[*ref] || (*ref == '.')))
        {
            return false;
        }
        ++ref;
    }
    return true;
}

/*****************************************************************************/
static
std::vector<bool> getValidFloatingPointDigits()
{
    std::vector<bool> validFloat(256, false);
    for ( char c = '0'; c <= '9'; c++ )
    {
        validFloat[c] = true;
    }
    validFloat['.'] = true;
    validFloat['e'] = true;
    validFloat['E'] = true;
    return validFloat;
}

/*****************************************************************************/
static
bool validateQual(const char* ref)
{
    static const std::vector<bool> validFloatingPointDigits
        = getValidFloatingPointDigits();

    while (*ref)
    {
        if (!validFloatingPointDigits[*ref])
        {
            return false;
        }
        ++ref;
    }

    return true;
}

/*****************************************************************************/
static
char null2tab(char c)
{
    return ('\0'==c) ? '\t' : c;
}

/*****************************************************************************/
static
void prettyPrint(const std::string& line)
{
    std::string prettyLine;
    prettyLine.resize(line.length(),' ');
    std::transform( line.begin(), line.end(), prettyLine.begin(), null2tab);
    log_os << "*** raw data:   [ " << prettyLine << " ]" << std::endl;
}

/*****************************************************************************/
namespace
{
template< int N >
void fixedException(VcfLine* vcf, const std::string& field,
                    const std::string& line)
{
    std::string message = (boost::format( "Undefined %s key '%s'\t(%s:%i)")
                           % std::string( VcfFields::FIXED[N] )
                           % field
                           % vcf->getChromosome()
                           % vcf->getPosition() ).str();
    if (vcf->strictParsing())
    {
        prettyPrint(line);
        BOOST_THROW_EXCEPTION(VcfException( message ));
    }
    else
    {
        static bool isWarn(true);
        if (isWarn)
        {
            isWarn = false;
            log_os << "WARNING: " << message << "\n";
        }
    }
}

/*****************************************************************************/

template< int N >
void genotypeException(VcfLine* vcf, const std::string& field,
                       const std::string& line)
{
    std::string message = (boost::format( "Undefined %s key '%s'\t(%s:%i)")
                           % std::string( VcfFields::GENOTYPE[N] )
                           % field
                           % vcf->getChromosome()
                           % vcf->getPosition() ).str();
    if (vcf->strictParsing())
    {
        prettyPrint(line);
        BOOST_THROW_EXCEPTION(VcfException( message ));
    }
    else
    {
        static bool isWarn(true);
        if (isWarn)
        {
            isWarn = false;
            log_os << "WARNING: " << message << "\n";
        }
    }
}

} // anom namespace

/*****************************************************************************/

std::ostream& operator<<(std::ostream& ostrm, const Sample& sample)
{
    ostrm << VcfLine::cStringJoin(sample.fieldVec_, ":");

    return ostrm;
}

/*****************************************************************************/
// VcfLine
/*****************************************************************************/

const char* const VcfLine::MISSING = ".";

/// QUAL is (phread-based) logarithmic, so we use -inf to represent the missing value
const double VcfLine::MISSING_QUAL = 0.0;

const char* const VcfLine::PASS = "PASS\t";
const char* const VcfLine::FAIL = "FAIL\t";

const char* const VcfLine::TRUE = "TRUE";

/*****************************************************************************/

std::istream& VcfLine::read(std::istream& is)
{
    if (readChrom(is))
    {
        if (!readPos(is))
        {
            prettyPrint(unparsed_);
            BOOST_THROW_EXCEPTION(VcfException( (boost::format("Failed to parse  (read) VCF position field:%s:%i")
                                                 % std::string( this->getChrom() ) % position_).str() ));
        }
        unparsed_.clear();
        if (!getline(is, unparsed_))
        {
            prettyPrint(unparsed_);
            BOOST_THROW_EXCEPTION(VcfException( (boost::format("Failed to read VCF line:%s:%i")
                                                 % std::string( this->getChrom() ) % position_ ).str() ));
        }
        // for consistency in the end of each field
        unparsed_.push_back('\t');
        id_.clear();
        ref_ = NULL;
        alt_.clear();
        filter_.clear();
        info_.clear();
        format_.clear();
        sampleVec_.clear();
    }
    return is;
}

/*****************************************************************************/

VcfLine& VcfLine::swapRaw(VcfLine& rhs)
{
    // FileElement::operator=(dynamic_cast<const FileElement&>(rhs));
    VcfLocus<unsigned int, size_t>::operator=( dynamic_cast< const VcfLocus<unsigned int, size_t>& >(rhs) );
    this->unparsed_.swap(rhs.unparsed_);
    this->setHeader( rhs.vcfHeader_ );            // no need to swap, as all lines should point to the same header
    this->spuriousHeader_ = rhs.spuriousHeader_;  // no need to swap, due to OR logic
    return *this;
}

/*****************************************************************************/

VcfLine& VcfLine::copyRaw(const VcfLine& rhs)
{
    // FileElement::operator=(dynamic_cast<const FileElement&>(rhs));
    VcfLocus<unsigned int, size_t>::operator=( dynamic_cast< const VcfLocus<unsigned int, size_t>& >(rhs) );
    this->unparsed_.resize( rhs.unparsed_.size() );
    std::copy( rhs.unparsed_.begin(), rhs.unparsed_.end(), this->unparsed_.begin() );
    this->setHeader( rhs.vcfHeader_ );
    this->spuriousHeader_ = rhs.spuriousHeader_;
    return *this;
}

/*****************************************************************************/

bool VcfLine::parse()
{
    std::string::iterator begin = unparsed_.begin();
    const std::string::iterator end = unparsed_.end() - 1;

    for ( size_t i=0; i < READ.size(); i++ )
    {
        if (!(this->*READ[i])(begin, end))
        {
            prettyPrint(unparsed_);
            BOOST_THROW_EXCEPTION(VcfException( (boost::format("Failed to parse mandatory %s field: (%s:%i)")
                                                 % VcfFields::FIXED[ i+2 ] % std::string( this->getChrom() ) % position_).str() ));
        }
    }
    if ( readInfo(begin,end) )
    {
        assert(vcfHeader_);
        if (readFormat(begin, end) && !readSamples(begin, end))
        {
            prettyPrint(unparsed_);
            BOOST_THROW_EXCEPTION(VcfException( (boost::format( "Failed to parse optional VCF fields: (%s:%i)")
                                                 % std::string( this->getChrom() ) % position_ ).str() ));
        }
        if (vcfHeader_->getFormatCount() != sampleVec_[0].fieldVec_.size())
        {
            prettyPrint(unparsed_);
            BOOST_THROW_EXCEPTION(VcfException( (boost::format( "Number of samples inconsistent with number of formats:%i (expected %i): (%s:%i)")
                                                 % sampleVec_[0].fieldVec_.size() % vcfHeader_->getFormatCount() % std::string( this->getChrom() ) % position_ ).str() ));
        }
        if (end != begin)
        {
            prettyPrint(unparsed_);
            BOOST_THROW_EXCEPTION(VcfException( (boost::format( "Spurious VCF fields: Found more than 10 fields in (%s:%i). Multi-sample not implemented yet.")
                                                 % std::string( this->getChrom() ) % position_ ).str() ));
        }

        if ( vcfHeader_->isBlockCompressed() )
        {
            const char* endTag = this->getInfo( vcfHeader_->getInfoIndex("END") );
            if (NULL != endTag)
            {
                unsigned long endPos = static_cast<unsigned long>( atoi(endTag) );
                if ( this->position_ < endPos )
                {
                    this->setRepeat( endPos - this->position_);
                }
                else
                {
                    BOOST_THROW_EXCEPTION(VcfException( (boost::format( "Wrong END value: "
                                                                        "Expected a value greater than POS (%lu), got %lu")
                                                         % this->position_ % endPos ).str() ));
                }
            }
        }
    }
    return true;
}

/*****************************************************************************/

void VcfLine::validate(std::vector<bool>& warn) const
{
    warn.at(0) = warn.at(0) || spuriousHeader_;
}

/*****************************************************************************/
// used by VcfLine::readChrom

static std::istream& skipCommentsAndEmptyLines(std::istream& is, char& c,
                                               bool& sH)
{
    sH = false;
    while (is.get(c)) // exit on break when a non empty line is found
    {
        if ('#' == c)
        {
            sH = is.get(c) && ('#' == c);
            // skip until the end of the line
            while (is.get(c) && ('\n' != c))
            {
                // just skip
            }
        }
        if ('\n' != c)
        {
            break;
        }
    }
    return is;
}

/*****************************************************************************/

std::istream& VcfLine::readChrom(std::istream& is)
{
    char c;
    std::vector< char > chrom(0);
    if (skipCommentsAndEmptyLines(is, c, spuriousHeader_))
    {
        while (is && '\t' != c)
        {
            chrom.push_back(c);
            is.get(c);
        }
        assert(is && ('\t' == c));
        chrom.push_back('\0');
//        size_t len = chrom.size()-1;
        assert(vcfHeader_);

        // DEBUG
        // std::cerr << "vcfHeader_->hasContigList() : "
        //           << vcfHeader_->hasContigList() << std::endl;

        chromosome_ = vcfHeader_->hasContigList()
                      ? vcfHeader_->getContigIndex( &chrom.front() )
                      : ContigList::get_mutable_instance().getIndex( &chrom.front() );
    }
    return is;
}

/*****************************************************************************/

std::istream& VcfLine::readPos(std::istream& is)
{
    char c;
    if (is.get(c) && ('.' == c))
    {
        position_ = UndefinedPos();
        is.get(c);
    }
    else
    {
        position_ = 0;
        while (is && ('\t' != c))
        {
            assert (('0' <= c) && ('9' >= c));
            position_ *= 10;
            position_ += (c - '0');
            is.get(c);
        }
    }
    assert(is && ('\t' == c));
    return is;
}

/*****************************************************************************/

bool VcfLine::readId(std::string::iterator& begin,
                     const std::string::iterator& end)
{
    // DEBUG
    // std::cerr << "readId" << std::endl;

    return vcfString_.readList( begin, end, id_, &VcfLine::direct, ';',
                                &fixedException< 2 >,
                                false );
}

/*****************************************************************************/

bool VcfLine::readRef(std::string::iterator& begin,
                      const std::string::iterator& end)
{
    ref_ = vcfString_.read(begin, end);
    return ref_ && validateRef(ref_);
}

/*****************************************************************************/
// may be used in readAlt

bool VcfLine::replaceAltIds()
{
    for (size_t i = 0; alt_.size(); ++i)
    {
        if (!validateRef(alt_[i]))
        {
            // TODO implement proper tests instead of asserts
            assert('<' == alt_[i][0]);
            assert('>' == alt_[i][strlen(alt_[i]) - 1]);
//            alt_[i][strlen(alt_[i]) - 1] = 0;
            //alt_[i][0] = 0;
//            alt_[i] = vcfHeader_->getAlt(alt_[i] + 1).getKey();
        }
    }
    return true;
}

/*****************************************************************************/

bool VcfLine::readAlt(std::string::iterator& begin,
                      const std::string::iterator& end)
{
    // DEBUG
    // std::cerr << "readAlt" << std::endl;

    return vcfString_.readList( begin, end, alt_, &VcfLine::direct, ',',
                                &fixedException< 4 >,
                                false );  // && replaceAltIds();
}

/*****************************************************************************/

bool VcfLine::readQual(std::string::iterator& begin,
                       const std::string::iterator& end)
{
    qual_ = vcfString_.read(begin, end);
    return qual_ && validateQual(qual_);
}

/*****************************************************************************/

bool VcfLine::readFilter(std::string::iterator& begin,
                         const std::string::iterator& end)
{
    pass_ = (end != begin) && !strncmp(PASS,&(*begin),5);
    if ( pass_ )
    {
        filter_.clear();
        *(begin+4) = '\0';
        begin += 5;
    }
    else
    {
        assert(vcfHeader_);

        // DEBUG
        // std::cerr << "readFilter" << std::endl;

        return vcfBool_.readList(begin, end, filter_, &VcfLine::tautology, ';',
                                 &fixedException< 6 >,
                                 true,
                                 vcfHeader_->getFilterCount() );
    }
    return true;
}

/*****************************************************************************/

bool VcfLine::readInfo(std::string::iterator& begin,
                       const std::string::iterator& end)
{
    assert(vcfHeader_);

    // DEBUG
    // std::cerr << "readInfo : vcfHeader_->getInfoCount() "
    //           << vcfHeader_->getInfoCount() << std::endl;

    return vcfFullString_.readList( begin, end, info_,
                                    &VcfLine::pseudoTautology, ';',
                                    &fixedException< 7 >,
                                    true,
                                    vcfHeader_->getInfoCount(),
                                    '=' );
}


/*****************************************************************************/

bool VcfLine::readFormat(std::string::iterator& begin,
                         const std::string::iterator& end)
{
    assert(vcfHeader_);

    return vcfLong_.readList( begin, end, format_,
                              &VcfLine::getFormatIndex, ':',
                              &genotypeException< 0 >,
                              false,
                              vcfHeader_->getFormatCount() );
}

/*****************************************************************************/

bool VcfLine::readSample(unsigned int sampleInd,
                         std::string::iterator& begin,
                         const std::string::iterator& end)
{
    std::vector<const char*>& fieldVec(sampleVec_.at(sampleInd).fieldVec_);
    fieldVec.clear();

    if (end == begin)
    {
        return false;
    }

    const std::string::iterator origin = begin;
    begin = std::find(begin, end, vcfString_.delimiter());
    // TODO: implement proper check for the end
    *begin = '\0';

    if (format_.empty())
    {
        // TODO: properly check that empty or missing value
        assert((origin == begin) || ((origin + 1 == begin) && (*MISSING == *origin)));
    }
    else
    {
        assert(vcfHeader_);
        fieldVec.resize(vcfHeader_->getFormatCount(), 0);
        std::string::iterator sampleBegin = origin;
        BOOST_FOREACH(const size_t i, format_)
        {
            const char delimiter1(':');
            fieldVec.at(i) = vcfString_.read(sampleBegin, begin, delimiter1);

            // DEBUG
            // std::cerr << "readSample : fieldVec.at(" << i << ") : ["
            //           << fieldVec.at(i) << "]" << std::endl;

            // TODO: implement proper check for the end
            assert(begin >= sampleBegin);
            //++sampleBegin;
        }
        // TODO: implement proper check for the end
        assert(begin == sampleBegin);
    }

    return true;
}

/*****************************************************************************/

bool VcfLine::readSamples(std::string::iterator& begin,
                          const std::string::iterator& end)
{
    const unsigned int numSamples(vcfHeader_->numSamples());
    sampleVec_.resize(numSamples);

    for (unsigned int sampleInd(0); sampleInd < numSamples; ++sampleInd)
    {
        // DEBUG
        // std::cerr << "readSamples : sampleInd " << sampleInd << std::endl;

        readSample(sampleInd, begin, end);

        if (sampleInd < (numSamples - 1))
        {
            ++begin;
        }
    }

    assert(end == begin);

    return true;
}

/*****************************************************************************/

std::string VcfLine::cStringJoin(const std::vector<const char*>& v,
                                 const char* separator)
{
    if (v.empty())
    {
        return std::string(MISSING);
    }
    std::vector<const char*>::const_iterator current = v.begin();
    std::string result = std::string( *current );
    while (v.end() != ++current)
    {
        result += (std::string(separator) + std::string(*current));
    }
    return result;
}

/*****************************************************************************/

std::string VcfLine::joinByIndex(const std::vector<size_t>& idx,
                                 const std::vector<VcfMetaInformation>& v,
                                 const char* separator) const
{
    if (idx.empty())
    {
        return std::string(MISSING);
    }
    std::vector<size_t>::const_iterator current = idx.begin();
    std::string result = v[ *current ].getId();
    while (idx.end() != ++current)
    {
        result += (std::string(separator) + v[ *current ].getId() );
    }
    return result;
}

/*****************************************************************************/

std::string VcfLine::joinIfTrue(const std::vector<bool>& idx,
                                const std::vector<VcfMetaInformation>& v,
                                const char* separator) const
{
    if (idx.empty())
    {
        return std::string(MISSING);
    }
    // TODO: throw proper exception
    assert( idx.size() == v.size() );
//    std::vector<bool>::const_iterator current = idx.begin();
    std::string result = idx[0] ? v[ 0 ].getId() : "";
    for (unsigned int i=1; idx.size() > i; i++)
    {
        if (idx[i])
        {
            result += (std::string(separator) + v[ i ].getId() );
        }
    }
    return result;
}

/*****************************************************************************/

std::string VcfLine::createMap(const std::vector<const char*>& info,
                               const std::vector<VcfMetaInformation>& v,
                               const char* separator) const
{
    if (info.empty())
    {
        return std::string(MISSING);
    }

    assert( info.size() == v.size() );
    std::string result;
    bool first(true);

    for (unsigned int i=0; info.size() > i; i++)
    {
        if (!info[i]) continue;

        if (first)
        {
            first = false;
        }
        else
        {
            result += std::string(separator);
        }

        result += v[ i ].getId();

        if (TRUE != info[i])
        {
            result += ("=" + std::string(info[i]) );
        }
    }

    return result;
}

/*****************************************************************************/

std::istream& operator>>(std::istream& is, VcfLine& vcfLine)
{
    if (vcfLine.read(is) && !vcfLine.parse())
    {
        is.setstate(std::ios::badbit);
    }
    return is;
}

/*****************************************************************************/

std::ostream& operator<<(std::ostream& os, const VcfLine& vcfLine)
{
    os << vcfLine.getChrom() << "\t"
       << vcfLine.getPos() << "\t"
       << vcfLine.getId() << "\t"
       << vcfLine.getRef() << "\t"
       << vcfLine.getAlt() << "\t"
       << vcfLine.getQual() << "\t"
       << vcfLine.getFilter() << "\t"
       << vcfLine.getInfo() << "\t"
       << vcfLine.getFormat();

    BOOST_FOREACH(const Sample& sample, vcfLine.sampleVec_)
    {
        os << "\t" << sample;
    }

    return os;
}

/*****************************************************************************/

