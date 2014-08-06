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

#include "VcfLocus.hh"
#include "VcfHeader.hh"
#include "VcfMetaInformation.hh"
#include "ContigList.hh"

/*****************************************************************************/

/**
 * \brief format-less element in a generic file
 */
class FileElement
{
public:
    FileElement() : fileId_(0), repeat_(0) {}

    int getFileId() const
    {
        return fileId_;
    }
    void setFileId( int fileId)
    {
        fileId_ = fileId;
    }
    size_t getRepeat() const
    {
        return repeat_;
    }
    void setRepeat( size_t repeat)
    {
        repeat_ = repeat;
    }
    void incRepeat()
    {
        ++repeat_;
    }
    void decRepeat()
    {
        --repeat_;
    }

    virtual void init() = 0;
    virtual ~FileElement() {}

    void swap(FileElement& rhs)
    {
        std::swap( fileId_, rhs.fileId_ );
        std::swap( repeat_, rhs.repeat_ );
    }

    FileElement& operator=(const FileElement& rhs)
    {
        fileId_ = rhs.fileId_;
        repeat_ = rhs.repeat_;
        return *this;
    }
protected:
    int fileId_;
    size_t repeat_;
};

/*****************************************************************************/

class Sample
{
public:
    Sample() { };

    const char* field(size_t fieldInd) const
    {
        return fieldVec_.at(fieldInd);
    }

    std::vector<const char*> fieldVec_;
};

std::ostream& operator<<(std::ostream& ostrm, const Sample& sample);

/*****************************************************************************/

class VcfLine : public FileElement, public VcfLocus<unsigned int, size_t>
{
    typedef bool (VcfLine::*ReadField)(std::string::iterator& begin,
                                       const std::string::iterator& end);
    std::vector< ReadField > READ;

public:
    VcfLine& self()
    {
        return *this;
    }
    static const char* const MISSING;

    // for QUAL
    static const double MISSING_QUAL;
    // for FILTER
    static const char* const PASS;
    static const char* const FAIL;
    // for INFO flags
    static const char* const TRUE;

public:
    // default constructor
    VcfLine()
        : FileElement(), VcfLocus<unsigned int, size_t>(), vcfHeader_(0),
          ref_(0), qual_(0), pass_(false), spuriousHeader_(false)
    {
        ;
    }

    bool strictParsing() const
    {
        assert(vcfHeader_);
        return vcfHeader_->isStrict();
    }
    /// \return CHROM: An identifier (i.e., chromosome name) from the reference genome
    const char* getChrom() const
    {
        assert(vcfHeader_);
        return vcfHeader_->hasContigList()
               ? vcfHeader_->getContig( VcfLocus<unsigned int, size_t>::getChromosome() ).getKey()
               : ContigList::get_const_instance().getContig( VcfLocus<unsigned int, size_t>::getChromosome() ).getKey();
    }

    /// \return POS: The reference position, or MAX_POS if missing. (the 1st base having position 1)
    unsigned long getPos() const
    {
        return VcfLocus<unsigned int, size_t>::getPosition();
    }

//    /// \return ID: vector of pointers to unique identifiers where available
//    const void getId(std::vector<const char *> &id) const {id = id_;}
    /// \return ID: Semi-colon separated list of unique identifiers where available
    std::string getId() const
    {
        return cStringJoin(id_, ";");
    }

    /// \return REF: reference base(s)
    std::string getRef() const
    {
        return ref_;
    }

//    /// \return ALT: vector of pointers to alternate non-reference alleles called on at least one of the samples
//    const void getAlt(std::vector<const char *> &alt) const {alt = alt_;}
    /// \return ALT: comma separated list of alternate non-reference alleles called on at least one of the samples
    std::string getAlt() const
    {
        return cStringJoin(alt_, ",");
    }

    /// \return QUAL: phred-scaled quality score for the assertion made in ALT, or -1 if missing
    double getQual() const
    {
        char* endptr = NULL;
        errno = 0;
        return (MISSING != qual_) ? strtod(qual_, &endptr) : MISSING_QUAL;
    }

    /**
     ** \brief FILTER: whether a particular filter passes or not
     ** \return Either "." (missing), or "PASS", or "FAIL".
     **/
    const char* getFilter(size_t i) const
    {
        return  ( pass_ ? PASS    :
                  ( filter_.empty() ? MISSING :
                    ( filter_.at(i) ? FAIL : PASS )));
    }

    /**
     ** \brief FILTER
     ** \return Entire FILTER string
     **/
    std::string getFilter() const
    {
        assert(vcfHeader_);
        return  ( pass_ ? std::string(PASS,PASS+4)
                  : joinIfTrue( filter_, vcfHeader_->getFilterList(),";") );
    }

    /**
     ** \brief FILTER==PASS
     ** \return Whether the VCF line passes all filters or not
     **/
    bool pass() const
    {
        return pass_;
    }

    /**
     ** \brief INFO: additional information
     ** \return the text of the corresponding INFO field, if present; NULL otherwise. Note that for flags, the return value is VcfLine::TRUE.
     ** \throw out_of_range
     **/
    const char* getInfo(size_t i) const
    {
        return info_.empty() ? NULL : info_.at(i);
    }

    /**
     ** \brief FORMAT & SAMPLE
     ** \return the text of the corresponding SAMPLE, based on the FORMAT field, if present; NULL otherwise.
     ** \throw out_of_range
     **/
    const char* getSampleField(size_t sampleInd, size_t fieldInd) const
    {
        return format_.empty() ? NULL : sampleVec_.at(sampleInd).field(fieldInd);
    }

    /// \return INFO
    std::string getInfo() const
    {
        assert(vcfHeader_);
        return createMap(info_, vcfHeader_->getInfoList(), ";");
    }

    /// \return FORMAT
    std::string getFormat() const
    {
        assert(vcfHeader_);
        return joinByIndex(format_, vcfHeader_->getFormatList(), ":");
    }

    /// \return SAMPLE
    std::string getSample(size_t sampleInd) const
    {
        return cStringJoin(sampleVec_.at(sampleInd).fieldVec_, ":");
    }

    // load the raw data from the stream. NOTE: invalidates the VcfLine
    std::istream& read(std::istream& is);

    // copy raw data between lines
    VcfLine& swapRaw(VcfLine& rhs);
    VcfLine& copyRaw(const VcfLine& rhs);

    // parse the encoded_ data. returns true on success
    bool parse();
    void validate(std::vector<bool>& warn) const;

//    const std::string getSample(const std::string &format) const;
    //char * const getSample(const std::string &format) const;

    void init()
    {
        vcfString_.setLine( &(*this) );
        vcfBool_.setLine( &(*this) );
        vcfFullString_.setLine( &(*this) );
        vcfLong_.setLine( &(*this) );
        READ = { &VcfLine::readId,
               &VcfLine::readRef,
               &VcfLine::readAlt,
               &VcfLine::readQual,
               &VcfLine::readFilter};
        spuriousHeader_ = false;
    }

    void setHeader(const VcfHeader* vcfHeader)
    {
        assert(vcfHeader);
        vcfHeader_ = vcfHeader;
        vcfString_.setHeader( &(*vcfHeader_) );
        vcfBool_.setHeader( &(*vcfHeader_) );
        vcfBool_.setIndex( &VcfHeader::getFilterIndex );
        vcfFullString_.setHeader( &(*vcfHeader_) );
        vcfFullString_.setIndex( &VcfHeader::getInfoIndex );
        vcfLong_.setHeader( &(*vcfHeader_) );
    }
//    void swapHeader(VcfLine &vcfLine)
//    {
//        const VcfHeader *tmpHeader = vcfLine.vcfHeader_;
//        vcfLine.setHeader( &(*vcfHeader_) );
//        this->setHeader( &(*tmpHeader) );
//    }

private:
    const VcfHeader* vcfHeader_;
    Dsv< VcfHeader, VcfLine, const char*, '\t' >    vcfString_;
    Dsv< VcfHeader, VcfLine, bool,         '\t' >    vcfBool_;
    Dsv< VcfHeader, VcfLine, const char*, '\t' >    vcfFullString_;
    Dsv< VcfHeader, VcfLine, size_t,       '\t' >    vcfLong_;

    const char* direct(const char* arg) const
    {
        return arg;
    }
    bool tautology(const char* /*arg*/) const
    {
        return true;
    }
    const char* pseudoTautology(const char* /*arg*/) const
    {
        return TRUE;
    }
    size_t getFormatIndex(const char* arg) const
    {
        assert(vcfHeader_);
        return vcfHeader_->getFormatIndex(arg);
    }

    std::string unparsed_;

    /**
     ** \brief 3. ID semi-colon separated list of unique identifiers where available.
     ** If this is a dbSNP variant it is encouraged to use the rs number(s). No
     ** identifier should be present in more than one data record. If there is
     ** no identifier available, then the missing value should be used. (Alphanumeric String)
     */
    std::vector<const char*> id_;
    /**
     ** \brief 4. REF reference base(s): Each base must be one of A,C,G,T,N.
     ** Bases should be in uppercase. Multiple bases are permitted.
     ** The value in the POS field refers to the position of the first base in
     ** the String.
     ** For InDels, the reference String must include the base before the
     ** event (which must be reflected in the POS field). (String, Required).
     **/
    char* ref_;
    /**
     ** \brief 5. ALT comma separated list of alternate non-reference alleles called on at least one of the samples.
     ** Options are base Strings made up of the bases A,C,G,T,N, or an
     ** angle-bracketed ID String ("<ID>"). If there are no alternative
     ** alleles, then the missing value should be used. Bases should be in
     ** uppercase. (Alphanumeric String; no whitespace, commas, or
     ** angle-brackets are permitted in the ID String itself)
     **/
    std::vector<const char*> alt_;
    /**
     ** \brief replace all ID strings '<ID>' in the list by the corresponding sequence from the VCF header.
     **/
    bool replaceAltIds();
    /**
     ** \brief 6. QUAL phred-scaled quality score for the assertion made in ALT. i.e. give -10log_10 prob(call in ALT is wrong).
     ** If ALT is "." (no variant) then this is -10log_10 p(variant), and if
     ** ALT is not "." this is -10log_10 p(no variant). High QUAL scores
     ** indicate high confidence calls. Although traditionally people use
     ** integer phred scores, this field is permitted to be a floating point
     ** to enable higher resolution for low confidence calls if desired.
     ** (Numeric)
     **/
    char* qual_;
    // double qual_; // we'll leave QUAL as string and lazily convert to double on demand, for efficiency reasons
    /**
     ** \brief 7. FILTER filter:
     ** PASS if this position has passed all filters, i.e. a call is made at this position. Otherwise, if the site
     ** has not passed all filters, a semicolon-separated list of codes for filters that fail. e.g. "q10;s50" might
     ** indicate that at this site the quality is below 10 and the number of samples with data is below 50% of the
     ** total number of samples. "0" is reserved and should not be used as a filter String. If filters have not been
     ** applied, then this field should be set to the missing value (Alphanumeric String)
     **/
    std::vector<bool> filter_;
    bool pass_;
    /**
     ** \brief 8. INFO additional information: (Alphanumeric String)
     ** INFO fields are encoded as a semicolon-separated series of short keys with optional values in the format:
     **   \<key\>=\<data\>[,data].
     ** Arbitrary keys are permitted, although the following sub-fields are reserved (albeit optional):
     ** \li AA ancestral allele
     ** \li AC allele count in genotypes, for each ALT allele, in the same order as listed
     ** \li AF allele frequency for each ALT allele in the same order as listed: use this when estimated from primary data, not called genotypes
     ** \li AN total number of alleles in called genotypes
     ** \li BQ RMS base quality at this position
     ** \li CIGAR cigar string describing how to align an alternate allele to the reference allele
     ** \li DB dbSNP membership
     ** \li DP combined depth across samples, e.g. DP=154
     ** \li END end position of the variant described in this record (esp. for CNVs)
     ** \li H2 membership in hapmap2
     ** \li MQ RMS mapping quality, e.g. MQ=52
     ** \li MQ0 Number of MAPQ == 0 reads covering this record
     ** \li NS Number of samples with data
     ** \li SB strand bias at this position
     ** \li SOMATIC indicates that the record is a somatic mutation, for cancer genomics
     ** \li VALIDATED validated by follow-up experiment
     ** etc. The exact format of each INFO sub-field should be specified in the meta-information (as described above).
     ** Example for an INFO field: DP=154;MQ=52;H2. Keys without corresponding values are allowed in order to indicate group membership (e.g. H2
     ** indicates the SNP is found in HapMap 2). It is not necessary to list all the properties that a site does NOT have, by e.g. H2=0.
     **
     ** Implementation note: for FLAGS, 'true' is the pointer "VcfLine::TRUE", 'false' is NULL.
     **/
    std::vector<const char*> info_;
    /**
     ** \brief FORMAT
     ** If genotype information is present, then the same types of data must be present for all samples. First a FORMAT field is given specifying
     ** the data types and order. This is followed by one field per sample, with the colon-separated data in this field corresponding to the types
     ** specified in the format. The first sub-field must always be the genotype (GT).
     **
     ** As with the INFO field, there are several common, reserved keywords that are standards across the community:
     ** \li GT genotype, encoded as alleles values separated by either of "/" or "|", e.g. The allele values are 0 for the reference allele (what is in the reference sequence), 1 for the first allele listed in ALT, 2 for the second allele list in ALT and so on. For diploid calls examples could be 0/1 or 1|0 etc. For haploid calls, e.g. on Y, male X, mitochondrion, only one allele value should be given. All samples must have GT call information; if a call cannot be made for a sample at a given locus, "." must be specified for each missing allele in the GT field (for example ./. for a diploid). The meanings of the separators are: "/ : genotype unphased" and "| : genotype phased"
     ** \li DP read depth at this position for this sample (Integer)
     ** \li FT sample genotype filter indicating if this genotype was "called" (similar in concept to the FILTER field). Again, use PASS to indicate that all filters have been passed, a semi-colon separated list of codes for filters that fail, or "." to indicate that filters have not been applied. These values should be described in the meta-information in the same way as FILTERs (Alphanumeric String)
     ** \li GL : three floating point log10-scaled likelihoods for AA,AB,BB genotypes where A=ref and B=alt; not applicable if site is not biallelic. For example: GT:GL 0/1:-323.03,-99.29,-802.53 (Numeric)
     ** \li GQ genotype quality, encoded as a phred quality -10log_10p(genotype call is wrong) (Numeric)
     ** \li HQ haplotype qualities, two phred qualities comma separated (Numeric)
     **
     ** If any of the fields is missing, it is replaced with the missing value. For example if the format is GT:GQ:DP:HQ then A|A:.:23:23,34 indicates
     ** that GQ is missing. Trailing fields can be dropped (with the exception of the GT field, which should always be present).
     **
     ** Additional Genotype fields can be defined in the meta-information. However, software support for such fields is not guaranteed.
     **
     ** Implementation note: the number of elements in the format_ is exactly the number of fields in the FORMAT of the current line. Each element is the index of the corresponding field in the header.
     **/
    std::vector<size_t> format_;
    /**
     ** \brief SAMPLE
     **
     ** Implementation note: the number of elements in the sample_ is exactly the number of FORMAT fields in the header. Fields that are not present in the current line are NULL.
     **/
    std::vector<Sample> sampleVec_;

    bool spuriousHeader_;

    /// skips empty lines and comments and reads the chromosome
    std::istream& readChrom(std::istream& is);
    /// read the position and convert it to an integer
    std::istream& readPos(std::istream& is);

    bool readId(std::string::iterator& begin,
                const std::string::iterator& end);
    bool readRef(std::string::iterator& begin,
                 const std::string::iterator& end);
    bool readAlt(std::string::iterator& begin,
                 const std::string::iterator& end);
    bool readQual(std::string::iterator& begin,
                  const std::string::iterator& end);
    bool readFilter(std::string::iterator& begin,
                    const std::string::iterator& end);
    // last mandatory field, can end with a '\n' or '\t'
    bool readInfo(std::string::iterator& begin,
                  const std::string::iterator& end);
    //std::istream &readInfo(std::string::const_iterator &begin, const std::string::const_iterator end);
    // first optional column
    bool readFormat(std::string::iterator& begin,
                    const std::string::iterator& end);
    bool readSample(unsigned int sampleInd,
                    std::string::iterator& begin,
                    const std::string::iterator& end);
    bool readSamples(std::string::iterator& begin,
                     const std::string::iterator& end);

public:
    // TODO: these should actually belong to Dsv
    // boost::join doesn't support vectors of 'char *'
    static std::string cStringJoin(const std::vector<const char*>& v, const char* separator);

private:
    std::string joinByIndex(const std::vector<size_t>& idx,
                            const std::vector<VcfMetaInformation>& v,
                            const char* separator) const;
    std::string joinIfTrue( const std::vector<bool>& idx,
                            const std::vector<VcfMetaInformation>& v,
                            const char* separator) const;
    std::string createMap(  const std::vector<const char*>& info,
                            const std::vector<VcfMetaInformation>& v,
                            const char* separator) const;

public:
    bool operator<( const VcfLine& rhs ) const
    {
        return VcfLocus<unsigned int, size_t>::operator<( dynamic_cast< const VcfLocus<unsigned int, size_t>& >(rhs) );
    }
    friend std::istream& operator>>(std::istream& is, VcfLine& vcfLine);
    friend std::ostream& operator<<(std::ostream& os, const VcfLine& vcfLine);
};

/*****************************************************************************/

std::istream& operator>>(std::istream& is, VcfLine& vcfLine);
std::ostream& operator<<(std::ostream& os, const VcfLine& vcfLine);

/// specialize DualStream for VCF use case
//typedef PairedStream<boost::iostreams::file_source,VcfLine> AnnotatedVcfLine;
// typedef AutoStream<VcfLine> AnnotatedVcfLine;

/// specialize ArrayStream for VCF use case
// typedef ArrayStream<AnnotatedVcfLine> VcfStreams;

/*****************************************************************************/
