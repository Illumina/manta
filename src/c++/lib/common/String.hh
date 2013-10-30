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

#pragma once

#include <string>
#include <vector>
#include <cassert>

#include "common/Exceptions.hh"

#include <boost/algorithm/string.hpp>
// #include <boost/range.hpp>
// #include <boost/foreach.hpp>

// #include <algorithm>


/*****************************************************************************/

/**
 * \brief Template-based (optionally) lazily-evaluated splitting of strings
 */
template <char D, bool LAZY=false>
class SplitString
{
public:
    typedef  std::string  type;
public:
    SplitString()
        : str_( "" )
        , sizeAlreadyEvaluated_( false )
        , size_(0)
        , tokenPos_( 1, std::string::npos )
        , tokens_(0)
    {}
    SplitString( const std::string& str )
        : str_( str )
        , sizeAlreadyEvaluated_( false )
        , size_(1)
        , tokenPos_( 2, std::string::npos )
    {
        if (! LAZY )  nonLazySplit();
    }

    SplitString& operator=( const std::string& str )
    {
        str_ = str;
        sizeAlreadyEvaluated_ = false;
        size_ = 1;
        if (! LAZY )  nonLazySplit();
        return *this;
    }

    size_t size()
    {
        if ( !sizeAlreadyEvaluated_ )
        {
            if ( LAZY )
            {
                if ( str_.empty() )
                {
                    return (size_t)0;
                }
                size_ = std::count( str_.begin(), str_.end(), D ) + 1;
                tokens_.resize(   size_,   std::string("") );
                tokenPos_.resize( size_+1, std::string::npos );
                tokenPos_[0] = 0;
                tokenPos_[size_] = str_.length()+1;
            }
            else
            {
                if ( tokens_.empty() )
                {
                    nonLazySplit();
                }
                size_ = tokens_.size();
            }
            sizeAlreadyEvaluated_ = true;
        }
        return size_;
    }

    bool empty()
    {
        return 0 == this->size();
    }

    void mergeAfter(size_t i)
    {
        if ( this->size() > i )
        {
            if ( LAZY )
            {
                tokenPos_.erase( tokenPos_.begin() + i, tokenPos_.begin() + size_ );
            }
            else
            {
                while ( size_ > i )
                {
                    tokens_[size_-2] += std::string(1,D) + tokens_[size_-1];
                    tokens_.pop_back();
                    --size_;
                }
            }
        }
    }

    std::string operator[]( const size_t index )
    {
        if ( this->empty() )  return std::string("");
        if ( LAZY && tokens_.at(index).empty() )
        {
            size_t nextTokenPos = getTokenPos(index+1);
            size_t thisTokenPos = getTokenPos(index);
            tokens_[index] = str_.substr( thisTokenPos, nextTokenPos - thisTokenPos - 1 );
        }
        return tokens_[index];
    }

    std::string& ref()
    {
        return str_;
    }
    const std::string& toString() const
    {
        return str_;
    }
private:
    size_t getTokenPos( size_t index ) const
    {
        size_t thisTokenPos = tokenPos_.at(index);
        if (thisTokenPos != std::string::npos)
        {
            return thisTokenPos;
        }
        size_t previousTokenPos = getTokenPos( index-1 );
        thisTokenPos = str_.find( D, previousTokenPos );
        assert( thisTokenPos != std::string::npos );
        ++thisTokenPos;
        return thisTokenPos;
    }

    void nonLazySplit()
    {
        boost::split(tokens_, str_, boost::is_any_of( std::string(1,D) ));
    }

    /**
     * \brief Internal storage
     */
    std::string str_;
    /**
     * \brief Lazy calculation of size()
     */
    bool sizeAlreadyEvaluated_;
    /**
     * 0 => empty str_
     * 1 => string has not yet been split
     * >1 => amount of chunks
     */
    size_t size_;
    /**
     * \brief Vector of positions for the beginning of each chunk
     */
    std::vector<size_t> tokenPos_;
    /**
     * \brief Cached sub-strings
     */
    std::vector<std::string> tokens_;
    //    std::vector< boost::iterator_range< std::string::const_iterator > > tokens_;
};


/**
 * \brief Fast, delimiter-separated-values manager.
 *
 * \tparam Head  the type of header
 * \tparam Line  the type of element
 * \tparam T     the type of each token
 * \tparam D     the separator
 */
template<typename Head, typename Line, typename T, char D>
class Dsv
{
    typedef T      (Line::*CallBack)(const char* name)const;
    typedef size_t (Head::*GetIndex)(const char* name)const;
    typedef void   (*Throw)(Line* vcf, const std::string& field, const std::string& line);

public:
    Dsv()
        : header_(0), line_(0), getIndex_(0), missing_(0), missingSize_(0)
    {
        ;
    }

    typedef T type;
    char delimiter() const
    {
        return D;
    }

    void setHeader(const Head* head)
    {
        header_   = &(*head);
    }
    void setLine(Line* line)
    {
        line_        = &(*line);
        missing_     = line_->MISSING;
        missingSize_ = (missing_) ? (sizeof(missing_) / sizeof(missing_[0])) : 0;
    }
    void setIndex(const GetIndex& getIndex)
    {
        getIndex_ = getIndex;
    }

    /** \brief Copy the string until the delimiter or until the end
     ** If found, the delimiter is discarded. the 'begin' iterator is moved
     ** either to the first position after the delimiter or to the 'end'.
     ** Parse the string directly, without using istreams for efficiency.
     ** \return pointer to the sub-string
     **/
    char* read(std::string::iterator& begin, const std::string::iterator& end,
               char delim = '\0') const
    {
        assert( header_ && line_ );
        if (! delim)
        {
            delim = D;
        }
        const std::string::iterator origin = begin;
        begin = std::find(begin, end, delim);
        *begin = '\0';
        if ((origin == begin) || isMissing(origin,begin) )
        {
            if (end != begin) ++begin;
            return const_cast<char*>(missing_);
        }
        else
        {
            if (end != begin) ++begin;
            return &(*origin);
        }
    }

    bool readList( std::string::iterator& begin,
                   const std::string::iterator& end,
                   std::vector<T>& tokenList,
                   CallBack callBack,
                   char delimiter1,
                   Throw dsvException,
                   bool fixSize,
                   size_t count = 0,
                   char delimiter2 = '\0')
    {
        assert( header_ && line_ );
        tokenList.clear();
        if (end == begin)
        {
            return false;
        }
        const std::string::iterator origin = begin;
        begin = std::find(begin, end, D);
        *begin = '\0';
        if ((origin == begin) || isMissing(origin,begin) )
        {
            // missing value
        }
        else
        {
            if (fixSize)
            {
                assert( NULL != getIndex_ );  // TODO: proper exception
                assert( count );              // TODO: proper exception
                // vector is initialized with whatever results from type-casting 0
                tokenList.resize( count, static_cast<T>(0) );
            }
            else
            {
                if ( count )  tokenList.reserve( count );
            }

            std::string::iterator stringBegin = origin;
            while (begin != stringBegin)
            {
                const std::string::iterator stringEnd = std::find(stringBegin, begin, delimiter1);
                *stringEnd = '\0';
                std::string::iterator equalSign = stringEnd;
                if (delimiter2)
                {
                    // when there is an '=' sign, the elements in vector tokenList need to be 'const char *'
                    assert( typeid(const char*) == typeid(T) );
                    // TODO: check that the ivalue is missing iif the field is expected to be a flag
                    equalSign = std::find(stringBegin, stringEnd, delimiter2);
                    *equalSign = '\0';
                }
                try
                {
                    T value = (line_->*callBack)(&(*stringBegin));
                    if (fixSize)
                    {
                        const size_t idx  = (header_->*getIndex_)(&(*stringBegin));
                        //const size_t idx = boost::bind(getIndex_, header_, _1)(&(*stringBegin));
                        tokenList.at(idx) = (stringEnd == equalSign) ? value : ((T) &(*(equalSign + 1)));
                    }
                    else
                    {
                        tokenList.push_back( value );
                    }
                }
                catch (illumina::common::OutOfBoundsException& ex)
                {
                    // if there isn't any '=', then we are overwriting 0 with 0 at the end
                    (*dsvException)( line_,
                                     std::string(stringBegin,equalSign),
                                     std::string(origin,end) );
                }
                catch (std::out_of_range& ex)
                {
                    throw;
                }
                stringBegin = (begin != stringEnd ? (stringEnd + 1) : begin);
            }
        }
        if (end == begin)
        {
            return false;
        }
        ++begin;
        return true;
    }

private:
    bool isMissing(const std::string::iterator& O, const std::string::iterator& B) const
    {
        // when missing_ is only 1 char we do it more efficiently
        return (1 == missingSize_) ? ((O + 1 == B) && *missing_ == *O)
               : !strcmp(&(*missing_), &(*O));
    }

    const Head* header_;
    Line* line_;
    GetIndex getIndex_;
//    char delimiter_;
    const char* missing_;
    size_t missingSize_;
};

/*
 * The key will be used as a persistent storage for id as a 0-terminated c-string
 *
 * It does create a string, so be careful with multi-threading
 */
class FastString
{
public:
    FastString() : str_(""), key_(0)   {}
    FastString(const FastString& rhs)
    {
        this->assign( rhs.string() );
    }
    FastString(const std::string& rhs)
    {
        this->assign( rhs );
    }

    const std::string& string() const
    {
        return str_;
    }
    std::string& string()
    {
        return str_;
    }
    const char* c_str() const
    {
        return key_.empty() ? NULL : &key_.front();
    }
    char* c_str()
    {
        return key_.empty() ? NULL : &key_.front();
    }

    bool empty() const
    {
        return key_.empty();
    }

    // not necessary to assert here, as assign() already keeps both in sync
    // size_t size() const                { size_t ret = key_.size() - 1; assert( str_.length() == ret ); return ret; }
    size_t size() const
    {
        return key_.size() - 1;
    }
    FastString& operator=(const FastString& rhs)
    {
        this->assign(rhs.string());
        return *this;
    }
    FastString& operator=(const std::string& rhs)
    {
        this->assign(rhs);
        return *this;
    }

private:
    // TODO check that the id is not empty
    void assign(const std::string& id)
    {
        str_ = id;
        key_.resize( id.size() + 1, '\0' );
        std::copy( id.begin(), id.end(), key_.begin() );
    }

    std::string str_;
    /// storage for the string as a "char *"
    std::vector<char> key_;
};

/*****************************************************************************/

