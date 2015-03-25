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

/// \file
///
/// an efficient (and slightly unsafe) class for basic tab-delimited files, etc...
///

/// \author Chris Saunders
///

#pragma once

#include <iosfwd>


struct istream_line_splitter
{

    istream_line_splitter(std::istream& is,
                          const unsigned line_buf_size=8*1024,
                          const char word_seperator='\t',
                          const unsigned max_word=0)
        : _is(is)
        , _line_no(0)
        , _n_word(0)
        , _buf_size(line_buf_size)
        , _sep(word_seperator)
        , _max_word(max_word)
        , _buf(new char[_buf_size])
    {

        if ((0==_max_word) || (MAX_WORD_COUNT < _max_word))
        {
            _max_word=MAX_WORD_COUNT;
        }
    }

    ~istream_line_splitter()
    {
        if (nullptr != _buf)
        {
            delete [] _buf;
            _buf=nullptr;
        }
    }

    unsigned
    n_word() const
    {
        return _n_word;
    }

    /// returns false for regular end of input:
    bool
    parse_line();

    // recreates the line before parsing
    void
    write_line(std::ostream& os) const;

    // debug output, which provides line number and other info before calling write_line
    void
    dump(std::ostream& os) const;


    enum { MAX_WORD_COUNT = 50 };
    char* word[MAX_WORD_COUNT];
private:

    void
    increase_buffer_size();

    std::istream& _is;
    unsigned _line_no;
    unsigned _n_word;
    unsigned _buf_size;
    char _sep;
    unsigned _max_word;
    char* _buf;
};



#if 0
{
    //usage example:
    istream_line_splitter dparse(data_is);

    while (dparse.parse_line())
    {
        static const unsigned col_count(46);
        if (dparse.n_word()!=col_count)
        {
            std::ostringstream oss;
            oss << "ERROR: unexpected number of columns in paired export line:\n\n";
            dparse.dump(oss);
            throw blt_exception(oss.str().c_str());
        }

        for (unsigned i(1); (i+1)<col_count; ++i)
        {
            dparse.word[i][strlen(dparse.word[i])] = sep;
        }
        const char* nocompress_segment(dparse.word[0]);
        const char* compress_segment(dparse.word[1]);

        /// ....etc
    }
}
#endif


