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
/// \author Chris Saunders
///

#pragma once

#include "blt_util/blt_exception.hh"

#include <algorithm>
#include <sstream>
#include <vector>


/// two predefined options for the ValClear type parameter to RangeMap:
template <typename T>
struct ZeroT
{
    void
    operator()(T& val) const
    {
        val = 0;
    }
};

template <typename T>
struct ClearT
{
    void
    operator()(T& val) const
    {
        val.clear();
    }
};


/// this object provides map-like storage for a set of positions which are
/// assumed to cluster in a small range
///
/// in practice this is very similar to an unbounded circular buffer, but heavily customized
/// for our specific application
///
template <typename KeyType, typename ValType, typename ValClear = ZeroT<ValType>>
struct RangeMap
{
    ///\TODO automate this w/ static assert/concepts:
    //Keytype must implement operator < +/-

    RangeMap() :
        _isEmpty(true),
        _minKeyIndex(0),
        _data(_minChunk),
        _occup(_minChunk,false)
    {}

    void
    clear()
    {
        _isEmpty=true;
        _minKeyIndex=0;
        std::fill(_occup.begin(),_occup.end(),false);
    }

    bool
    empty() const
    {
        return _isEmpty;
    }

    bool
    isKeyPresent(
        const KeyType& k) const
    {
        return (! ((_isEmpty) || (k < _minKey) || (k > _maxKey) || (! _occup[getKeyIndex(k)])));
    }

    ValType&
    getRef(
        const KeyType& k)
    {
        if (_isEmpty)
        {
            _minKey=k;
        }
        else if (k < _minKey)
        {
            expand(_maxKey-k);
            const unsigned dataSize(_data.size());
            _minKeyIndex = ((_minKeyIndex + dataSize)-(_minKey-k) ) % dataSize;
            _minKey = k;
        }

        if ((_isEmpty) || (k > _maxKey))
        {
            expand(k-_minKey);
            _maxKey = k;
            _isEmpty=false;
        }


        const unsigned kindex(getKeyIndex(k));
        if (! _occup[kindex])
        {
            _clearFunc(_data[kindex]);
            _occup[kindex]=true;
        }

        return _data[kindex];
    }

    const ValType&
    getConstRef(
        const KeyType& k) const
    {
        enforceKeyPresent(k);
        return _data[getKeyIndex(k)];
    }

    const ValType&
    getConstRefDefault(
        const KeyType& k,
        const ValType& defaultVal) const
    {
        if (! isKeyPresent(k)) return defaultVal;
        return _data[getKeyIndex(k)];
    }

    void
    erase(
        const KeyType& k)
    {
        enforceKeyPresent(k);
        const unsigned kindex(getKeyIndex(k));
        _occup[kindex]=false;

        if (k != _minKey) return;

        // we have to shift minKey up to the next valid value:
        const unsigned keySize(_maxKey-_minKey);
        const unsigned dataSize(_data.size());
        for (unsigned offset(1); offset<=keySize; ++offset)
        {
            if (! _occup[(_minKeyIndex+offset) % dataSize]) continue;
            _minKeyIndex=(_minKeyIndex+offset) % dataSize;
            _minKey += offset;
            return;
        }

        _isEmpty=true;
    }

    /// erase all contents with keys sorting less than or equal to k
    void
    eraseTo(
        const KeyType& k)
    {
        // special cases:
        if (_isEmpty) return;
        if (_minKey > k) return;
        if (_maxKey <= k)
        {
            clear();
            return;
        }

        // clear _occup values up to k:
        const unsigned kindex(getKeyIndex(k));
        if (_minKeyIndex<=kindex)
        {
            std::fill(_occup.begin()+_minKeyIndex,_occup.begin()+kindex+1,false);
        }
        else
        {
            std::fill(_occup.begin(),_occup.begin()+kindex+1,false);
            std::fill(_occup.begin()+_minKeyIndex,_occup.end(),false);
        }


        // we have to shift minKey up to the next valid value:
        const unsigned keySize(_maxKey-_minKey);
        const unsigned dataSize(_data.size());
        for (unsigned offset(1); offset<=keySize; ++offset)
        {
            if (! _occup[(_minKeyIndex+offset) % dataSize]) continue;
            _minKeyIndex=(_minKeyIndex+offset) % dataSize;
            _minKey += offset;
            return;
        }
        _isEmpty=true;
    }

#ifdef DEBUG_RMAP
    /// debug dumper:
    void
    dump(const char* msg, std::ostream& os) const
    {
        os << "rangeMap dump: " << msg << "\n"
           << "\tempty: " << _isEmpty << "\n"
           << "\tminKeyIndex: " << _minKeyIndex << "\n"
           << "\tminKey: " << _minKey << "\n"
           << "\tmaxKey: " << _maxKey << "\n"
           << "\tdatasize: " << _data.size() << "\n";
    }
#endif

private:
    /// assumes key has already been validated!
    unsigned
    getKeyIndex(
        const KeyType& k) const
    {
        return ((_minKeyIndex + (k-_minKey)) % _data.size());
    }

    /// rotate data so that minKeyIndex is 0
    void
    normRotate()
    {
        std::rotate(_data.begin(),_data.begin()+_minKeyIndex,_data.end());
        std::rotate(_occup.begin(),_occup.begin()+_minKeyIndex,_occup.end());
        _minKeyIndex=0;
    }

    // expand to larger of 2x current size or minSize+minChunk:
    void
    expand(
        const unsigned minSize)
    {
        if (minSize <= _data.size()) return;
        const unsigned newSize(std::max(static_cast<unsigned>(2*_data.size()),minSize+_minChunk));
        normRotate();
        _data.resize(newSize);
        _occup.resize(newSize,false);
    }

    void
    enforceKeyPresent(
        const KeyType& k) const
    {
        if (isKeyPresent(k)) return;
        std::ostringstream oss;
        oss << "Attempting to retrieve an invalid key '" << k << "'\n";
        throw blt_exception(oss.str().c_str());
    }

    static constexpr unsigned _minChunk = 1024;

    bool _isEmpty;
    unsigned _minKeyIndex;
    KeyType _minKey;
    KeyType _maxKey;
    std::vector<ValType> _data;
    std::vector<bool> _occup;
    ValClear _clearFunc;
};
