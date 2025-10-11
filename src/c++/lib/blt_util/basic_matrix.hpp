//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2025 Illumina, Inc.
//
// This program is licensed under the terms of the Polyform strict license
//
// ***As far as the law allows, the software comes as is, without
// any warranty or condition, and the licensor will not be liable
// to you for any damages arising out of these terms or the use
// or nature of the software, under any kind of legal claim.***
//
// You should have received a copy of the PolyForm Strict License 1.0.0
// along with this program.  If not, see <https://polyformproject.org/licenses/strict/1.0.0>.
//
//

#pragma once

#include <vector>

/// very simple matrix implementation, row major
template <typename T>
struct basic_matrix {
  typedef typename std::vector<T>         data_t;
  typedef typename data_t::iterator       iterator;
  typedef typename data_t::const_iterator const_iterator;

  basic_matrix(const unsigned rowCount = 0, const unsigned colCount = 0)
    : _colCount(colCount), _data(rowCount * colCount)
  {
  }

  void resize(const unsigned rowCount, const unsigned colCount)
  {
    _colCount = colCount;
    _data.resize(rowCount * colCount);
  }

  T& val(const unsigned row, const unsigned col) { return _data[(row * _colCount + col)]; }

  const T& val(const unsigned row, const unsigned col) const { return _data[(row * _colCount + col)]; }

  bool empty() { return _data.empty(); }

  size_t size() { return _data.size(); }

  iterator begin() { return _data.begin(); }

  const_iterator begin() const { return _data.begin(); }

  iterator end() { return _data.end(); }

  const_iterator end() const { return _data.end(); }

private:
  unsigned       _colCount;
  std::vector<T> _data;
};
