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

/// \file
/// \author Chris Saunders
///

#pragma once

#include <set>

/// implements B -= A
template <typename T>
void inplaceSetSubtract(const std::set<T>& A, std::set<T>& B)
{
  auto ait(A.cbegin());
  auto ait_end(A.cend());
  auto bit(B.cbegin());
  auto bit_end(B.cend());
  while ((bit != bit_end) && (ait != ait_end)) {
    if (*ait < *bit) {
      ++ait;
    } else {
      if (*ait == *bit) {
        const auto blast(bit);
        ++bit;
        B.erase(blast);
      } else {
        ++bit;
      }
    }
  }
}
