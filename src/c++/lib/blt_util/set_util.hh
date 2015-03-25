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
/// \author Chris Saunders
///

#pragma once

#include <set>


/// implements B -= A
template <typename T>
void
inplaceSetSubtract(
    const std::set<T>& A,
    std::set<T>& B)
{
    typedef typename std::set<T>::const_iterator scit;
    scit ait(A.begin()), ait_end(A.end());
    scit bit(B.begin()), bit_end(B.end());
    while ( (bit != bit_end) && (ait != ait_end) )
    {
        if (*ait < *bit)
        {
            ++ait;
        }
        else
        {
            if (*ait == *bit)
            {
                const scit blast(bit);
                ++bit;
                B.erase(blast);
            }
            else
            {
                ++bit;
            }
        }
    }
}
