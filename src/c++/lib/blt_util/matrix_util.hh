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
// <https://github.com/downloads/sequencing/licenses/>.
//

// $Id: matrix_util.hh 743 2007-08-14 15:47:12Z ctsa $

/// \file

template <typename FloatType>
void
matrix_state_reduction(FloatType* reduced_mat,
                       const unsigned reduced_size,
                       const FloatType* full_mat,
                       const unsigned full_size,
                       const unsigned* reduction_map,
                       const bool is_average,
                       const bool is_exclude_diag)
{

    array_zero(reduced_mat,reduced_size*reduced_size);

    unsigned* reduced_add_count(0);
    if (is_average)
    {
        reduced_add_count = new unsigned[reduced_size*reduced_size];
        array_zero(reduced_add_count,reduced_size*reduced_size);
    }

    for (unsigned i(0); i<full_size; ++i)
    {
        const unsigned ri = reduction_map[i];
        for (unsigned j(0); j<full_size; ++j)
        {
            const unsigned rj = reduction_map[j];

            if (is_exclude_diag && i==j) continue;

            reduced_mat[rj+ri*reduced_size] += full_mat[j+i*full_size];
            if (is_average) reduced_add_count[rj+ri*reduced_size]++;
        }
    }

    if (is_average)
    {
        for (unsigned ri(0); ri<reduced_size; ++ri)
        {
            for (unsigned rj(0); rj<reduced_size; ++rj)
            {
                const unsigned val = reduced_add_count[rj+ri*reduced_size];
                if (val) reduced_mat[rj+ri*reduced_size] /= static_cast<FloatType>(val);
            }
        }
        delete [] reduced_add_count;
    }
}
