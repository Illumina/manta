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

///
/// \author Chris Saunders
///

#include <iterator>


namespace ALIGNPATH
{


/// convert the input path use seq match '=' and mismatch 'X'
///
template <typename symIter>
void
apath_add_seqmatch(
    symIter queryBegin,
    symIter queryEnd,
    symIter refBegin,
    symIter refEnd,
    path_t& apath)
{
    path_t apath2;

    const unsigned as(apath.size());
    for (unsigned segmentIndex(0); segmentIndex<as; ++segmentIndex)
    {
        const path_segment& ps(apath[segmentIndex]);
        if (is_segment_align_match(ps.type))
        {
            for (unsigned segmentPos(0); segmentPos<ps.length; ++segmentPos)
            {
                assert(queryBegin != queryEnd);
                assert(refBegin != refEnd);

                const bool isSeqMatch((*queryBegin) == (*refBegin));
                apath_append(apath2, ( isSeqMatch ? SEQ_MATCH : SEQ_MISMATCH));

                ++queryBegin;
                ++refBegin;
            }
        }
        else
        {
            apath2.push_back(ps);

            if (is_segment_type_read_length(ps.type)) std::advance(queryBegin,ps.length);
            if (is_segment_type_ref_length(ps.type)) std::advance(refBegin,ps.length);
        }
    }

    apath = apath2;
}

}
