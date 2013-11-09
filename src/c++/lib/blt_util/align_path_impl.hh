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
/// \author Chris Saunders
///

#include <iterator>

#include "blt_util/blt_exception.hh"


namespace ALIGNPATH
{


/// convert the input path use seq match '=' and mismatch 'X'
///
template <typename symIter>
void
apath_add_seqmatch(
    const symIter queryBegin,
    const symIter queryEnd,
    const symIter refBegin,
    const symIter refEnd,
    path_t& apath)
{
    path_t apath2;

    symIter queryIndex(queryBegin);
    symIter refIndex(refBegin);

    const unsigned as(apath.size());
    for (unsigned segmentIndex(0); segmentIndex<as; ++segmentIndex)
    {
        const path_segment& ps(apath[segmentIndex]);
        if (is_segment_align_match(ps.type))
        {
            for (unsigned segmentPos(0); segmentPos<ps.length; ++segmentPos)
            {
                if (queryIndex >= queryEnd)
                {
                    throw blt_exception("apath_add_seqmatch: past end of query\n");
                }

                if (refIndex >= refEnd)
                {
                    throw blt_exception("apath_add_seqmatch: past end of reference\n");
                }

                bool isSeqMatch((*queryIndex) == (*refIndex));
                // N always counts as match
                if (*queryIndex == 'N') isSeqMatch = true;
                apath_append(apath2, ( isSeqMatch ? SEQ_MATCH : SEQ_MISMATCH));

                ++queryIndex;
                ++refIndex;
            }
        }
        else
        {
            apath2.push_back(ps);

            if (is_segment_type_read_length(ps.type)) std::advance(queryIndex,ps.length);
            if (is_segment_type_ref_length(ps.type)) std::advance(refIndex,ps.length);
        }
    }

    apath = apath2;
}

}
