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

#include "manta/SVLocusSet.hh"

#include <algorithm>
#include <iostream>



void
SVLocusSet::
merge(SVLocus& locus)
{
    // test each node for intersection and insert/join to existing node:
    const unsigned n_nodes(locus.graph.size());
    for(unsigned i(0);i<n_nodes;++i)
    {
        SVLocusNode* nodePtr(locus.graph[i].get());

        typedef in_type::const_iterator in_citer;
        in_citer it(std::lower_bound(_inodes.begin(), _inodes.end(), nodePtr, nodeCompare()));

        std::vector<SVLocusNode*> intersect;

        // first look forward and extend to find all nodes which this node intersects:
        for(in_citer it_fwd(it); it_fwd !=_inodes.end(); ++it_fwd)
        {
            SVLocusNode* inodePtr(*it_fwd);
            if(! nodePtr->interval.range.is_range_intersect(inodePtr->interval.range)) break;
            intersect.push_back(inodePtr);
        }

        // now find all intersecting nodes in reverse direction:
        for(in_citer it_rev(it); it_rev !=_inodes.begin(); )
        {
            --it_rev;
            SVLocusNode* inodePtr(*it_rev);
            if(! nodePtr->interval.range.is_range_intersect(inodePtr->interval.range)) break;
            intersect.push_back(inodePtr);
        }
    }
}



void
SVLocusSet::
write(std::ostream& os) const
{
    os << "FOO";
}
