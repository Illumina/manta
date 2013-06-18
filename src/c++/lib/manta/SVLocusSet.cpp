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
    unsigned locusIndex(_loci.size());

    // test each node for intersection and insert/join to existing node:
    const unsigned n_nodes(locus.graph.size());
    for(unsigned i(0);i<n_nodes;++i)
    {
        SVLocusNode* nodePtr(locus.graph[i].get());

        typedef ins_type::const_iterator in_citer;

        std::vector<inode> intersect;

        // get all existing nodes which intersect with this one:
        {
            inode search;
            search.nodePtr=nodePtr;
            in_citer it(std::lower_bound(_inodes.begin(), _inodes.end(), search, inodeCompare()));

            // first look forward and extend to find all nodes which this node intersects:
            for(in_citer it_fwd(it); it_fwd !=_inodes.end(); ++it_fwd)
            {
                inode in(*it_fwd);
                if(! nodePtr->interval.range.is_range_intersect(in.nodePtr->interval.range)) break;
                intersect.push_back(in);
            }

            // now find all intersecting nodes in reverse direction:
            for(in_citer it_rev(it); it_rev !=_inodes.begin(); )
            {
                --it_rev;
                inode in(*it_rev);
                if(! nodePtr->interval.range.is_range_intersect(in.nodePtr->interval.range)) break;
                intersect.push_back(in);
            }
        }

        if(intersect.empty())
        {
            // if no nodes intersect, then insert node into locus set:
            if(locusIndex>=_loci.size())
            {
                _loci.resize(locusIndex+1);
            }

            _loci[locusIndex].graph.push_back(locus.graph[i]);

            inode in;
            in.index=locusIndex;
            in.nodePtr=nodePtr;
            _inodes.insert(in);
        }
        else
        {
            // otherwise, merge node to the first intersecting node, then merge remaining nodes with each other
        }
    }
}



void
SVLocusSet::
write(std::ostream& os) const
{
    os << "FOO";
}
