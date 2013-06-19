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

#include "blt_util/log.hh"
#include "common/Exceptions.hh"
#include "manta/SVLocusSet.hh"

#include "boost/foreach.hpp"

#include <algorithm>
#include <iostream>
#include <sstream>



void
SVLocusSet::
merge(SVLocus& inputLocus)
{

#ifdef DEBUG_SVL
    checkState();
    log_os << "SVLocusSet::merge inputLocus: " << inputLocus;
#endif

    const unsigned startLocusIndex(_loci.size());
    unsigned locusIndex(startLocusIndex);


    // test each node for intersection and insert/join to existing node:
    BOOST_FOREACH(SVLocusNode* inputNodePtr, inputLocus)
    {
        ins_type intersect;
        getNodeIntersect(inputNodePtr,intersect);

#ifdef DEBUG_SVL
        log_os << "SVLocusSet::merge inputNode: " << *inputNodePtr;
        log_os << "insersect_size: " << intersect.size() << "\n";
        BOOST_FOREACH(const ins_type::value_type& val, intersect)
        {
            log_os << "i-index: " << val.second << " " << *val.first << "\n";
        }
#endif

        if (intersect.empty())
        {
            // if no nodes intersect, then insert inputNode into a new locus:
            if (locusIndex>=_loci.size())
            {
                _loci.resize(locusIndex+1);
            }

            insertLocusNode(locusIndex,inputLocus,inputNodePtr);
        }
        else
        {
            // otherwise, merge inputNode into an existing locus

            {
                //
                // first step is to assign all intersect clusters to the lowest index number
                //

                // get lowest index number:
                BOOST_FOREACH(const ins_type::value_type& val, intersect)
                {
                    if (val.second < locusIndex)
                    {
                        locusIndex = val.second;
                    }
                }

                // assign all intersect clusters to the lowest index number
                combineLoci(startLocusIndex,locusIndex);
                BOOST_FOREACH(const ins_type::value_type& val, intersect)
                {
                    combineLoci(val.second,locusIndex);
                }

#ifdef DEBUG_SVL
                log_os << "Reassigned all intersecting nodes to index: " << locusIndex << " sli:" << startLocusIndex << "\n";
                checkState();
#endif

            }

            {
                //
                // second step is to add this inputNode into the graph
                //
                insertLocusNode(locusIndex,inputLocus,inputNodePtr);
#ifdef DEBUG_SVL
                log_os << "InsertedNewNode  in locus " << locusIndex << "\n";
                checkState();
#endif
            }

            {
                //
                // third step is to merge this inputNode with each intersecting inputNode,
                // and eliminate the intersecting node:
                //
                BOOST_FOREACH(const ins_type::value_type& val, intersect)
                {
#ifdef DEBUG_SVL
                    log_os << "MergeAndRemove: " << *val.first << "\n";
#endif
                    inputNodePtr->mergeNode(*val.first);
                    removeNode(val.first);
#ifdef DEBUG_SVL
                    log_os << "Finished: " << *val.first << "\n";
                    checkState();
#endif
                }

            }
        }
    }

#ifdef DEBUG_SVL
    checkState();
#endif
}



void
SVLocusSet::
getNodeIntersect(
    SVLocusNode* inputNodePtr,
    ins_type& intersect)
{
    typedef ins_type::iterator in_iter;

    intersect.clear();

#ifdef DEBUG_SVL
    log_os << "SVLocusSet::getNodeIntersect inputNode: " << *inputNodePtr;
    checkState();
#endif

    // get all existing nodes which intersect with this one:
    in_iter it(_inodes.lower_bound(inputNodePtr));

    // first look forward and extend to find all nodes which this inputNode intersects:
    for (in_iter it_fwd(it); it_fwd !=_inodes.end(); ++it_fwd)
    {
        SVLocusNode* nodePtr(it_fwd->first);
#ifdef DEBUG_SVL
        log_os << "FWD test: " << *nodePtr;
#endif
        if (! inputNodePtr->interval.isIntersect(nodePtr->interval)) break;
        intersect[nodePtr] = it_fwd->second;
    }

    // now find all intersecting nodes in reverse direction:
    for (in_iter it_rev(it); it_rev !=_inodes.begin(); )
    {
        --it_rev;
        SVLocusNode* nodePtr(it_rev->first);
#ifdef DEBUG_SVL
        log_os << "REV test: " << *nodePtr;
#endif
        if (! inputNodePtr->interval.isIntersect(nodePtr->interval)) break;
        intersect[nodePtr] = it_rev->second;
    }
}



void
SVLocusSet::
combineLoci(
    const unsigned fromIndex,
    const unsigned toIndex)
{
    assert(toIndex<_loci.size());

    if (fromIndex == toIndex) return;
    if (fromIndex>=_loci.size()) return;

    SVLocus& fromLocus(_loci[fromIndex]);
    if (fromLocus.empty()) return;

    SVLocus& toLocus(_loci[toIndex]);
    BOOST_FOREACH(SVLocusNode* fromNodePtr, fromLocus)
    {
        toLocus.copyNode(fromLocus,fromNodePtr);

        // update indexing structure:
        assert(_inodes[fromNodePtr] == fromIndex);
        _inodes[fromNodePtr] = toIndex;
    }
    fromLocus.clear();
}



void
SVLocusSet::
write(std::ostream& os) const
{
    os << "LOCUSSET_START\n";
    BOOST_FOREACH(const SVLocus& locus, _loci)
    {
        os << locus;
    }
    os << "LOCUSSET_END\n";
}



void
SVLocusSet::
checkState() const
{
    using namespace illumina::common;

    ins_type repeatCheck;

    unsigned locusIndex(0);
    unsigned nodeSize(0);
    BOOST_FOREACH(const SVLocus& locus, _loci)
    {
        nodeSize += locus.size();

        BOOST_FOREACH(SVLocusNode* nodePtr, locus)
        {
            assert(repeatCheck.count(nodePtr) == 0);

            ins_type::value_type val(std::make_pair(nodePtr,locusIndex));
            repeatCheck.insert(val);

            ins_type::const_iterator iter(_inodes.find(nodePtr));
            assert(iter != _inodes.end());
            assert(iter->second == locusIndex);
        }
        locusIndex++;
    }

    if (nodeSize != _inodes.size())
    {
        using namespace illumina::common;
        std::ostringstream oss;
        oss << "ERROR: SVLocusSet conflicting internal node counts. nodeSize: " << nodeSize << " inodeSize: " << _inodes.size() << "n";
        BOOST_THROW_EXCEPTION(PreConditionException(oss.str()));
    }
}

