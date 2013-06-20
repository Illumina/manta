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

#include "blt_util/log.hh"
#include "common/Exceptions.hh"
#include "manta/SVLocusSet.hh"

#include "boost/archive/binary_iarchive.hpp"
#include "boost/archive/binary_oarchive.hpp"
#include "boost/foreach.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>


unsigned
SVLocusSet::
getStartLocusIndex() const
{
    if(_emptyLoci.empty()) return _loci.size();
    return *_emptyLoci.begin();
}



void
SVLocusSet::
merge(SVLocus& inputLocus)
{

#ifdef DEBUG_SVL
    checkState();
    log_os << "SVLocusSet::merge inputLocus: " << inputLocus;
#endif

    const unsigned startLocusIndex(getStartLocusIndex());
    unsigned locusIndex(startLocusIndex);
    if(locusIndex<_loci.size())
    {
        assert(_loci[locusIndex].empty());
    }
    bool isFirst(true);

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
            moveIntersectToLowIndex(intersect,startLocusIndex,locusIndex);

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
                    mergeNodePtr(val.first,inputNodePtr);
                    removeNode(val.first);
#ifdef DEBUG_SVL
                    log_os << "Finished: " << *val.first << "\n";
                    checkState();
#endif
                }
            }
        }

        if(isFirst && (locusIndex == startLocusIndex))
        {
            _emptyLoci.erase(locusIndex);
            isFirst=false;
        }
    }

#ifdef DEBUG_SVL
    checkState();
#endif
}



void
SVLocusSet::
merge(SVLocusSet& inputSet)
{
    BOOST_FOREACH(SVLocus& locus, inputSet._loci)
    {
        merge(locus);
    }
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
moveIntersectToLowIndex(
        ins_type& intersect,
        const unsigned startLocusIndex,
        unsigned& locusIndex)
{
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
    _emptyLoci.insert(fromIndex);
}



void
SVLocusSet::
mergeNodePtr(SVLocusNode* fromPtr,
             SVLocusNode* toPtr)
{
    ins_type::iterator iter(_inodes.find(toPtr));
    assert(iter != _inodes.end());
    ins_type::value_type val(*iter);
    _inodes.erase(iter);
    toPtr->mergeNode(*fromPtr);
    _inodes.insert(val);
}




void
SVLocusSet::
dump(std::ostream& os) const
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
save(const char* filename) const
{
    using namespace boost::archive;

    assert(NULL != filename);
    std::ofstream ofs(filename, std::ios::binary);
    binary_oarchive oa(ofs);

    BOOST_FOREACH(const SVLocus& locus, _loci)
    {
        if(locus.empty()) continue;
        oa << locus;
    }
}



void
SVLocusSet::
load(const char* filename)
{
    using namespace boost::archive;

    clear();

    assert(NULL != filename);
    std::ifstream ifs(filename, std::ios::binary);
    binary_iarchive ia(ifs);

    SVLocus locus;
    while(ifs.peek() != EOF)
    {
        locus.clear();
        ia >> locus;
        if(locus.empty()) continue;
        _loci.push_back(locus);
    }

    reconstructIndex();
}



void
SVLocusSet::
reconstructIndex()
{
    _inodes.clear();
    _emptyLoci.clear();

    unsigned index(0);
    BOOST_FOREACH(SVLocus& locus, _loci)
    {
        BOOST_FOREACH(SVLocusNode* nodePtr, locus)
        {
            _inodes.insert(std::make_pair(nodePtr,index));
        }
        if(locus.empty()) _emptyLoci.insert(index);
        index++;
    }
}



void
SVLocusSet::
dumpIndex(std::ostream& os) const
{
    os << "SVLocusSet Index START\n";
    BOOST_FOREACH(const ins_type::value_type& in, _inodes)
    {
        os << "Index: " << in.second << " nodeptr: " << in.first << "\n";
    }
    os << "SVLocusSet Index END\n";
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
        locus.checkState();

        nodeSize += locus.size();

        BOOST_FOREACH(SVLocusNode * const nodePtr, locus)
        {
            {
                ins_type::const_iterator iter(repeatCheck.find(nodePtr));
                if(iter != repeatCheck.end())
                {
                    std::ostringstream oss;
                    oss << "ERROR: repeated node in SVLocusSet\n"
                        << "\tPrevObs index: " << iter->second << " node: " << *(iter->first)
                        << "\tThisObs index: " << locusIndex << " node: " << *nodePtr;
                    BOOST_THROW_EXCEPTION(PreConditionException(oss.str()));
                }
            }

            ins_type::value_type val(std::make_pair(nodePtr,locusIndex));
            repeatCheck.insert(val);

            ins_type::const_iterator citer(_inodes.find(nodePtr));
            if(citer == _inodes.end())
            {
                std::ostringstream oss;
                oss << "ERROR: locus node is missing from node index\n"
                    << "\tNode index: " << locusIndex << " node: " << *nodePtr;
                BOOST_THROW_EXCEPTION(PreConditionException(oss.str()));
            }
            if(citer->second != locusIndex)
            {
                std::ostringstream oss;
                oss << "ERROR: locus node is mis-assigned has conflicting index number in node index\n"
                    << "\tinode index_value: " << citer->second << "\n"
                    << "\tNode index: " << locusIndex << " node: " << *nodePtr;
                BOOST_THROW_EXCEPTION(PreConditionException(oss.str()));
            }
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

