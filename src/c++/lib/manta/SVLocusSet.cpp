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



std::ostream&
operator<<(std::ostream& os, const SVLocusSet::NodeAddressType& a)
{
    os << a.first << ":" << a.second;
    return os;
}



void
SVLocusSet::
merge(SVLocus& inputLocus)
{
    using namespace illumina::common;

#ifdef DEBUG_SVL
    checkState();
    log_os << "SVLocusSet::merge inputLocus: " << inputLocus;
#endif

    const LocusIndexType startLocusIndex(insertLocus(inputLocus));
    const SVLocus& startLocus(_loci[startLocusIndex]);
    LocusIndexType headLocusIndex(startLocusIndex);

    // test each node for intersection and insert/join to existing node:
    const unsigned nodeCount(startLocus.size());
    for(unsigned nodeIndex(0);nodeIndex<nodeCount;++nodeIndex)
    {
        LocusSetIndexerType intersect(*this);
        getNodeIntersect(startLocusIndex, nodeIndex, intersect);

#ifdef DEBUG_SVL
        log_os << "SVLocusSet::merge inputNode: " << startLocus.getNode(nodeIndex);
        log_os << "insersect_size: " << intersect.size() << "\n";
        BOOST_FOREACH(const LocusSetIndexerType::value_type& val, intersect)
        {
            log_os << "i-index: " << val << " node: " <<  getNode(val) << "\n";
        }
#endif

        if (headLocusIndex != startLocusIndex)
        {
            if(intersect.empty())
            {
                NodeAddressType val(std::make_pair(startLocusIndex,nodeIndex));
                std::ostringstream oss;
                oss << "ERROR: no intersecting nodes found during merge\n"
                    << "\tsearch node: " << val << " " << getNode(val)
                    << "\thli: " << headLocusIndex << "\n";
                BOOST_THROW_EXCEPTION(PreConditionException(oss.str()));
            }
            if(1==intersect.size()) continue;
        }
        else
        {
            if(intersect.empty()) continue;
        }

        // merge inputNode into an existing locus
        moveIntersectToLowIndex(intersect,startLocusIndex,headLocusIndex);

        intersect.clear();
        getNodeIntersect(startLocusIndex, nodeIndex, intersect);

        assert(! intersect.empty());

#ifdef DEBUG_SVL
        log_os << "intersect2_size: " << intersect.size() << "\n";
        BOOST_FOREACH(const LocusSetIndexerType::value_type& val, intersect)
        {
            log_os << "i2-index: " << val << " node: " <<  getNode(val) << "\n";
        }
#endif

        // merge overlapping nodes in order from highest nodeid to lowest, so that the
        // merge process does not invalid nodeids of higher value
        //
        // to do this, we first need to find a node corresponding the input node, and sort
        // the remaining nodes by nodeIndex:
        //
        NodeAddressType inputSuperAddy;
        std::vector<NodeAddressType> nodeIndices;
        {
            bool isInputSuperFound(false);
            const known_pos_range& inputRange(getLocus(startLocusIndex).getNode(nodeIndex).interval.range);

            BOOST_FOREACH(const LocusSetIndexerType::value_type& val, intersect)
            {
                assert(val.first==headLocusIndex);

                // one node must be a superset of the input node, find this and store separately:
                if((! isInputSuperFound) && getNode(val).interval.range.is_superset_of(inputRange))
                {
                    inputSuperAddy=val;
                    isInputSuperFound=true;
                    continue;
                }
                nodeIndices.push_back(val);
            }
            assert(isInputSuperFound);
            std::sort(nodeIndices.rbegin(),nodeIndices.rend());
        }

        // merge this inputNode with each intersecting inputNode,
        // and eliminate the intersecting node:
        //
        NodeAddressType mergeTargetAddy(inputSuperAddy);
        BOOST_FOREACH(NodeAddressType nodeAddy, nodeIndices)
        {
            if(nodeAddy<mergeTargetAddy) std::swap(nodeAddy,mergeTargetAddy);
#ifdef DEBUG_SVL
            log_os << "MergeAndRemove: " << nodeAddy << "\n";
#endif
            mergeNodePtr(nodeAddy,mergeTargetAddy);
            removeNode(nodeAddy);
#ifdef DEBUG_SVL
            log_os << "Finished: " << nodeAddy << "\n";
            checkState();
#endif
        }
    }

    if(startLocusIndex != headLocusIndex)
    {
#ifdef DEBUG_SVL
        log_os << "clearLocusIndex: " << startLocusIndex << "\n";
#endif

        clearLocus(startLocusIndex);
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
        const LocusIndexType locusIndex,
        const NodeIndexType nodeIndex,
        LocusSetIndexerType& intersect)
{
    typedef LocusSetIndexerType::iterator in_iter;

    intersect.clear();

#ifdef DEBUG_SVL
    log_os << "SVLocusSet::getNodeIntersect inputNode: " << locusIndex << ":" << nodeIndex << " " << getNode(std::make_pair(locusIndex,nodeIndex));
    checkState();
#endif

    // get all existing nodes which intersect with this one:
    const NodeAddressType inputAddy(std::make_pair(locusIndex,nodeIndex));
    in_iter it(_inodes.lower_bound(inputAddy));
    const SVLocusNode& inputNode(getNode(inputAddy));

    // first look forward and extend to find all nodes which this inputNode intersects:
    for (in_iter it_fwd(it); it_fwd !=_inodes.end(); ++it_fwd)
    {
        if(it_fwd->first == locusIndex) continue;
#ifdef DEBUG_SVL
        log_os << "FWD test: " << (*it_fwd) << " " << getNode(*it_fwd);
#endif
        if (! inputNode.interval.isIntersect(getNode(*it_fwd).interval)) break;
        intersect.insert(*it_fwd);
#ifdef DEBUG_SVL
        log_os << "FWD insert: " << (*it_fwd) << "\n";
#endif
    }

    // now find all intersecting nodes in reverse direction:
    for (in_iter it_rev(it); it_rev !=_inodes.begin(); )
    {
        --it_rev;
        if(it_rev->first == locusIndex) continue;
#ifdef DEBUG_SVL
        log_os << "REV test: " << (*it_rev) << " " << getNode(*it_rev);
#endif
        if (! inputNode.interval.isIntersect(getNode(*it_rev).interval)) break;
        intersect.insert(*it_rev);
#ifdef DEBUG_SVL
        log_os << "REV insert: " << (*it_rev) << "\n";
#endif
    }
}



void
SVLocusSet::
getRegionIntersect(
        const int32_t tid,
        const int32_t beginPos,
        const int32_t endPos,
        LocusSetIndexerType& intersect)
{
    const LocusIndexType startLocusIndex(insertLocus(SVLocus()));
    const NodeIndexType nodeIndex = getLocus(startLocusIndex).addNode(tid,beginPos,endPos);

    getNodeIntersect(startLocusIndex,nodeIndex,intersect);

    clearLocus(startLocusIndex);
}



void
SVLocusSet::
moveIntersectToLowIndex(
        LocusSetIndexerType& intersect,
        const LocusIndexType startLocusIndex,
        LocusIndexType& locusIndex)
{
    const unsigned startHeadLocusIndex(locusIndex);

    // assign all intersect clusters to the lowest index number
    const bool isClearSource(startLocusIndex!=startHeadLocusIndex);

    // get lowest index number that is not startLocusIndex:
    bool isFirst(true);
    BOOST_FOREACH(const LocusSetIndexerType::value_type& val, intersect)
    {
        if ((!isFirst) && (val.first >= locusIndex)) continue;
        locusIndex = val.first;
        isFirst=false;
    }

    combineLoci(startHeadLocusIndex,locusIndex,isClearSource);
    BOOST_FOREACH(const LocusSetIndexerType::value_type& val, intersect)
    {
        combineLoci(val.first,locusIndex);
    }

#ifdef DEBUG_SVL
    log_os << "Reassigned all intersecting nodes to index: " << locusIndex << " shli: " << startHeadLocusIndex << " sli:" << startLocusIndex << "\n";
    checkState();
#endif
}



void
SVLocusSet::
combineLoci(
    const LocusIndexType fromIndex,
    const LocusIndexType toIndex,
    const bool isClearSource)
{
    assert(toIndex<_loci.size());

#ifdef DEBUG_SVL
    log_os << "combineLoci: from: " << fromIndex << " toIndex: " << toIndex << " isClear:" << isClearSource << "\n";
#endif

    if (fromIndex == toIndex) return;
    if (fromIndex>=_loci.size()) return;

    SVLocus& fromLocus(_loci[fromIndex]);
    if (fromLocus.empty()) return;

    SVLocus& toLocus(_loci[toIndex]);
    toLocus.copyLocus(fromLocus);
    if(isClearSource) clearLocus(fromIndex);
}



LocusIndexType
SVLocusSet::
insertLocus(
    const SVLocus& inputLocus)
{
    LocusIndexType locusIndex(0);
    if(_emptyLoci.empty())
    {
        static const unsigned maxIndex(std::numeric_limits<LocusIndexType>::max());
        locusIndex=_loci.size();
        assert(locusIndex<maxIndex);
        _loci.resize(locusIndex+1);
    }
    else
    {
        locusIndex=(*_emptyLoci.begin());
        assert(_loci[locusIndex].empty());
        _emptyLoci.erase(locusIndex);
    }

    SVLocus& locus(_loci[locusIndex]);
    observe_notifier(locus);
    locus.updateIndex(locusIndex);
    locus.copyLocus(inputLocus);
    return locusIndex;
}



void
SVLocusSet::
mergeNodePtr(NodeAddressType fromPtr,
             NodeAddressType toPtr)
{
#ifdef DEBUG_SVL
    log_os << "MergeNode: from: " << fromPtr << " to: " << toPtr << " fromLocusSize: " << getLocus(fromPtr.first).size() << "\n";
#endif
    LocusSetIndexerType::iterator iter(_inodes.find(toPtr));
    assert(iter != _inodes.end());
    assert(fromPtr.first == toPtr.first);
    getLocus(fromPtr.first).mergeNode(fromPtr.second,toPtr.second);
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
dumpRegion(std::ostream& os,
        const int32_t tid,
        const int32_t beginPos,
        const int32_t endPos)
{
    LocusSetIndexerType intersect(*this);
    getRegionIntersect(tid,beginPos,endPos,intersect);
    BOOST_FOREACH(const LocusSetIndexerType::value_type& val, intersect)
    {
        os << "SVNode from LocusIndex: " << val.second << " "<< val.first;
    }
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
        const LocusIndexType locusIndex(size());
        _loci.push_back(locus);
        SVLocus& locusCopy(_loci.back());
        observe_notifier(locusCopy);
        locusCopy.updateIndex(locusIndex);
    }

    reconstructIndex();
    checkState();
}



void
SVLocusSet::
reconstructIndex()
{
    _inodes.clear();
    _emptyLoci.clear();

    LocusIndexType locusIndex(0);
    BOOST_FOREACH(SVLocus& locus, _loci)
    {
        const unsigned nodeCount(locus.size());
        for(NodeIndexType nodeIndex(0);nodeIndex<nodeCount;++nodeIndex)
        {
            _inodes.insert(std::make_pair(locusIndex,nodeIndex));
        }
        if(locus.empty()) _emptyLoci.insert(locusIndex);
        locusIndex++;
    }
}



void
SVLocusSet::
dumpIndex(std::ostream& os) const
{
    os << "SVLocusSet Index START\n";
    BOOST_FOREACH(const LocusSetIndexerType::value_type& in, _inodes)
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

    unsigned locusIndex(0);
    unsigned totalNodeCount(0);
    BOOST_FOREACH(const SVLocus& locus, _loci)
    {
        locus.checkState();

        const unsigned nodeCount(locus.size());
        totalNodeCount += nodeCount;

        for(NodeIndexType nodeIndex(0);nodeIndex<nodeCount;++nodeIndex)
        {
            LocusSetIndexerType::const_iterator citer(_inodes.find(std::make_pair(locusIndex,nodeIndex)));
            if(citer == _inodes.end())
            {
                std::ostringstream oss;
                oss << "ERROR: locus node is missing from node index\n"
                    << "\tNode index: " << locusIndex << " node: " << getNode(std::make_pair(locusIndex,nodeIndex));
                BOOST_THROW_EXCEPTION(PreConditionException(oss.str()));
            }
            if((citer->first != locusIndex) || (citer->second != nodeIndex))
            {
                std::ostringstream oss;
                oss << "ERROR: locus node is has conflicting index number in node index\n"
                    << "\tinode index_value: " << citer->first << ":" << citer->second << "\n"
                    << "\tNode index: " << locusIndex << ":" << locusIndex << " node: " << getNode(std::make_pair(locusIndex,nodeIndex));
                BOOST_THROW_EXCEPTION(PreConditionException(oss.str()));
            }
        }
        locusIndex++;
    }

    if (totalNodeCount != _inodes.size())
    {
        using namespace illumina::common;
        std::ostringstream oss;
        oss << "ERROR: SVLocusSet conflicting internal node counts. totalNodeCount: " << totalNodeCount << " inodeSize: " << _inodes.size() << "n";
        BOOST_THROW_EXCEPTION(PreConditionException(oss.str()));
    }
}

