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

/**
 ** \brief Checks variants being processed against a truth set.
 **
 ** Checks variants being processed against a truth set.
 **
 ** \author Richard Shaw
 **/

#pragma once

#include <string>
#include <vector>
#include <map>

#include "svgraph/EdgeInfo.hh"
#include "svgraph/SVLocusSet.hh"

#include "format/Variant.hh"

/*****************************************************************************/

class EdgeInfoRecord
{
public:
    EdgeInfoRecord(const EdgeInfo& edgeInfo);
    const EdgeInfo& getEdgeInfo() const;

    enum DiscardReason { KEPT, NO_DERIVED_SVS };
    bool discard(DiscardReason discardReason);
    bool discarded(DiscardReason& discardReason) const;

    bool operator<(const EdgeInfoRecord& edgeInfoRecordB) const;
    bool operator==(const EdgeInfoRecord& edgeInfoRecordB) const;

private:
    EdgeInfo _edgeInfo;
    DiscardReason _discardReason;
};

/*****************************************************************************/

typedef std::vector<EdgeInfoRecord> EdgeInfoRecordVec;
typedef EdgeInfoRecordVec::iterator EdgeInfoRecordVecIter;

/*****************************************************************************/

class TruthTracker
{
public:
    TruthTracker();
    bool loadTruth(const std::string& vcfFilePath,
                   const bam_header_info& bamHeaderInfo);
    bool evalLocusSet(const SVLocusSet& svLocusSet);
    bool addEdge(const EdgeInfo& edge, const SVLocusSet& cset);
    bool discardEdge(const EdgeInfo& edge,
                     EdgeInfoRecord::DiscardReason reason);

    bool dumpStats();

private:
    bool _hasTruth;
    std::string _vcfFilePath;
    VariantVec _trueVariantVec;

    typedef std::map<std::string, int32_t> ChromNameTidMap;
    ChromNameTidMap _chromNameTidMap;

    typedef unsigned int VariantKey;
    typedef std::vector<VariantKey> VariantKeyVec;
    typedef VariantKeyVec::iterator VariantKeyVecIter;

    typedef std::map<VariantKey, EdgeInfoRecordVec> TruthEdgeVecMap;
    typedef TruthEdgeVecMap::iterator TruthEdgeVecMapIter;
    typedef TruthEdgeVecMap::value_type TruthEdgeVecMapEle;
    TruthEdgeVecMap _truthEdgeVecMap;

    // Really want just EdgeInfo but that does not sort.
    typedef std::map<EdgeInfoRecord, VariantKeyVec> EdgeTruthVecMap;
    typedef EdgeTruthVecMap::iterator EdgeTruthVecMapIter;
    typedef EdgeTruthVecMap::value_type EdgeTruthVecMapEle;
    EdgeTruthVecMap _edgeTruthVecMap;

    unsigned int _numEdgesSeen;
};

/*****************************************************************************/
