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
    EdgeInfoRecord();
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

    friend std::ostream& operator<<(std::ostream& ostrm,
                                    const EdgeInfoRecord& record);
};

std::ostream& operator<<(std::ostream& ostrm, const EdgeInfoRecord& record);

/*****************************************************************************/

typedef std::vector<EdgeInfoRecord> EdgeInfoRecordVec;
typedef EdgeInfoRecordVec::iterator EdgeInfoRecordVecIter;

/*****************************************************************************/

class SVLog
{
public:
    enum Outcome { UNKNOWN, WRITTEN, SPAWNED, IMPRECISE_SELF_EDGE,
                   LOW_PAIR_COUNT_SELF_EDGE, LOW_SOMATIC_SCORE };

    SVLog(unsigned int index);
    unsigned int ind() const;
    void reportOutcome(const Outcome outcome);
    const Outcome outcome() const;

protected:
    unsigned int _ind;
    Outcome _outcome;
};

std::ostream& operator<<(std::ostream& ostrm, const SVLog::Outcome outcome);

/*****************************************************************************/

class AssembledSVLog : public SVLog
{
public:
    AssembledSVLog(unsigned int ind);

    friend std::ostream& operator<<(std::ostream& ostrm,
                                    const AssembledSVLog& log);
};

std::ostream& operator<<(std::ostream& ostrm, const AssembledSVLog& log);

/*****************************************************************************/

class CandSVLog : public SVLog
{
public:
    CandSVLog(unsigned int ind);
    void reportNumAssembled(unsigned int numAssembled); // just a cross-check
    void addAssembledSV();
    void reportOutcome(const Outcome outcome);
    unsigned int numAssembledSVs() const;
    unsigned int numWrittenAssembledSVs() const;

private:
    unsigned int _expectedNumAssembled;
    unsigned int _numWrittenAssembled;

    typedef std::vector<AssembledSVLog> AssembledSVLogVec;
    AssembledSVLogVec _assembledSVLogVec;

    friend std::ostream& operator<<(std::ostream& ostrm, const CandSVLog& log);
};

std::ostream& operator<<(std::ostream& ostrm, const CandSVLog& log);

/*****************************************************************************/

class EdgeLog : public EdgeInfoRecord, public SVLog
{
public:
    EdgeLog(unsigned int ind, const EdgeInfo& edgeInfo);

    typedef unsigned int VariantKey;
    typedef std::vector<VariantKey> VariantKeyVec;
    typedef VariantKeyVec::iterator VariantKeyVecIter;
    void addVariantKey(VariantKey variantKey);
    const VariantKeyVec& variantKeyVec() const;

    void reportNumCands(unsigned int numCands); // just a cross-check
    void addCandSV();
    void reportNumAssembled(unsigned int numAssembled); // just a cross-check
    void addAssembledSV();
    void reportOutcome(const Outcome outcome);
    unsigned int numCandSVs() const;
    unsigned int numWrittenCandSVs() const;

private:
    VariantKeyVec _variantKeyVec;

    unsigned int _expectedNumCands;

    // Num candidate SVs that either were written 
    // or spawned written assembled SVs
    unsigned int _numWrittenCands;

    typedef std::vector<CandSVLog> CandSVLogVec;
    CandSVLogVec _candSVLogVec;

    friend std::ostream& operator<<(std::ostream& ostrm, const EdgeLog& log);
};

std::ostream& operator<<(std::ostream& ostrm, const EdgeLog& log);

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

    void reportNumCands(unsigned int numCands);
    void addCandSV();
    void reportNumAssembled(unsigned int numAssembled);
    void addAssembledSV();
    void reportOutcome(const SVLog::Outcome outcome);

    bool dumpMatchedEdgeInfo();
    bool dumpStats();

private:
    bool _hasTruth;
    std::string _vcfFilePath;
    VariantVec _trueVariantVec;

    typedef std::map<std::string, int32_t> ChromNameTidMap;
    ChromNameTidMap _chromNameTidMap;

    typedef std::map<EdgeLog::VariantKey, EdgeInfoRecordVec> TruthEdgeVecMap;
    typedef TruthEdgeVecMap::iterator TruthEdgeVecMapIter;
    typedef TruthEdgeVecMap::value_type TruthEdgeVecMapEle;
    TruthEdgeVecMap _truthEdgeVecMap;

    // Really want just EdgeInfo but that does not sort.
    typedef std::map<EdgeInfoRecord, EdgeLog> EdgeLogMap;
    typedef EdgeLogMap::iterator EdgeLogMapIter;
    typedef EdgeLogMap::value_type EdgeLogMapEle;
    typedef std::pair<EdgeLogMapIter, bool> EdgeLogMapIterBoolPr;
    EdgeLogMap _edgeLogMap;
    EdgeLogMapIter _lastAddedEdgeLogMapEleIter;

    // Vector of _edgeLogMap ele iterators in order added
    // -> can dump in edge index order without needing a sort.
    typedef std::vector<EdgeLogMapIter> EdgeLogMapIterVec;
    EdgeLogMapIterVec _edgeLogMapIndexVec;

    unsigned int _numEdgesSeen;
};

/*****************************************************************************/
