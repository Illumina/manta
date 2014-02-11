// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2014 Illumina, Inc.
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

#include "common/Variant.hh"

#include "manta/SVCandidate.hh"
#include "svgraph/EdgeInfo.hh"
#include "svgraph/SVLocusSet.hh"

/*****************************************************************************/

class ObsRecord
{
public:
    ObsRecord();
    ObsRecord(const SVObservation& svObs);
    const SVObservation& getObs() const;

    enum DiscardReason { KEPT };
    bool discard(DiscardReason discardReason);
    bool discarded(DiscardReason& discardReason) const;

    bool operator<(const ObsRecord& obsRecordB) const;
    bool operator==(const ObsRecord& obsRecordB) const;

private:
    SVObservation _svObs;
    DiscardReason _discardReason;

    friend std::ostream& operator<<(std::ostream& ostrm,
                                    const ObsRecord& record);
};

std::ostream& operator<<(std::ostream& ostrm, const ObsRecord& record);

/*****************************************************************************/

typedef std::vector<ObsRecord> ObsRecordVec;
typedef ObsRecordVec::iterator ObsRecordVecIter;

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
    enum Outcome { UNKNOWN, WRITTEN, SPAWNED, IMPRECISE_NON_SPANNING,
                   LOW_SPANNING_COUNT, LOW_SOMATIC_SCORE
                 };

    SVLog(unsigned int index);
    unsigned int ind() const;
    void reportOutcome(const Outcome outcome);
    Outcome outcome() const;

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
    TruthTracker(
        const std::string& vcfFilePath,
        const bam_header_info& bamHeaderInfo);

    TruthTracker(
        const std::string& vcfFilePath,
        const SVLocusSet& cset);

private:
    bool loadTruth(const std::string& vcfFilePath,
                   const bam_header_info& bamHeaderInfo);

#ifdef EASY_ITER_OVER_NODE_EDGES
    bool evalLocusSet(const SVLocusSet& svLocusSet);
#endif

    void reportNumCands(unsigned int numCands);

public:
    // For ESL
    void addObservation(const SVObservation& svObs);

    // For GSC
    bool addEdge(const EdgeInfo& edge);
    bool discardEdge(const EdgeInfo& edge,
                     EdgeInfoRecord::DiscardReason reason);

    void reportNumCands(
        const unsigned numCands,
        const EdgeInfo& edge)
    {
        if (0 == numCands)
        {
            discardEdge(edge, EdgeInfoRecord::NO_DERIVED_SVS);
        }
        else
        {
            reportNumCands(numCands);
        }
    }

    void addCandSV();
    void reportNumAssembled(unsigned int numAssembled);
    void addAssembledSV();
    void reportOutcome(const SVLog::Outcome outcome);

    bool dumpObs();
    bool dumpMatchedEdgeInfo();
    bool dumpStats();

    void dumpAll()
    {
        dumpObs();
        dumpMatchedEdgeInfo();
        dumpStats();
    }

private:
    SVLocusSet _emptyLocusSet;
    const SVLocusSet& _cset;
    bool _hasTruth;
    VariantVec _trueVariantVec;

    typedef std::map<std::string, int32_t> ChromNameTidMap;
    ChromNameTidMap _chromNameTidMap;

    // SVObservation tracking

    typedef std::map<EdgeLog::VariantKey, ObsRecordVec> TruthObsVecMap;
    typedef TruthObsVecMap::iterator TruthObsVecMapIter;
    typedef TruthObsVecMap::value_type TruthObsVecMapEle;
    TruthObsVecMap _truthObsVecMap;

    unsigned long _numObs;


    // Edge tracking

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
