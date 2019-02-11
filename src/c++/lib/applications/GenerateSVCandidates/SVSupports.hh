//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2019 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

/// \file
/// \author Xiaoyu Chen
///

#pragma once

#include <set>
#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include <memory>

#include "manta/SVCandidateSetData.hh"
#include "htsapi/bam_streamer.hh"
#include "htsapi/bam_dumper.hh"
#include "htsapi/bam_record_util.hh"


typedef std::shared_ptr<bam_record> bam_record_ptr;
typedef std::shared_ptr<bam_streamer> bam_streamer_ptr;
typedef std::shared_ptr<bam_dumper> bam_dumper_ptr;


/// Records a single read that support one or more SVs for evidence-BAM output
struct SupportRead
{
    typedef std::map<std::string, std::set<std::string>> SV_supportType_t;

    void addNewSV(
        const std::string& svId, const std::string& supportType)
    {
        if (SVs.find(svId) == SVs.end())
        {
            std::set<std::string> supportTypeSet;
            SVs[svId] = supportTypeSet;
        }

        SVs[svId].insert(supportType);
    }

    bool
    operator<(
        const SupportRead& rhs) const
    {
        if (tid < rhs.tid) return true;
        if (tid == rhs.tid)
        {
            return (pos < rhs.pos);
        }
        return false;
    }



    int tid = -1;
    int pos = 0;
    SV_supportType_t SVs;
};

std::ostream&
operator<<( std::ostream& os, const SupportRead& suppRd);


/// Records a single fragment (read1 & read2)
/// that support one or more SVs for evidence-BAM output,
/// indicating the evidence type
struct SupportFragment
{

    void setReads(
        const bam_record& bamRead)
    {
        if (bamRead.is_first())
        {
            read1.tid = bamRead.target_id();
            read1.pos = bamRead.pos();
            read2.tid = bamRead.mate_target_id();
            read2.pos = bamRead.mate_pos();
        }
        else if (bamRead.is_second())
        {
            read1.tid = bamRead.mate_target_id();
            read1.pos = bamRead.mate_pos();
            read2.tid = bamRead.target_id();
            read2.pos = bamRead.pos();
        }
    }

    void addSpanningSupport(
        const std::string& svID)
    {
        read1.addNewSV(svID, "PR");
        read2.addNewSV(svID, "PR");
    }

    void addSplitSupport(
        const bool isRead1,
        const std::string& svID)
    {
        if (isRead1)
        {
            read1.addNewSV(svID, "SR");
            read2.addNewSV(svID, "SRM");
        }
        else
        {
            read2.addNewSV(svID, "SR");
            read1.addNewSV(svID, "SRM");
        }
    }


    SupportRead read1;
    SupportRead read2;
};

std::ostream&
operator<<( std::ostream& os, const SupportFragment& suppFrg);

typedef std::map<std::string, SupportFragment> support_fragments_t;


/// Records all supporting fragments
/// that support one or more SVs for evidence-BAM output
struct SupportFragments
{
    SupportFragment& getSupportFragment(
        const bam_record& bamRead)
    {
        const std::string qname(bamRead.qname());

        // create a new entry in the map
        if (supportFrags.find(qname) == supportFrags.end())
        {
            SupportFragment newFrag;
            newFrag.setReads(bamRead);
            supportFrags[qname] = newFrag;
        }

        return supportFrags[qname];
    }


    SupportFragment& getSupportFragment(
        const SVCandidateSetSequenceFragment& seqFrag)
    {
        // Tentatively add an assertion
        // \TODO add the logic if only supplementary (NOT primary) reads being set
        assert(seqFrag.read1.isSet() || seqFrag.read2.isSet());

        const SVCandidateSetRead& read(seqFrag.read1.isSet() ? seqFrag.read1 : seqFrag.read2);

        const std::string qname(read.bamrec.qname());

        // create a new entry in the map
        if (supportFrags.find(qname) == supportFrags.end())
        {
            SupportFragment newFrag;
            newFrag.setReads(read.bamrec);
            supportFrags[qname] = newFrag;
        }

        return supportFrags[qname];
    }

    support_fragments_t supportFrags;
};

std::ostream&
operator<<( std::ostream& os, const SupportFragments& suppFrgs);


struct SupportSamples
{
    SupportFragments& getSupportFragments(
        const unsigned index)
    {
        assert(index < supportSamples.size());
        return supportSamples[index];
    }

    std::vector<SupportFragments> supportSamples;
};

std::ostream&
operator<<( std::ostream& os, const SupportSamples& suppSmps);

void
processBamRecords(bam_streamer& origBamStream,
                  const GenomeInterval& interval,
                  const support_fragments_t& supportFrags,
                  bam_dumper& bamDumper);

void
writeSupportBam(const bam_streamer_ptr& origBamStream,
                const SupportFragments& svSupportFrags,
                const bam_dumper_ptr& supportBamDumper);

