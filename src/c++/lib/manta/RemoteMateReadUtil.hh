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

#pragma once

#include "htsapi/bam_record.hh"

#include <string>


bool
isMateInsertionEvidenceCandidate(
    const bam_record& bamRead,
    const unsigned minMapq);


bool
isMateInsertionEvidenceCandidate2(
    const bam_record& bamRead,
    const bool isSearchForLeftOpen,
    const bool isSearchForRightOpen);


/// information recorded for reads where we need to grab the mate from a remote locus
///
/// typically these are chimeras with a MAPQ0 mate used to assemble a large insertion
///
struct RemoteReadInfo
{
    RemoteReadInfo(
        const bam_record& bamRead)
        : qname(bamRead.qname()),
          readNo(bamRead.read_no()==1 ? 2 : 1),
          tid(bamRead.mate_target_id()),
          pos(bamRead.mate_pos() - 1),
          localPos(bamRead.pos() - 1),
          readSize(bamRead.read_size()),
          isLocalFwd(bamRead.is_fwd_strand()),
          isFound(false),
          isUsed(false)
    {}

    bool
    operator<(
        const RemoteReadInfo& rhs) const
    {
        if (tid < rhs.tid) return true;
        if (tid == rhs.tid)
        {
            return (pos < rhs.pos);
        }
        return false;
    }

    std::string qname;
    int readNo; // this is read number of the target
    int tid;
    int pos;
    int localPos;
    int readSize;
    bool isLocalFwd;
    bool isFound;
    bool isUsed;
};
