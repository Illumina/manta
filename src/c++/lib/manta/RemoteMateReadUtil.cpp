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

///
/// \author Ole Schulz-Trieglaff
/// \author Chris Saunders
///

#include "RemoteMateReadUtil.hh"



bool
isMateInsertionEvidenceCandidate(
    const bam_record& bamRead,
    const unsigned minMapq)
{
    if (! bamRead.is_paired()) return false;
    if (bamRead.is_unmapped() || bamRead.is_mate_unmapped()) return false;

    if (bamRead.map_qual() < minMapq) return false;

    if (bamRead.target_id() < 0) return false;
    if (bamRead.mate_target_id() < 0) return false;

    if (bamRead.target_id() != bamRead.mate_target_id()) return true;

    /// TODO: better candidate definition based on fragment size distro:
    static const int minSize(10000);
    return (std::abs(bamRead.pos()-bamRead.mate_pos()) >= minSize);
}



bool
isMateInsertionEvidenceCandidate2(
    const bam_record& bamRead,
    const bool isSearchForLeftOpen,
    const bool isSearchForRightOpen)
{
    if ((! isSearchForLeftOpen) && (! bamRead.is_fwd_strand())) return false;
    if ((! isSearchForRightOpen) && bamRead.is_fwd_strand()) return false;
    return true;
}
