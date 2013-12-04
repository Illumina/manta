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
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#pragma once

#include "blt_util/bam_record.hh"
#include "svgraph/GenomeInterval.hh"


/// This enables specification of different methods which
/// must traverse a range of reads in a bam file. By abstracting
/// multiple methods to this interface, we can accomplish multiple
/// tasks over a single pass of the BAM records while maintaining
/// isolation of methods
///
struct BamRegionProcessor
{
    virtual
    ~BamRegionProcessor() {}

    /// provide the index of the next bam file, must be called before switching files/samples
    ///
    /// for each bam index, return the requested interval for this operation,
    /// operations with closely related intervals will be compbined
    /// and the union of intervals will be processed
    virtual
    const GenomeInterval&
    nextBamIndex(
        const unsigned bamIndex) = 0;

    /// provide the next bam record
    virtual
    void
    processRecord(
        const bam_record& bamRead) = 0;
};
