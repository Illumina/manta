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

#pragma once

#include "manta/SVLocus.hh"

#include <iosfwd>
#include <set>
#include <vector>


// A set of non-overlapping SVLocus objects
//
struct SVLocusSet {

    /// merge new locus into the set:
    ///
    /// locus is destroyed in this process
    ///
    void
    merge(SVLocus& locus);

    void
    write(std::ostream& os) const;

private:
    // contains the full set of loci
    std::vector<SVLocus> _loci;

    struct nodeCompare {

        bool
        operator()(const SVLocusNode* a,
                   const SVLocusNode* b) const
        {
            return ((a->interval)<(b->interval));
        }
    };

    // provides an intersection search of non-overlapping nodes:
    typedef std::set<SVLocusNode*, nodeCompare> in_type;
    in_type _inodes;
};
