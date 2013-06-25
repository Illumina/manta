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

/// \author Chris Saunders
///

#pragma once

#include "blt_util/bam_util.hh"

#include "boost/serialization/string.hpp"
#include "boost/serialization/vector.hpp"

#include <iosfwd>
#include <string>
#include <vector>



/// minimal c++ bam header info
///
/// this class replicates the minimum information
/// from the bam header required to parse regions
/// (ie. chr1:20-30). It is friendlier to mem management
/// and serialization than using the samtools struct
///
struct bam_header_info
{

    bam_header_info()
    {}

    bam_header_info(const bam_header_t& header);

    bool
    operator==(const bam_header_info& rhs) const
    {
        const unsigned data_size(chrom_data.size());
        if (chrom_data.size() != rhs.chrom_data.size()) return false;
        for (unsigned i(0); i<data_size; ++i)
        {
            if (chrom_data[i] == rhs.chrom_data[i]) continue;
            return false;
        }
        return true;
    }

    template<class Archive>
    void serialize(Archive& ar, const unsigned /* version */)
    {
        ar& chrom_data;
    }

    struct chrom_info
    {
        chrom_info(
            const char* init_label = NULL,
            const unsigned init_length = 0) :
            label((NULL==init_label) ? "" : init_label ),
            length(init_length)
        {}

        bool
        operator==(const chrom_info& rhs) const
        {
            return ((label == rhs.label) && (length == rhs.length));
        }

        template<class Archive>
        void serialize(Archive& ar, const unsigned /* version */)
        {
            ar& label& length;
        }

        std::string label;
        unsigned length;
    };
    std::vector<chrom_info> chrom_data;
};


std::ostream&
operator<<(std::ostream& os, const bam_header_info& bhi);
