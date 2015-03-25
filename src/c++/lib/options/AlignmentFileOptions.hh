// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2015 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//


#pragma once

#include <string>
#include <vector>


/// bam input file object shared by all programs which require these as input
struct AlignmentFileOptions
{
    std::vector<std::string> alignmentFilename;
    std::vector<bool> isAlignmentTumor; ///< indicates which positions in the alignmnetFilename correspond to tumor
};
