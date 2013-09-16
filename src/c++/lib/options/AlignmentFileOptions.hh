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


#pragma once

#include <string>
#include <vector>


struct AlignmentFileOptions
{
    std::vector<std::string> alignmentFilename;
    std::vector<bool> isAlignmentTumor; ///< indicates which positions in the alignmnetFilename correspond to tumor
};
