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

/*****************************************************************************/

#include "VcfFields.hh"

/*****************************************************************************/

const char* const VcfFields::FIXED[]
    = {"CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"};
const char* const VcfFields::GENOTYPE[] = {"FORMAT","SAMPLE"};

/*****************************************************************************/

unsigned int VcfFields::numFixedFields()
{
    return (sizeof(VcfFields::FIXED) / sizeof(VcfFields::FIXED[0]));
}

/*****************************************************************************/


