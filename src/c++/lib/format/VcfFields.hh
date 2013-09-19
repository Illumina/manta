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

/**
 ** The classes in this file support VCF v4.0, as defined in the
 ** 1000 Genomes project:
 **
 ** http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40
 **/

///
/// \author Refactored by Richard Shaw from PUMA by Mauricio Varea
///

/*****************************************************************************/

#pragma once

/*****************************************************************************/

class VcfFields
{
public:
    static const char * const FIXED[];
    static const char * const GENOTYPE[];
};

/*****************************************************************************/
