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

///
/// \author Chris Saunders
///

#include "WriteSVToVcf.hh"

#include <iostream>



void
writeSVVcfHeaderPrefix(
    const char* referenceFilename,
    const char* version,
    std::ostream& os)
{
    os << "##fileformat=VCFv4.1\n";
    os << "##source=manta " << version << "\n";
    os << "##reference=file://" << referenceFilename << "\n";
    os << "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n";
    os << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n";
    os << "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">\n";
    os << "##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakend\">\n";
    os << "##INFO=<ID=BND_PAIR_SUPPORT,Number=1,Type=Integer,Description=\"Confidently mapped reads supporting this variant at this breakend (mapping may not be confident at remote breakend)\">\n";
    os << "##INFO=<ID=PAIR_SUPPORT,Number=1,Type=Integer,Description=\"Read pairs supporting this variant where both reads are confidently mapped\">\n";

    os << "##ALT=<ID=BND,Description=\"Translocation Breakend\">\n";
}



void
writeSVVcfHeaderSuffix(std::ostream& os)
{
    os << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
}


