// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
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

#include "blt_util/log.hh"
#include "htsapi/bam_dumper.hh"

#include <cstdlib>

#include <iostream>



bam_dumper::
bam_dumper(const char* filename,
           const bam_header_t* header)
{
    _bfp = samopen(filename, "wb", header);

    if (NULL == _bfp)
    {
        log_os << "ERROR: Failed to open output BAM file: " << filename << "\n";
        exit(EXIT_FAILURE);
    }
}
