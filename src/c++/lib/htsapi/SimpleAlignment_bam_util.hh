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
///
///

#pragma once

#include "blt_util/SimpleAlignment.hh"
#include "htsapi/bam_record.hh"


SimpleAlignment
getAlignment(const bam_record& bamRead);
