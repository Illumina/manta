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

///
/// \author Refactored by Richard Shaw from PUMA by Mauricio Varea
///

/*****************************************************************************/

#pragma once

#include "blt_util/thirdparty_push.h"

#include "boost/serialization/singleton.hpp"

#include "blt_util/thirdparty_pop.h"

#include "ChromosomeMetadata.hh"

/*****************************************************************************/

class ContigList : public boost::serialization::singleton<ContigList>
{
public:
    size_t getIndex(const char* name);
    const ChromosomeMetadata& getContig(size_t i) const
    {
        return chrList_.at(i);
    }

private:
    std::vector<ChromosomeMetadata> chrList_;
};

/*****************************************************************************/
