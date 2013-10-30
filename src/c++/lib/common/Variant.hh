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

///
/// \author Richard Shaw
///

/*****************************************************************************/

#pragma once

#include <string>
#include <map>
#include <iostream>
#include <vector>

#include "svgraph/GenomeInterval.hh"

/*****************************************************************************/

class Variant
{
public:
    enum Type { UNKNOWN, INTERTRANSLOC, INVERSION, INSERTION, DELETION,
                TANDUP, COMPLEX
              };

    Variant();
    Variant(const Variant::Type& typeVal,
            const GenomeInterval& brkptAObj, bool isFwdAFlag,
            const GenomeInterval& brkptBObj, bool isFwdBFlag);

    Type type() const;
    const GenomeInterval& brkptA() const;
    bool isFwdA() const;
    const GenomeInterval& brkptB() const;
    bool isFwdB() const;

    static bool isStdType(const Type variantType);
    static bool isSingleChromType(const Type variantType);

private:
    Type myType;
    GenomeInterval myBrkptA;
    bool myIsFwdAFlag;
    GenomeInterval myBrkptB;
    bool myIsFwdBFlag;
};

/*****************************************************************************/

std::ostream& operator<<(std::ostream& ostrm, Variant::Type variantType);

std::ostream& operator<<(std::ostream& ostrm, Variant variant);

/*****************************************************************************/

typedef std::vector<Variant> VariantVec;
typedef VariantVec::const_iterator VariantVecCIter;

/*****************************************************************************/
