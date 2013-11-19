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

#include "Variant.hh"

/*****************************************************************************/

std::ostream& operator<<(std::ostream& ostrm, Variant::Type variantType)
{
    switch (variantType)
    {
    case Variant::UNKNOWN:
        ostrm << "Unknown";
        break;
    case Variant::INTERTRANSLOC:
        ostrm << "translocation";
        break;
    case Variant::INVERSION:
        ostrm << "inversion";
        break;
    case Variant::INSERTION:
        ostrm << "insertion";
        break;
    case Variant::DELETION:
        ostrm << "deletion";
        break;
    case Variant::TANDUP:
        ostrm << "duplication";
        break;
    case Variant::COMPLEX:
        ostrm << "compound";
        break;
    }

    return ostrm;
}

/*****************************************************************************/

Variant::Variant()
    : myType(Variant::UNKNOWN), myIsFwdAFlag(true), myIsFwdBFlag(true)

{
    ;
}

/*****************************************************************************/

Variant::Variant(const std::string& idStr, const Variant::Type& typeVal,
                 const GenomeInterval& brkptAObj, bool isFwdAFlag,
                 const GenomeInterval& brkptBObj, bool isFwdBFlag)
    : myId(idStr), myType(typeVal),
      myBrkptA(brkptAObj), myIsFwdAFlag(isFwdAFlag),
      myBrkptB(brkptBObj), myIsFwdBFlag(isFwdBFlag)
{
    ;
}

/*****************************************************************************/

const std::string& Variant::id() const
{
    return myId;
}

/*****************************************************************************/

Variant::Type Variant::type() const
{
    return myType;
}

/*****************************************************************************/

const GenomeInterval& Variant::brkptA() const
{
    return myBrkptA;
}

/*****************************************************************************/

bool Variant::isFwdA() const
{
    return myIsFwdAFlag;
}

/*****************************************************************************/

const GenomeInterval& Variant::brkptB() const
{
    return myBrkptB;
}

/*****************************************************************************/

bool Variant::isFwdB() const
{
    return myIsFwdBFlag;
}

/*****************************************************************************/

bool Variant::isStdType(const Type variantType)
{
    return (isSingleChromType(variantType)
            || (variantType == INTERTRANSLOC));
}

/*****************************************************************************/

bool Variant::isSingleChromType(const Type variantType)
{
    switch (variantType)
    {
    case INVERSION:
    case INSERTION:
    case DELETION:
    case TANDUP:
        return true;

    default:
        break;
    }

    return false;
}

/*****************************************************************************/

std::ostream& operator<<(std::ostream& ostrm, Variant variant)
{
    ostrm << variant.id() << ' ' << variant.type()
          << ' ' << variant.brkptA() << ' ' << (variant.isFwdA() ? '+' : '-')
          << ' ' << variant.brkptB() << ' ' << (variant.isFwdB() ? '+' : '-');

    return ostrm;
}

/*****************************************************************************/
