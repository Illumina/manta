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

#include <boost/format.hpp>
#include <boost/throw_exception.hpp>
#include <boost/assign/list_of.hpp>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include "common/Exceptions.hh"
#include "String.hh"

#include "VcfFile.hh"

/*****************************************************************************/

VcfFile::VcfFile(const std::string& pathStr,
                 const std::map<std::string, int32_t>& chromNameTidMap)
    : myVcfTypeMap(),
      myChromNameTidMap(chromNameTidMap),
      myVcfHeaderLoadedFlag(false),
      mySvTypeFieldId(-1), mySvLenFieldId(-1), myMateIdFieldId(-1)
{
    myVcfTypeMap.clear(); //cppcheck distraction
    myVcfTypeMap = boost::assign::map_list_of<std::string, Variant::Type>("BND", Variant::INTERTRANSLOC)("INV", Variant::INVERSION)("INS", Variant::INSERTION)("DEL", Variant::DELETION)("DUP:TANDEM", Variant::TANDUP)("compound", Variant::COMPLEX);

    myStrm.open(pathStr.c_str());

    if (!myStrm)
    {
        BOOST_THROW_EXCEPTION(InvalidOptionException(
                                  std::string("Failed to open `")
                                  + pathStr + "'"));
    }

    myPathStr = pathStr;
}

/*****************************************************************************/
static
bool parseTranslocAltStr(const std::string& altStr,
                         std::string& chromNameB, unsigned long& posB,
                         bool& isFwdA, bool& isFwdB)
{
    // const char leftBrkt('[');
    // const char rightBrkt(']');
    // const char chromPosSepar(':');

    // Hard-code to keep cppcheck happy.
    // SplitString<leftBrkt> leftBrktSplit(altStr);
    // SplitString<rightBrkt> rightBrktSplit(altStr);
    SplitString<'['> leftBrktSplit(altStr);
    SplitString<']'> rightBrktSplit(altStr);

    const bool brktLnotR(leftBrktSplit.size() == 3);

    if (!(brktLnotR || (rightBrktSplit.size() == 3)))
    {
        return false;
    }

    const bool baseIsAtStart(brktLnotR
                             ? (!leftBrktSplit[0].empty())
                             : (!rightBrktSplit[0].empty()));

    // Hard-code to keep cppcheck happy.
    // SplitString<chromPosSepar> posParts(brktLnotR
    //                                     ? leftBrktSplit[1]
    //                                     : rightBrktSplit[1]);
    SplitString<':'> posParts(brktLnotR
                              ? leftBrktSplit[1]
                              : rightBrktSplit[1]);

    if (posParts.size() != 2)
    {
        return false;
    }

    chromNameB = posParts[0];
    posB = boost::lexical_cast<unsigned long>(posParts[1]);

    isFwdA = baseIsAtStart;
    isFwdB = brktLnotR;

    if (brktLnotR && !baseIsAtStart)
    {
        isFwdA = true;
        isFwdB = false;
    }

    return true;
}

/*****************************************************************************/
static
void findFlankBounds(const unsigned long pos, const unsigned long flankSize,
                     pos_t& leftBound, pos_t& rightBound)
{
    leftBound = (pos > flankSize) ? (pos - flankSize) : 1;
    rightBound = pos + flankSize;
}

/*****************************************************************************/

bool VcfFile::getVariant(Variant& variant, bool& wasLast)
{
    if (!myVcfHeaderLoadedFlag)
    {
        myVcfHeaderLoadedFlag = loadHeader();
        mySvTypeFieldId = myVcfHeader.getInfoIndex("SVTYPE");
        mySvLenFieldId = myVcfHeader.getInfoIndex("SVLEN");
        myMateIdFieldId = myVcfHeader.getInfoIndex("MATEID");
    }

    VcfLine vcfLine;
    vcfLine.setFileId(0);
    vcfLine.setHeader(&myVcfHeader);
    vcfLine.init();

    bool novelVariant(false);
    Variant::Type variantType(Variant::UNKNOWN);

    while (!novelVariant)
    {
        myStrm >> vcfLine;

        if (myStrm.fail())
        {
            return false;
        }

        wasLast = myStrm.eof();

        const std::string typeStr(vcfLine.getInfo(mySvTypeFieldId));
        VcfTypeMapCIter vcfTypeCIter = myVcfTypeMap.find(typeStr);

        if (vcfTypeCIter != myVcfTypeMap.end())
        {
            variantType = vcfTypeCIter->second;
        }

        if (variantType == Variant::INTERTRANSLOC)
        {
            // Check if first line (proceed) or second line (discard).
            MateSeenMapIter mateSeenIter(myMateSeenMap.find(vcfLine.getId()));

            if (mateSeenIter == myMateSeenMap.end())
            {
                novelVariant = true;
            }
            else
            {
                mateSeenIter->second = true;
            }
        }
        else if (Variant::isStdType(variantType))
        {
            novelVariant = true;
            // i.e discard unknown, compound
        }
    }

    const std::string chromNameA(vcfLine.getChrom());
    const int32_t chromTidA(myChromNameTidMap.at(chromNameA));
    const unsigned long posA(vcfLine.getPos());

    int32_t chromTidB(chromTidA);
    unsigned long posB(0);
    bool isFwdA(true);
    bool isFwdB(true);

    if (Variant::isSingleChromType(variantType))
    {
        int svLen(boost::lexical_cast<int>(vcfLine.getInfo(mySvLenFieldId)));

        if (variantType == Variant::DELETION)
        {
            svLen = -svLen;
        }

        posB = posA + static_cast<unsigned long>(svLen);
    }
    else // Variant::INTERTRANSLOC
    {
        std::string chromNameB;

        if (!parseTranslocAltStr(vcfLine.getAlt(), chromNameB, posB,
                                 isFwdA, isFwdB))
        {
            BOOST_THROW_EXCEPTION(VcfException( (boost::format("Failed to parse transloc ALT `%s'") % vcfLine.getAlt()).str()));
        }

        chromTidB = myChromNameTidMap.at(chromNameB);

        // We have read the first line of a transloc & want to skip the other.
        const std::string mateIdStr(vcfLine.getInfo(myMateIdFieldId));
        myMateSeenMap[mateIdStr] = false;
    }

    const unsigned long flankSize(20);
    pos_t leftBoundA(0), rightBoundA(0), leftBoundB(0), rightBoundB(0);

    findFlankBounds(posA, flankSize, leftBoundA, rightBoundA);
    findFlankBounds(posB, flankSize, leftBoundB, rightBoundB);

    variant
        = Variant(variantType,
                  GenomeInterval(chromTidA, leftBoundA, rightBoundA), isFwdA,
                  GenomeInterval(chromTidB, leftBoundB, rightBoundB), isFwdB);

    // DEBUG
    // std::cerr << vcfLine << std::endl;
    // std::cerr << variant << std::endl;


    return true;
}

/*****************************************************************************/

bool VcfFile::getVariantVec(VariantVec& variantVec)
{
    variantVec.clear();
    Variant variant;
    bool wasLast(false);

    while (getVariant(variant, wasLast))
    {
        variantVec.push_back(variant);

        if (wasLast)
        {
            break;
        }
    }

    // Check that all transloc line mates found
    foreach (MateSeenMap::value_type mateSeenMapEle, myMateSeenMap)
    {
        if (!mateSeenMapEle.second)
        {
            std::cerr << "Failed to find mate with ID " << mateSeenMapEle.first
                      << std::endl;
        }
    }

    return true;
}

/*****************************************************************************/

bool VcfFile::loadHeader()
{
    myStrm >> myVcfHeader;

    // FIXME : error checking

    return true;
}

/*****************************************************************************/
