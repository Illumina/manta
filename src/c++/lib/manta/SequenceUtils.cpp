
#include "common/SequenceUtils.hh"

#include "boost/scoped_array.hpp"
#include "boost/lexical_cast.hpp"

#include <cstdlib>
#include <iostream>

const double invlog10(1./(std::log(10.)));

/// returns log(1+x), switches to special libc function when abs(x) is small
///
static
double
log1p_switch(const double x) {
    // better number??
    static const double smallx_thresh(0.01);
    if (std::abs(x)<smallx_thresh) {
        return log1p(x);
    } else {
        return std::log(1+x);
    }
}


SequenceUtils::
ReadScorer::
ReadScorer()
    : _qmin(phredScoreOffset) {

#ifdef DEBUG_SU
    std::cout << "Filling logpcorrectratio table" << std::endl;
#endif
    _logpcorrectratio[_qmin] = 0;
    for (int i(_qmin+1); i<MAX_Q; ++i) {
        const double eprob(SequenceUtils::convertPhredToProbError(i-phredScoreOffset));
#ifdef DEBUG_SU
        std::cout << "i=" << i << " " << log1p_switch(-eprob) << " " << std::log(eprob/3.) << std::endl;
#endif
        _logpcorrectratio[i] = log1p_switch(-eprob) - std::log(eprob/3.);
    }

#ifdef DEBUG_SU
    std::cout << "Readscorer dumping _logpcorrectratio : " << std::endl;
    for (int i(_qmin); i<MAX_Q; ++i) {
        std::cout << "_logpcorrectratio[" <<  i << "] = " << _logpcorrectratio[i] << std::endl;
    }
#endif
}


double
SequenceUtils::
getBreakPointHypothesisScore(const Alignment& alignment,
                             const int breakPointPos,
                             const bool enteringBreakPoint) {
    return ReadScorer::get().
           getBreakPointHypothesisScore(alignment, breakPointPos,
                                        enteringBreakPoint);
}

/*****************************************************************************/

int
SequenceUtils::
getMostProbableBreakPoint(const Alignment& alignment,
                          bool& poorAlignIsAtStart)
{
    return ReadScorer::get().
           getMostProbableBreakPoint(alignment, poorAlignIsAtStart);
}

/*****************************************************************************/

double
SequenceUtils::
ReadScorer::
getAlignmentScore(const char* qualityString,
                  const char* type) const {
    int posInRead = 0;
    int tmp = 0;
    double alignScore(0.);
    bool inGap(false);

    if (!isValidQualString(qualityString)) {
#ifdef DEBUG_SU
        std::cout << "!isValidQualString(qualityString)" << std::endl;
#endif
        return alignScore;
    }

#ifdef DEBUG_SU
    std::cout << "getAlignmentScore type=" << type << std::endl;
    std::cout << "getAlignmentScore qual=" << qualityString << std::endl;
#endif
    const unsigned nt(strlen(type));
    for (unsigned i=0; i<nt; ++i)
    {
#ifdef DEBUG_SU
        std::cout << "getAlignmentScore : i = " << i << " alignScore = " << alignScore << std::endl;
#endif
        // ASCII: 48 = 0, 57=9
        if ((int) type[i] <= 57 && (int) type[i] >= 48)
        {
#ifdef DEBUG_SU
            std::cout << "getAlignmentScore : i = " << i << " alignScore = " << alignScore << std::endl;
#endif
            tmp *= 10;
            tmp += (type[i] - '0');
        }
        else if (type[i] == 'A' || type[i] == 'C' || type[i] == 'G' || type[i] == 'T')
        {
            posInRead += tmp;
            tmp=0;
            // Bases in gap -> dels : ignoring for scoring purposes
            //                        & do not advance the read position.
            if (!inGap) {
                // mismatch
#ifdef DEBUG_SU
                std::cerr << "getAlignmentScore: " << posInRead << " " << _logpcorrectratio[static_cast<int>(qualityString[posInRead])]
                          << " " << qualityString[posInRead] << " " << static_cast<int>(qualityString[posInRead]) << std::endl;
#endif
                alignScore +=  _logpcorrectratio[static_cast<int>(qualityString[posInRead])];
                ++posInRead;
            }
        } else if (type[i] == 'N') {
            posInRead += tmp;
            tmp=0;
            // what should N in the reference be?? for now it counts as a match
            ++posInRead;
        } else if (type[i]=='^') {
            posInRead += tmp;
            tmp=0;
            inGap = true;
        } else if (type[i] == '$') {
            // If tmp is non-zero here, it indicates an insertion.
            // FIXME : probably should update score using open & extend penalties.
            // Currently just advancing the read position.
            posInRead += tmp;
            tmp=0;
            inGap = false;
        } else {
            std::cerr << "ERROR: Unexpected match descriptor: '" << type << "'\n";
            exit(EXIT_CODE_FAILURE);
        }
    }
    return alignScore;
}






double
SequenceUtils::
ReadScorer::
getSemiAlignedReadMetric(const Alignment& alignment) const {
    return ((getAlignmentScore(alignment.getQuality().c_str(),
                               alignment.getMatchDescriptor().c_str()))
            * invlog10);
}

double
SequenceUtils::
ReadScorer::
getSemiAlignedReadMetric(const BamTools::BamAlignment& alignment) const {
    std::string matchDescriptor;
    alignment.GetTag("XD", matchDescriptor);
    return ((getAlignmentScore(alignment.Qualities.c_str(),matchDescriptor.c_str()))*invlog10);
}


double
SequenceUtils::
getSemiAlignedReadMetric(const Alignment& alignment) {
    return ReadScorer::get().getSemiAlignedReadMetric(alignment);
}

double
SequenceUtils::
getSemiAlignedReadMetric(const BamTools::BamAlignment& alignment) {
    return ReadScorer::get().getSemiAlignedReadMetric(alignment);
}

/*****************************************************************************/

bool
SequenceUtils::
getReadPairInfo(const Alignment& alignment,
                const unsigned readLen,
                const unsigned partnerReadLen,
                bool& isChimeric,
                long& leftPos,
                unsigned int& insertSize,
                RelOrient& relOrient) {
#ifdef DEBUG_SU
    std::cerr << "# Getting Read Pair Info!" << std::endl;
#endif

    // const std::string& chrom(alignment.getChromosome());
    const Match& partnerMatch(alignment.getPartnerMatch());
    const std::string& partnerChrom(partnerMatch.getChromosome());

    const long pos(alignment.getPosition());
    long partnerPos(partnerMatch.getPosition());

    Match::Strand strand(alignment.getStrand());
    Match::Strand partnerStrand(partnerMatch.getStrand());

    const bool isRead1(alignment.getReadNumber() == 1);

    if (!(((strand == Match::Forward)
           || (strand == Match::Reverse))
          && ((partnerStrand == Match::Forward)
              || (partnerStrand == Match::Reverse)))) {
        // Interested only in pairs where both partners are aligned.
        return false;
    }
#ifdef DEBUG_SU
    std::cerr << "# Before chimeric check!" << std::endl;
#endif
    // Possible chimeric translocation.
    if (!partnerChrom.empty()) {
        // Partner Chrom is blank if the two chromosomes are the same.
        isChimeric = true;
        leftPos = pos;
        insertSize = 0;
#ifdef DEBUG_SU
        std::cerr << "# IS chimeric in get read pair info!" << std::endl;
#endif
        return true;
    }
#ifdef DEBUG_SU
    std::cerr << "# NOT chimeric in get read pair info!" << std::endl;
#endif
    isChimeric = false;

    // Same chromosome -> partnerPos is offset -> convert to absolute position.
    partnerPos += pos;

    if (partnerPos == pos) {
        // Not interested in dimers/palindromes but values as per SVfinder.pl.
        leftPos = pos;
        insertSize = (isRead1 ? readLen : partnerReadLen);
        return false;
    }

    const bool isLeftRead(pos < partnerPos);
    leftPos = isLeftRead ? pos : partnerPos;

    relOrient = RelOrient(isRead1 ? pos : partnerPos,
                          isRead1 ? strand : partnerStrand,
                          isRead1 ? partnerPos : pos,
                          isRead1 ? partnerStrand : strand);

    insertSize = (isLeftRead
                  ? (partnerPos - pos + partnerReadLen)
                  : (pos - partnerPos + readLen));

    return true;
}

/*****************************************************************************/

bool
SequenceUtils::
ReadScorer::
isValidQualString(const char* const qual) const {

    const char* q(qual);
    while (*q != '\0') {
        const int q33val(*q);
        if (q33val < _qmin or q33val >= MAX_Q) {
            std::cerr << "ERROR:: Invalid quality value: " << static_cast<int>(*q)-phredScoreOffset
                      << " in position: " << (q-qual+1) << " of quality string: " << qual << "\n";
            return false;
        }
        ++q;
    }
    return true;
}

double SequenceUtils::
ReadScorer::
getBreakPointHypothesisScore(const Alignment& alignment, const unsigned breakPointPos, const bool enteringBreakPoint) const {

    unsigned j=0;
    const char* qualityString = alignment.getQuality().c_str();
    const char* type = alignment.getMatchDescriptor().c_str();
    unsigned lengthOfAlignment = strlen (qualityString);
    double alignScore = 0.0;
#ifdef DEBUG_SU
    std::cerr <<  "In getBreakPointHypothesisScore: qualityString " << qualityString << "\t" << "type " <<  type << "\t" << "breakPointPos " << breakPointPos << "\n";
#endif
    if (breakPointPos==0) {
        if (enteringBreakPoint)
            alignScore = ReadScorer::getAlignmentScoreAfterBreakPoint(qualityString,type);
        else
            alignScore = ReadScorer::getAlignmentScore(qualityString,type);
#ifdef DEBUG_SU
        std::cerr << "alignScore = " << 	alignScore << std::endl;
#endif
        return alignScore;
    } else if (breakPointPos>=lengthOfAlignment) {
        if (enteringBreakPoint)
            alignScore = ReadScorer::getAlignmentScore(qualityString,type);
        else
            alignScore = ReadScorer::getAlignmentScoreAfterBreakPoint(qualityString,type);
#ifdef DEBUG_SU
        std::cerr << "alignScore = " << alignScore << std::endl;
#endif
        return alignScore;
    }

    //if all matches
    const unsigned nt(strlen(type));
    bool allMatches = true;
#ifdef DEBUG_SU
    std::cerr << "NT: " << nt << ", type: " << type << std::endl;
#endif
    for (unsigned i=0; i<nt; ++i)
    {
        if ((int) type[i] > 57 ||  (int) type[i] < 48) { //not between 0 and 9
            allMatches = false;
            break;
        }
    }
#ifdef DEBUG_SU
    std::cerr << " lengthOfAlignment = " << lengthOfAlignment << std::endl;
#endif
    const std::string qualityStringBeforeBreakPoint(qualityString,qualityString+breakPointPos);
    const std::string qualityStringAfterBreakPoint(qualityString+breakPointPos,qualityString+lengthOfAlignment);

    if (allMatches) {
        int numAfterBreakPoint = lengthOfAlignment-breakPointPos;
        char numBeforeBreakPointStr[5];
        char numAfterBreakPointStr[5];

        sprintf(numBeforeBreakPointStr, "%d", breakPointPos);
        sprintf(numAfterBreakPointStr, "%d", numAfterBreakPoint);

        if (enteringBreakPoint)	 {
            alignScore = ReadScorer::getAlignmentScoreAfterBreakPoint(qualityStringAfterBreakPoint.c_str(),numAfterBreakPointStr);
        } else {
            alignScore = ReadScorer::getAlignmentScore(qualityStringBeforeBreakPoint.c_str(),numBeforeBreakPointStr);
        }
#ifdef DEBUG_SU
        std::cerr << "Returning here: " << alignScore << std::endl;
#endif
        return (alignScore);
    } //all matches, like 100


    //figure out where to break the type string
    j=0;
    int posInRead = 0;
    //const unsigned nt(strlen(type));
    unsigned i=0, tmp=0;
    boost::scoped_array<char> typePtr(new char[2*(nt+1)]);
    char* typeBeforeBreakPoint(typePtr.get());
    char* typeAfterBreakPoint(typePtr.get()+(nt+1));
    bool breakPointFound = false;
    unsigned typeCutPosBefore, typeCutPosAfter;
    //char* typeComplete = (char*) malloc (lengthOfAlignment+1);

    for (i=0; i<nt; ++i) {
        if ((int) type[i] <= 57 && (int) type[i] >= 48)  {
            tmp *= 10;
            tmp += (type[i] - '0');
            //if (posInRead+tmp>=lengthOfAlignment) {break;}
        }  else {

            if (posInRead+tmp==breakPointPos) {
                unsigned typeCutPos=i-1;
                for (j=0; j<nt; j++) {
                    if (j<=typeCutPos) {
                        typeBeforeBreakPoint[j]=type[j];
                    } else {
                        typeAfterBreakPoint[j-typeCutPos-1]=type[j];
                    }
                }
                typeBeforeBreakPoint[typeCutPos+1]= '\0';
                typeAfterBreakPoint[j-typeCutPos-1]='\0';
                breakPointFound=true;
                break;
            } else if (posInRead+tmp>breakPointPos) {

                unsigned numBaseBeforeBreakPoint = breakPointPos-posInRead;
                unsigned numBaseAfterBreakPoint = posInRead+tmp-breakPointPos;

                unsigned startPosAfterBreakPoint = 0;
                unsigned digitsOfTmp = (int)(log(tmp)/log(10))+1;

                typeCutPosBefore=i-digitsOfTmp;

                //if (tmp>9) typeCutPosBefore=i-2; else typeCutPosBefore=i-1;
                typeCutPosAfter = i;

                if (numBaseAfterBreakPoint>9) {
                    for (j=0; j<digitsOfTmp; j++) {
                        typeAfterBreakPoint[j]=(int)(numBaseAfterBreakPoint/pow(10,(digitsOfTmp-j-1)))+48;
                        numBaseAfterBreakPoint -= boost::lexical_cast<unsigned>(floor(numBaseAfterBreakPoint/pow(10,(digitsOfTmp-j-1)))*pow(10,(digitsOfTmp-j-1)));
                    }
                    startPosAfterBreakPoint += digitsOfTmp;
                    //typeAfterBreakPoint[0]=(int)(numBaseAfterBreakPoint/10)+48;
                    //typeAfterBreakPoint[1]=(int)(numBaseAfterBreakPoint%10)+48;
                    //startPosAfterBreakPoint+=2;
                } else {
                    typeAfterBreakPoint[0]=(int)(numBaseAfterBreakPoint)+48;
                    startPosAfterBreakPoint++;
                }
                for (j=0; j<nt; j++) {
                    if (j<typeCutPosBefore) {
                        typeBeforeBreakPoint[j]=type[j];
                    } else if (j>=typeCutPosAfter) {
                        typeAfterBreakPoint[startPosAfterBreakPoint+j-typeCutPosAfter]=type[j];
                    }
                }
                typeAfterBreakPoint[startPosAfterBreakPoint+j-typeCutPosAfter]='\0';
                if (numBaseBeforeBreakPoint>9) {
                    for (j=0; j<digitsOfTmp; j++) {
                        typeBeforeBreakPoint[typeCutPosBefore+j] = (int)(numBaseBeforeBreakPoint/pow(10,digitsOfTmp-j-1))+48;
                        numBaseBeforeBreakPoint -= boost::lexical_cast<unsigned>(floor(numBaseBeforeBreakPoint/pow(10,digitsOfTmp-j-1))*pow(10,digitsOfTmp-j-1));
                    }
                    typeBeforeBreakPoint[typeCutPosBefore+digitsOfTmp]='\0';
                    //typeBeforeBreakPoint[typeCutPosBefore]=(int)(numBaseBeforeBreakPoint/10)+48;
                    //typeBeforeBreakPoint[typeCutPosBefore+1]=(int)(numBaseBeforeBreakPoint%10)+48;
                    //typeBeforeBreakPoint[typeCutPosBefore+2]='\0';
                } else {
                    typeBeforeBreakPoint[typeCutPosBefore]=(int)(numBaseBeforeBreakPoint)+48;
                    typeBeforeBreakPoint[typeCutPosBefore+1]='\0';
                }

                breakPointFound=true;
                break;
            }
            posInRead += tmp;
            posInRead++;
            tmp = 0;
#ifdef DEBUG_SU
            std::cerr << "posInRead= " << posInRead << std::endl;
#endif
        }
    }
#ifdef DEBUG_SU
    std::cerr << "BreakPointFound: " << breakPointFound << std::endl;
#endif

    if (!breakPointFound) {
        unsigned numBaseBeforeBreakPoint = breakPointPos-posInRead;
        unsigned numBaseAfterBreakPoint = posInRead+tmp-breakPointPos;
        unsigned typeCutPosBefore /*, typeCutPosAfter*/;
        unsigned startPosAfterBreakPoint = 0;

        unsigned digitsOfTmp = (int)(log(tmp)/log(10))+1;
        typeCutPosBefore=i-digitsOfTmp;
        //if (tmp>9) typeCutPosBefore=i-2; else typeCutPosBefore=i-1;
        //typeCutPosAfter = i;

        if (numBaseAfterBreakPoint>9) {
            for (j=0; j<digitsOfTmp; j++) {
                typeAfterBreakPoint[j]=(int)(numBaseAfterBreakPoint/pow(10,digitsOfTmp-j-1))+48;
                numBaseAfterBreakPoint -= boost::lexical_cast<unsigned>(floor(numBaseAfterBreakPoint/pow(10,digitsOfTmp-j-1))*pow(10,digitsOfTmp-j-1));
            }
            startPosAfterBreakPoint += digitsOfTmp;

            //typeAfterBreakPoint[0]=(int)(numBaseAfterBreakPoint/10) + 48;
            //typeAfterBreakPoint[1]=(int)(numBaseAfterBreakPoint%10) + 48;
            //startPosAfterBreakPoint+=2;
        } else {
            typeAfterBreakPoint[0]=(int)(numBaseAfterBreakPoint)+48;
            startPosAfterBreakPoint++;
        }
        for (j=0; j<nt; j++) {
            if (j<typeCutPosBefore) {
                typeBeforeBreakPoint[j]=type[j];
            }
        }
        if (numBaseBeforeBreakPoint>9) {
            for (j=0; j<digitsOfTmp; j++) {
                typeBeforeBreakPoint[typeCutPosBefore+j] = (int)(numBaseBeforeBreakPoint/pow(10,digitsOfTmp-j-1))+48;
                numBaseBeforeBreakPoint -= boost::lexical_cast<unsigned>(floor(numBaseBeforeBreakPoint/pow(10,digitsOfTmp-j-1))*pow(10,digitsOfTmp-j-1));
            }
            typeBeforeBreakPoint[typeCutPosBefore+digitsOfTmp]='\0';

            //typeBeforeBreakPoint[typeCutPosBefore]=(int)(numBaseBeforeBreakPoint/10)+48;
            //typeBeforeBreakPoint[typeCutPosBefore+1]=(int)(numBaseBeforeBreakPoint%10)+48;
            //typeBeforeBreakPoint[typeCutPosBefore+2]='\0';
        } else if (numBaseBeforeBreakPoint>0) {
            typeBeforeBreakPoint[typeCutPosBefore]=(int)(numBaseBeforeBreakPoint)+48;
            typeBeforeBreakPoint[typeCutPosBefore+1]='\0';
        } else {
            typeBeforeBreakPoint[typeCutPosBefore]='\0';
        }
        typeAfterBreakPoint[startPosAfterBreakPoint]='\0';

    }
#ifdef DEBUG_SU
    std::cerr << "Before return: " << breakPointPos << ", enteringBreakPoint: " << enteringBreakPoint << std::endl;
#endif
    if (breakPointPos>0) {
        if (enteringBreakPoint)
            alignScore += ReadScorer::getAlignmentScore(qualityStringBeforeBreakPoint.c_str(),typeBeforeBreakPoint);
        else
            alignScore += ReadScorer::getAlignmentScoreAfterBreakPoint(qualityStringBeforeBreakPoint.c_str(),typeBeforeBreakPoint);
    }
#ifdef DEBUG_SU
    std::cerr << "Before return2: " << breakPointPos << ", enteringBreakPoint: " << enteringBreakPoint << ", lengthOfAlignment: "
              << lengthOfAlignment << std::endl;
#endif
    if (breakPointPos < lengthOfAlignment) {
        if (enteringBreakPoint)
            alignScore += ReadScorer::getAlignmentScoreAfterBreakPoint(qualityStringAfterBreakPoint.c_str(),typeAfterBreakPoint);
        else
            alignScore += ReadScorer::getAlignmentScore(qualityStringAfterBreakPoint.c_str(),typeAfterBreakPoint);
    }

#ifdef DEBUG_SU
    std::cerr << "Finished, final score: " << alignScore << "!" << std::endl;
#endif

    return alignScore;
}

/*****************************************************************************/

int
SequenceUtils::
ReadScorer::
getMostProbableBreakPoint(const Alignment& alignment,
                          bool& poorAlignIsAtStart) const
{
    int breakPointPos = 0;

    const char* qualityString = alignment.getQuality().c_str();
    int lengthOfAlignment = strlen (qualityString);

    const char* type = alignment.getMatchDescriptor().c_str();
    bool noGaps = true;

    for (unsigned j=0; j<alignment.getMatchDescriptor().size(); ++j) {
//			cerr << j <<": type[j]= " << type[j] << endl;
        if (type[j] == '^' ) {
            noGaps = false;
            break;
        }
    }
//cerr <<"in 	getMostProbableBreakPoint, noGaps = " <<noGaps << endl;

    if (!noGaps) {
        int indelPos = -1;
        int exactRun=0;
        int firstMismatch=-1;
        int numLeadingMatches(0);

        for (unsigned int i(0); i<alignment.getMatchDescriptor().size(); i++)
        {
            if (isdigit(alignment.getMatchDescriptor()[i]))
            {
                exactRun*=10;
                exactRun+=
                    (static_cast<int>(alignment.getMatchDescriptor()[i])-48);
            } // if
            else
            {   // found a nondigit character, return leading values
                numLeadingMatches+=exactRun;
                exactRun=0;

                if (isalpha(alignment.getMatchDescriptor()[i]))
                {
                    if (firstMismatch==-1) firstMismatch=numLeadingMatches;
                    numLeadingMatches++;
                }
                else if (alignment.getMatchDescriptor()[i]=='^')
                {
                    if (indelPos==-1) indelPos=numLeadingMatches;
                    return (indelPos);
                }
            } // ~else
        } // ~for

    } //end !noGaps

    double minScore
        = SequenceUtils::getBreakPointHypothesisScore(alignment,
                                                      breakPointPos,
                                                      true);

    for (int i=1; i<lengthOfAlignment; i++) {
        double score
            = SequenceUtils::getBreakPointHypothesisScore(alignment, i,
                                                          true);
        if (score<minScore) {
            minScore = score;
            breakPointPos = i;
            poorAlignIsAtStart = false;
        }

        //cerr  << i << ":" << score << ":";
        score = SequenceUtils::getBreakPointHypothesisScore(alignment, i,
                                                            false);
        if (score<minScore) {
            minScore = score;
            breakPointPos = i;
            poorAlignIsAtStart = true;
        }

        //cerr << score << "\t";
    }

    if (alignment.getMatch().getStrand() == Match::Reverse) {
        //cerr << "breakPoint = " << breakPointPos << endl;
        breakPointPos--;
    }

    //cerr << "\n";
    return breakPointPos;
}

/*****************************************************************************/

double
SequenceUtils::
ReadScorer::
getAlignmentScoreAfterBreakPoint(const char* qualityString,
                                 const char* type) const {
    int posInRead = 0;
    unsigned tmp = 0;
    double alignScore(0.);
    const unsigned nt(strlen(type));
    if (!isValidQualString(qualityString)) {
        return alignScore;
    }

    for (unsigned i=0; i<nt; ++i)
    {
#ifdef DEBUG_SU
        std::cerr << "*** in getAlignmentScore: type[i]= " << type[i] << std::endl;
#endif
        if ((int) type[i] <= 57 && (int) type[i] >= 48) //between 0 and 9
        {
            tmp *= 10;
            tmp += (type[i] - '0');
        }
        else if (type[i] == 'A' || type[i] == 'C' || type[i] == 'G' || type[i] == 'T' || type[i]=='N')
        {
            // match but we penalize for match
            for (unsigned j=0; j<tmp; j++) {
                alignScore +=  _logpcorrectratio[static_cast<int>(qualityString[posInRead+j])];
            }

            posInRead += tmp;
            tmp=0;
            ++posInRead;
        }
        else if (type[i]=='^' || type[i] == '$' ) {
        }
        else {
            std::cerr << "[ERROR]: Unexpected match descriptor in getAlignmentScoreAfterBreakPoint: type " << type << "\n";
            std::cerr << "qualityString " << qualityString << "\n";
            std::cerr << "[ERROR]: Current value being checked: '" << type[i] << "'" << std::endl;
            exit(EXIT_CODE_FAILURE);
        }
    }
    //matches
    if (tmp>0) {
        for (unsigned j=0; j<tmp; j++) {
            alignScore +=  _logpcorrectratio[static_cast<int>(qualityString[posInRead+j])];
        }
    }

    return alignScore;
}

