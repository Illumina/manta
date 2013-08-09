
#pragma once

#include <iomanip>
#include <sstream>

#include <boost/utility.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

#include <cstring>

//#define DEBUG_SU


typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;

const int phredScoreOffset = 33;

// @class SequenceUtils
class SequenceUtils : private boost::noncopyable {
public:
    // compute the total read quality from its ASCII representation
    static
    int totalBaseQualityForRead(const char* readQuality) {
        int sum(0);
        const unsigned read_size(strlen(readQuality));
        for (unsigned int i(0); i<read_size; i++) {
            sum+=((int)readQuality[i])-phredScoreOffset;
        } // ~for
        return sum;
    } // ~totalBaseQualityForRead( const char* readQuality)


    static
    std::string getSampleName(const Alignment& alignment) {
        char buf[256];
        sprintf(buf, "%s %d %d", alignment.getSpot().getTile().getMachineName().c_str(), alignment.getSpot().getTile().getRunNumber(), alignment.getSpot().getTile().getLaneNumber());
        return (std::string)buf;
    } // ~getSampleName

    static
    std::string getLaneKey(const Alignment& alignment) {
        std::ostringstream ostrm;
        ostrm << alignment.getSpot().getTile().getMachineName() << "_"
              // << std::setfill('0') << std::setw(4)
              << alignment.getSpot().getTile().getRunNumber() << ":"
              << alignment.getSpot().getTile().getLaneNumber();
        return ostrm.str();
    }

    static
    std::string getReadId(const Alignment& alignment, bool ofPartner=false) {
        std::ostringstream ostrm;
        ostrm << getLaneKey(alignment) << ":"
              << alignment.getTileNumber() << ":"
              << alignment.getX() << ":"
              << alignment.getY();

        const std::string indexStr(alignment.getIndex());

        if (!indexStr.empty()) {
            ostrm << "#" << indexStr;
        }

        unsigned int readNumber(alignment.getReadNumber());

        if (readNumber != 0) {
            if (ofPartner) {
                readNumber = (readNumber == 1) ? 2 : 1;
            }

            ostrm << "/" << readNumber;
        }

        return ostrm.str();
    }

    static
    std::string getReadPartnerId(const std::string& readId) {
        boost::char_separator<char> separ("/");
        Tokenizer tokenizer(readId, separ);

        std::vector<std::string>
        partsStrVec(tokenizer.begin(), tokenizer.end());

        if (partsStrVec.size() < 2) { // single read readID -> just return it
            return readId;
        }

        // Assume 1 <-> 2
        const unsigned int
        readNum(boost::lexical_cast<unsigned int>(partsStrVec[1]));
        const std::string partnerReadId(partsStrVec[0] + "/"
                                        + boost::lexical_cast<std::string>(
                                            3 - readNum));
        return partnerReadId;
    }

    static
    bool getReadPairInfo(const Alignment& alignment,
                         const unsigned readLen,
                         const unsigned partnerReadLen,
                         bool& isChimeric,
                         long& leftPos,
                         unsigned int& insertSize,
                         RelOrient& relOrient);

    static
    double convertQualToProbCorrect(int qv) {
        return 1.0/(1+std::pow(10, (-qv/10.0)));
    }

    static
    double convertPhredToProbCorrect(int qv) {
        return (1.0-std::pow(10, (-qv/10.0)));
    }

    static
    double convertLogOddsToProbError(int qv) {
        return 1./(1.+std::pow(10.,(static_cast<double>(qv)/10.)));
    }

    static
    double convertPhredToProbError(int qv) {
        return std::min(1.,std::pow(10.,(-static_cast<double>(qv)/10.)));
    }

    static
    bool
    isValidQualString(const char* qual);

    static
    double
    getSemiAlignedReadMetric(const Alignment& alignment);

    static
    double
    getSemiAlignedReadMetric(const BamTools::BamAlignment& alignment);

    static
    double
    getBreakPointHypothesisScore(const Alignment& alignment,
                                 const int breakPointPos,
                                 const bool enteringBreakPoint);

    static
    int getMostProbableBreakPoint(const Alignment& alignment,
                                  bool& poorAlignIsAtStart);

private:
    struct ReadScorer : private boost::noncopyable {

        /** Instance getter
         *
        */
        static const ReadScorer& get() {
            static const ReadScorer rs;
            return rs;
        }

        bool
        isValidQualString(const char* qual) const;

        double getSemiAlignedReadMetric(const Alignment& alignment) const;

        // same check for BamAlignment, it is the responsibility of the caller to ensure that a matchdescriptor is set
        double getSemiAlignedReadMetric( const BamTools::BamAlignment& alignment ) const;

        double getAlignmentScoreAfterBreakPoint(const char* qualityString, const char* type) const;

        double getBreakPointHypothesisScore(const Alignment& alignment,
                                            const unsigned breakPointPos,
                                            const bool enteringBreakPoint) const;

        int getMostProbableBreakPoint(const Alignment& alignment,
                                      bool& poorAlignIsAtStart) const;

    private:
        explicit
        ReadScorer();
        ~ReadScorer() {}

        double
        getAlignmentScore(const char* qualityString,
                          const char* type) const;

        enum { MAX_Q = 128 };
        const int _qmin;
        double _logpcorrectratio[MAX_Q];
    };


    SequenceUtils() {}
    ~SequenceUtils() {}
};


