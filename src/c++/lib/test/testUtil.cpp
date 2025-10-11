//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2025 Illumina, Inc.
//
// This program is licensed under the terms of the Polyform strict license
//
// ***As far as the law allows, the software comes as is, without
// any warranty or condition, and the licensor will not be liable
// to you for any damages arising out of these terms or the use
// or nature of the software, under any kind of legal claim.***
//
// You should have received a copy of the PolyForm Strict License 1.0.0
// along with this program.  If not, see <https://polyformproject.org/licenses/strict/1.0.0>.
//
//

/// \file
/// \author Trevor Ramsay
///

#include "testUtil.hpp"

#include "blt_util/string_util.hpp"

#include "boost/filesystem.hpp"

#include <fstream>

std::string getNewTempFile()
{
  using namespace boost::filesystem;

  const path tempDir(temp_directory_path());
  while (true) {
    path tempFile = tempDir / unique_path();
    if (exists(tempFile)) continue;
    return tempFile.string();
  }
}

std::string getValueFromTSVKeyValFile(const std::string& tsvFile, const std::string& key)
{
  std::ifstream is(tsvFile);

  static const char        delimiter('\t');
  std::string              line;
  std::vector<std::string> words;
  while (std::getline(is, line)) {
    if (line.find(key) == std::string::npos) continue;
    split_string(line, delimiter, words);
    if (words.size() < 2) continue;
    if (words[0] == key) return words[1];
  }
  return "";
}

#if 0
#include "blt_util/align_path.hpp"
#include "htsapi/align_path_bam_util.hpp"

#include <memory>

#include "htsapi/SimpleAlignment_bam_util.hpp"
#include "htsapi/bam_dumper.hpp"
#include "htsapi/bam_record.hpp"
#include "htsapi/bam_record_util.hpp"
#include "htsapi/sam_util.hpp"
#include "manta/SVLocusScanner.hpp"

#include "boost/algorithm/string.hpp"
#include "boost/scoped_array.hpp"
#include "boost/test/unit_test.hpp"

#include <fstream>
#include <iostream>

inline
std::istream&
operator>>(std::istream& is, line& l)
{
    std::getline(is, l);
    return is;
}

inline
std::string
bamRecordToString(
    std::string readName,
    bam_record record,
    bam_header_info header)
{
    auto ConcatWithTab = [](std::string input)
    {
        return input + "\t";
    };

    std::string strOut = "";

    //QNAME FLAG RNAME POS MAPQ CIGAR MRNM/RNEXT MPOS/PNEXT ISIZE SEQ QUAL TAGs

    // Get chrom name from header info.
    std::string readTargetId = header.chrom_data[record.get_data()->core.tid].label;
    std::string readMateTargetId = header.chrom_data[record.get_data()->core.mtid].label;
    if ( readTargetId == readMateTargetId) readMateTargetId = "=";

    strOut.append(ConcatWithTab(readName));
    strOut.append(ConcatWithTab(std::to_string(record.get_data()->core.flag)));
    strOut.append(ConcatWithTab(readTargetId));
    strOut.append(ConcatWithTab(std::to_string(record.pos())));
    strOut.append(ConcatWithTab(std::to_string(record.map_qual())));

    ALIGNPATH::path_t apath;
    bam_cigar_to_apath(record.raw_cigar(), record.n_cigar(), apath);
    std::string cigarString = "";
    apath_to_cigar(apath, cigarString);
    strOut.append(ConcatWithTab(cigarString));

    strOut.append(ConcatWithTab(readMateTargetId));
    strOut.append(ConcatWithTab(std::to_string(record.mate_pos())));
    strOut.append(ConcatWithTab(std::to_string(record.get_data()->core.isize)));
    if ( record.get_bam_read().size() == 0)
    {
        strOut.append(ConcatWithTab("*"));
        strOut.append(ConcatWithTab("*"));
    }
    else
    {
        strOut.append(ConcatWithTab(record.get_bam_read().get_string()));

        for ( unsigned i = 0; i < record.get_bam_read().size(); i++ )
        {
            std::string  s(1, '0' + record.qual()[i]);
            strOut.append(s);
        }
        strOut.append("\t");
    }

    std::string auxInfo = "NH:i:1\tNM:i:0\tRG:Z:1";

    strOut.append(auxInfo);

    static const char rgTag[] = {'S','A'};
    const char* pTag(record.get_string_tag(rgTag));

    if ( NULL != pTag)
    {
        strOut.append("\t");
        std::string stringTag(( NULL == pTag) ? "" : pTag);
        strOut.append(stringTag);
    }

    return strOut;
}

#endif
