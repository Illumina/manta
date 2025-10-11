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
/// \author Chris Saunders
///

#include "svgraph/SVLocusSampleCounts.hpp"
#include "blt_util/io_util.hpp"

#include <iomanip>
#include <iostream>

static void writeLine(std::ostream& os, const char* label, const double val)
{
  static const char sep('\t');

  os << std::fixed;
  os << label << sep;
  os << std::setprecision(0);
  os << val << sep;
  os << "N/A" << '\n';
}

static void writeLine(std::ostream& os, const char* label, const double val, const double total)
{
  static const char sep('\t');

  os << std::fixed;
  os << label << sep;
  os << std::setprecision(0);
  os << val << sep;
  os << std::setprecision(4);
  os << val / total << '\n';
}

void SampleReadInputCounts::write(std::ostream& os) const
{
  const double dtotal(total());
  StreamScoper ss(os);
  writeLine(os, "MinMapqFiltered", minMapq, dtotal);
  writeLine(os, "NotFiltered", evidenceCount.total, dtotal);
  writeLine(os, "NotFilteredAndIgnored", evidenceCount.ignored, dtotal);
  writeLine(os, "NotFilteredAndAnomalousPair", evidenceCount.anom, dtotal);
  writeLine(os, "NotFilteredAndAnomalousPairRemotes", evidenceCount.remoteRecoveryCandidates, dtotal);
  writeLine(os, "NotFilteredAndSplitRead", evidenceCount.split, dtotal);
  writeLine(os, "NotFilteredAndSplitReadInAnomalousPair", evidenceCount.anomAndSplit, dtotal);
  writeLine(os, "NotFilteredAndSplitReadSupplementarySegments", evidenceCount.splitSupplementarySegment);
  writeLine(os, "NotFilteredAndLargeIndel", evidenceCount.indel, dtotal);
  writeLine(os, "NotFilteredAndSemiAligned", evidenceCount.assm, dtotal);
}

void SampleEvidenceCounts::write(std::ostream& os) const
{
  static const char sep('\t');

  double total(0);
  for (unsigned i(0); i < SVEvidenceType::SIZE; ++i) {
    total += eType[i];
  }

  StreamScoper ss(os);
  os << std::fixed << std::setprecision(4);
  for (unsigned i(0); i < SVEvidenceType::SIZE; ++i) {
    os << "EvidenceType_" << SVEvidenceType::label(i) << sep << eType[i] << sep << eType[i] / total << '\n';
  }
  os << "ClosePairs" << sep << closeCount << '\n';
}

void SampleReadCounts::write(std::ostream& os, const char* label) const
{
  os << "\n[" << label << "]\n";
  os << "Source\t" << sampleSource << "\n";
  input.write(os);
  evidence.write(os);
}

void AllSampleReadCounts::write(std::ostream& os, const std::vector<std::string>& sampleLabels) const
{
  assert(size() == sampleLabels.size());
  const unsigned s(size());
  for (unsigned i(0); i < s; ++i) {
    getSampleCounts(i).write(os, sampleLabels[i].c_str());
  }
}
