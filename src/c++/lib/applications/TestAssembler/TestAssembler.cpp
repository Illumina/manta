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

#include "TestAssembler.hpp"

static void runTestAssembler(const TestAssemblerOptions& opt)
{
  // check that we have write permission on the output file early:
  {
    OutStream outs(opt.outputFilename);
  }

  const ReadScannerOptions scanOpt;
  const AssemblerOptions   asmOpt;

  AssemblyReadInput reads;
  for (const std::string& alignmentFilename : opt.alignFileOpt.alignmentFilenames) {
    log_os << "[INFO] Extracting reads from file: '" << alignmentFilename << "'\n";

    extractAssemblyReadsFromBam(scanOpt, asmOpt, opt.referenceFilename, alignmentFilename, reads);
  }

  AssemblyReadOutput readInfo;
  Assembly           contigs;

  log_os << "[INFO] Assmbling read input.\n";
  runIterativeAssembler(asmOpt, reads, readInfo, contigs);

  OutStream     outs(opt.outputFilename);
  std::ostream& os(outs.getStream());

  const unsigned contigCount(contigs.size());
  log_os << "[INFO] Assembly complete. Contig count: " << contigCount << "\n";

  for (unsigned contigIndex(0); contigIndex < contigCount; ++contigIndex) {
    os << ">Contig" << contigIndex << "\n";
    os << contigs[contigIndex].seq << "\n";
  }
}

void TestAssembler::runInternal(int argc, char* argv[]) const
{
  TestAssemblerOptions opt;

  parseTestAssemblerOptions(*this, argc, argv, opt);
  runTestAssembler(opt);
}
