//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2019 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
