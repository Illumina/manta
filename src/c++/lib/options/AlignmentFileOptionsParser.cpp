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

#include "options/AlignmentFileOptionsParser.hpp"
#include "options/optionsUtil.hpp"

#include <set>

typedef std::vector<std::string> files_t;

boost::program_options::options_description getOptionsDescription(AlignmentFileOptions& /*opt*/)
{
  namespace po = boost::program_options;
  po::options_description desc("alignment-files");
  //clang-format off
  desc.add_options()(
      "align-file",
      po::value<files_t>(),
      "alignment file in BAM or CRAM format (may be specified multiple times, assumed to be non-tumor if tumor file(s) provided)")(
      "tumor-align-file",
      po::value<files_t>(),
      "tumor sample alignment file in BAM or CRAM format (may be specified multiple times)");
  // clang-format on

  return desc;
}

bool parseOptions(
    const boost::program_options::variables_map& vm, AlignmentFileOptions& opt, std::string& errorMsg)
{
  errorMsg.clear();

  // paste together tumor and normal:
  {
    files_t normal;
    files_t tumor;
    if (vm.count("align-file")) {
      normal = (boost::any_cast<files_t>(vm["align-file"].value()));
    }
    if (vm.count("tumor-align-file")) {
      tumor = (boost::any_cast<files_t>(vm["tumor-align-file"].value()));
    }
    opt.alignmentFilenames = normal;
    opt.alignmentFilenames.insert(opt.alignmentFilenames.end(), tumor.begin(), tumor.end());
    opt.isAlignmentTumor.clear();
    opt.isAlignmentTumor.resize(normal.size(), false);
    opt.isAlignmentTumor.resize(opt.alignmentFilenames.size(), true);
  }

  if (opt.alignmentFilenames.empty()) {
    errorMsg = "Must specify at least one input alignment file";
  } else {
    // check that alignment files exist, and names do not repeat
    std::set<std::string> nameCheck;
    for (std::string& afile : opt.alignmentFilenames) {
      if (checkAndStandardizeRequiredInputFilePath(afile, "alignment file", errorMsg)) break;
      if (nameCheck.count(afile)) {
        std::ostringstream oss;
        oss << "Repeated alignment filename: " << afile << "\n";
        errorMsg = oss.str();
        break;
      }
      nameCheck.insert(afile);
    }
  }

  return (!errorMsg.empty());
}
