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
// <https://github.com/downloads/sequencing/licenses/>.
//

///
/// \author Ole Schulz-Trieglaff
///

#include "AssembleSVBreakend.hh"

#include "applications/GenerateSVCandidates/GSCOptions.hh"
#include "applications/AssembleSVBreakend/ASBOptions.hh"

#include "manta/SVLocusAssembler.hh"
#include "manta/SVReferenceUtil.hh"
#include "assembly/AssembledContig.hh"

#include "blt_util/bam_header_util.hh"
#include "blt_util/input_stream_handler.hh"
#include "blt_util/log.hh"
#include "blt_util/io_util.hh"

#include "alignment/AlignmentUtil.hh"
#include "alignment/GlobalAligner.hh"
#include "alignment/GlobalJumpAligner.hh"

#include "common/OutStream.hh"

#include "boost/foreach.hpp"
#include "boost/shared_ptr.hpp"

#include <fstream>
#include <iostream>


void
analyzeAlignment(const JumpAlignmentResult<int>& res,
                 const pos_t refOffset1,
                 const pos_t refOffset2)
{
    //onst unsigned minAlignContext(4);

    std::cout << "alignment score = " << res.score << std::endl;
    //std::cout << "alignment jumpInsertSize = " << res.jumpInsertSize << std::endl;
    std::string cigar1;
    apath_to_cigar(res.align1.apath,cigar1);
    std::cout << "align1 start = " << res.align1.beginPos << std::endl;
    std::cout << "align1 cigar = " << cigar1 << std::endl;
    //std::cout << "align1 cigar prefix aligned " << (hasAlignedPrefix(res.align1,minAlignContext) > 0 ? "Yes" : "No") << std::endl;
    //std::cout << "align1 cigar suffix aligned " << (hasAlignedSuffix(res.align1,minAlignContext) > 0 ? "Yes" : "No") << std::endl;

    std::string cigar2;
    apath_to_cigar(res.align2.apath,cigar2);
    std::cout << "align2 start = " << res.align2.beginPos << std::endl;
    std::cout << "align2 cigar = " << cigar2 << std::endl;
    //std::cout << "align2 cigar prefix aligned " << (hasAlignedPrefix(res.align2,minAlignContext) > 0 ? "Yes" : "No") << std::endl;
    //std::cout << "align2 cigar suffix aligned " << (hasAlignedSuffix(res.align2,minAlignContext) > 0 ? "Yes" : "No") << std::endl;

    std::cout << "breakpoint estimate 1 " << estimateBreakPointPos(res.align1, refOffset1) << std::endl;
    std::cout << "breakpoint estimate 2 " << estimateBreakPointPos(res.align2, refOffset2) << std::endl;

    std::cout << "Alignment 1 length " << apath_read_length(res.align1.apath) << std::endl;
    std::cout << "Alignment 2 length " << apath_read_length(res.align2.apath) << std::endl;

    if (! (res.align1.isAligned() && res.align2.isAligned()) )
    {
        std::cout << "At least one not aligned." << std::endl;
    }

//    if (!isConsistentAlignment(res,minAlignContext))
//    {
//        std::cout << "Alignments not consistent." << std::endl;
//    }
}

static
void
runASB(const ASBOptions& opt)
{
    ///
    typedef boost::shared_ptr<bam_streamer> stream_ptr;
    std::vector<stream_ptr> bam_streams_bkpt1;
    std::vector<stream_ptr> bam_streams_bkpt2;

    // setup all data for main alignment loop:
    BOOST_FOREACH(const std::string& afile, opt.alignFileOpt.alignmentFilename)
    {
        stream_ptr tmp(new bam_streamer(afile.c_str(),opt.breakend1.c_str()));
        bam_streams_bkpt1.push_back(tmp);
    }
    BOOST_FOREACH(const std::string& afile, opt.alignFileOpt.alignmentFilename)
    {
        stream_ptr tmp(new bam_streamer(afile.c_str(),opt.breakend2.c_str()));
        bam_streams_bkpt2.push_back(tmp);
    }

    // TODO check header compatibility between all open bam streams
    const unsigned n_inputs(bam_streams_bkpt1.size());
    // assume headers compatible after this point....
    assert(0 != n_inputs);

    std::cout << "Assembling region 1 " << opt.breakend1 << std::endl;
    std::cout << "Assembling region 2 " << opt.breakend2 << std::endl;

    const bam_header_t& header(*(bam_streams_bkpt1[0]->get_header()));
    const bam_header_info bamHeader(header);

    int32_t tid1(0), beginPos1(0), endPos1(0);
    int32_t tid2(0), beginPos2(0), endPos2(0);
    parse_bam_region(bamHeader,opt.breakend1,tid1,beginPos1,endPos1);
    parse_bam_region(bamHeader,opt.breakend2,tid2,beginPos2,endPos2);

    const GenomeInterval breakendRegion1(tid1,beginPos1,endPos1);
    SVBreakend bp1;
    bp1.interval = breakendRegion1;
    const GenomeInterval breakendRegion2(tid2,beginPos2,endPos2);
    SVBreakend bp2;
    bp2.interval = breakendRegion2;

    std::cout << "Translating into " << breakendRegion1 << " and " << breakendRegion2 << std::endl;
    std::cout << "Translating into " << bp1 << " and " << bp2 << std::endl;

    SmallAssemblerOptions assembleOpt;
    SVLocusAssembler svla(opt.scanOpt, assembleOpt, opt.alignFileOpt, opt.statsFilename);
    Assembly a;
    svla.assembleSVBreakends(bp1,bp2, false, false, a);
    std::cout << "Assembled " << a.size() << " contig(s)." << std::endl;

    // how much additional reference sequence should we extract from around
    // each side of the breakend region?
    static const pos_t extraRefEdgeSize(300);

    reference_contig_segment bpref1;
    getIntervalReferenceSegment(opt.referenceFilename, bamHeader, extraRefEdgeSize, bp1.interval, bpref1);
    const std::string bpRefStr1(bpref1.seq());

    reference_contig_segment bpref2;
    getIntervalReferenceSegment(opt.referenceFilename, bamHeader, extraRefEdgeSize, bp2.interval, bpref2);
    const std::string bpRefStr2(bpref2.seq());

    const pos_t refOffset1(std::max(0, (bp1.interval.range.begin_pos()-extraRefEdgeSize)));
    const pos_t refOffset2(std::max(0, (bp2.interval.range.begin_pos()-extraRefEdgeSize)));

    int jumpScore(-3);
    GlobalJumpAligner<int> aligner(AlignmentScores<int>(1,-2,-6,-0,-2),jumpScore);

    std::ofstream os(opt.contigOutfile.c_str());

    unsigned minContigLength(75);

    unsigned numTotal(0);
    //unsigned numConsistentAndAligned(0);
    unsigned numTooShort(0);


    for (Assembly::const_iterator ct = a.begin();
         ct != a.end();
         ++ct)
    {
        if (ct->seq.size() < minContigLength)
        {
            ++numTooShort;
            continue;
        }

        std::cout << "ctg : " << ct->seq << std::endl;
        std::cout << "contig length " << ct->seq.size() << std::endl;
        std::cout << "seed read count " << ct->seedReadCount << std::endl;

        os << ">contig_" << numTotal << std::endl;
        os << ct->seq << std::endl;

        ++numTotal;

        std::cout << "FORWARD" << std::endl;
        JumpAlignmentResult<int> forwardRes;
        aligner.align(ct->seq.begin(),ct->seq.end(),
                      bpRefStr1.begin(),bpRefStr1.end(),
                      bpRefStr2.begin(),bpRefStr2.end(),
                      forwardRes);
        analyzeAlignment(forwardRes,refOffset1,refOffset2);

        std::cout << "BACKWARD" << std::endl;
        JumpAlignmentResult<int> backwardRes;
        aligner.align(ct->seq.begin(),ct->seq.end(),
                      bpRefStr2.begin(),bpRefStr2.end(),
                      bpRefStr1.begin(),bpRefStr1.end(),
                      backwardRes);
        analyzeAlignment(backwardRes,refOffset1,refOffset2);

        //++numConsistentAndAligned;
    }
    std::cout << "\nAssembled " << numTotal << " contigs." << std::endl;
    //std::cout << "numConsistentAndAligned = " << numConsistentAndAligned << std::endl;
    std::cout << "numTooShort = " << numTooShort << std::endl;
    os.close();
}


void
AssembleSVBreakend::
runInternal(int argc, char* argv[]) const
{
    ASBOptions opt;
    parseASBOptions(*this,argc,argv,opt);
    runASB(opt);
}


