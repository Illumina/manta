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
#include "manta/AssembledContig.hh"

#include "blt_util/bam_header_util.hh"
#include "blt_util/input_stream_handler.hh"
#include "blt_util/log.hh"
#include "blt_util/io_util.hh"

#include "common/OutStream.hh"

#include "boost/foreach.hpp"
#include "boost/shared_ptr.hpp"

#include <fstream>
#include <iostream>


static
void
runASB(const ASBOptions& opt)
{
	///
    typedef boost::shared_ptr<bam_streamer> stream_ptr;
    std::vector<stream_ptr> bam_streams;

    // setup all data for main alignment loop:
    BOOST_FOREACH(const std::string& afile, opt.alignmentFilename)
    {
        stream_ptr tmp(new bam_streamer(afile.c_str(),opt.breakend.c_str()));
        bam_streams.push_back(tmp);
    }

    // TODO check header compatibility between all open bam streams
    const unsigned n_inputs(bam_streams.size());

    std::cout << "Assembling region " << opt.breakend << std::endl;

    // assume headers compatible after this point....
    assert(0 != n_inputs);

    const bam_header_t& header(*(bam_streams[0]->get_header()));
    const bam_header_info bamHeader(header);

    int32_t tid(0), beginPos(0), endPos(0);
    parse_bam_region(bamHeader,opt.breakend,tid,beginPos,endPos);

    const GenomeInterval breakendRegion(tid,beginPos,endPos);
    SVBreakend bp;
    bp.interval = breakendRegion;

    std::cout << "Translating into " << breakendRegion << std::endl;
    std::cout << "Translating into " << bp << std::endl;

    SVLocusAssembler svla(opt);
    Assembly a;
    svla.assembleSVBreakend(bp,a);
    
    std::cout << "Assembled " << a.size() << " contig(s)." << std::endl;
    std::ofstream os(opt.contigOutfile.c_str());
    unsigned n(1);
    for (Assembly::const_iterator ct = a.begin();
         ct != a.end();
         ++ct)
    {
        //std::cout << ct->seq << std::endl;
        //std::cout << "seed read count " << ct->seedReadCount << std::endl;
        os << ">contig_" << n << std::endl;
        os << ct->seq << std::endl;
        ++n;
    }
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


