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
/// \author Chris Saunders
///

#include "GenerateSVCandidates.hh"
#include "GSCOptions.hh"
#include "EdgeRetriever.hh"

#include "blt_util/input_stream_handler.hh"
#include "blt_util/bam_streamer.hh"
#include "common/OutStream.hh"
#include "manta/ReadGroupStatsSet.hh"
#include "manta/SVLocusSet.hh"

#include "boost/shared_ptr.hpp"
#include "boost/foreach.hpp"

#include <iostream>



static
void
runGSC(const GSCOptions& opt)
{
    {
        // early test that we have permission to write to output file
        OutStream outs(opt.outputFilename);
    }

    // load in set:
    SVLocusSet set;
    set.load(opt.graphFilename.c_str());
    const SVLocusSet& cset(set);

    EdgeRetriever edger(opt.binIndex,opt.binCount,cset);

    while(edger.next())
    {
        const EdgeInfo& edge=edger.getEdge();

        ///BOGUS:
        if(edge.locusIndex==0) return;
    }





    typedef boost::shared_ptr<bam_streamer> stream_ptr;
    std::vector<stream_ptr> bam_streams;

    // setup all data for main analysis loop:
    BOOST_FOREACH(const std::string& afile, opt.alignmentFilename)
    {
        stream_ptr tmp(new bam_streamer(afile.c_str()));
        bam_streams.push_back(tmp);
    }

    // TODO check header compatibility between all open bam streams
    const unsigned n_inputs(bam_streams.size());

    // assume headers compatible after this point....

    assert(0 != n_inputs);

    input_stream_data sdata;
    for (unsigned i(0); i<n_inputs; ++i)
    {
        sdata.register_reads(*bam_streams[i],i);
    }

    input_stream_handler sinput(sdata);

}



void
GenerateSVCandidates::
runInternal(int argc, char* argv[]) const
{

    GSCOptions opt;

    parseGSCOptions(*this,argc,argv,opt);
    runGSC(opt);
}
