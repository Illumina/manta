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

#include "EstimateSVLoci.hh"
#include "ESLOptions.hh"
#include "SVLocusSetFinder.hh"

#include "blt_util/input_stream_handler.hh"
#include "blt_util/log.hh"
#include "common/OutStream.hh"
#include "manta/ReadGroupStatsSet.hh"

#include "boost/foreach.hpp"
#include "boost/shared_ptr.hpp"

#include <iostream>
#include <vector>



static
void
runESL(const ESLOptions& opt)
{

    {
        // early test that we have permission to write to output file
        OutStream outs(opt.outputFilename);
    }

    SVLocusSetFinder locusFinder(opt);

    typedef boost::shared_ptr<bam_streamer> stream_ptr;
    std::vector<stream_ptr> bam_streams;

    // setup all data for main alignment loop:
    BOOST_FOREACH(const std::string& afile, opt.alignmentFilename)
    {
        stream_ptr tmp(new bam_streamer(afile.c_str(),opt.region.c_str()));
        bam_streams.push_back(tmp);
    }

    // TODO check header compatibility between all open bam streams
    const unsigned n_inputs(bam_streams.size());

    // assume headers compatible after this point....

    if (n_inputs)
    {
        locusFinder.setBamHeader(*(bam_streams[0]->get_header()));
    }

    input_stream_data sdata;
    for (unsigned i(0); i<n_inputs; ++i)
    {
        sdata.register_reads(*bam_streams[i],i);
    }

    // loop through alignments:
    input_stream_handler sinput(sdata);
    while (sinput.next())
    {
        const input_record_info current(sinput.get_current());

        if       (current.itype != INPUT_TYPE::READ)
        {
            log_os << "ERROR: invalid input condition.\n";
            exit(EXIT_FAILURE);
        }

        const bam_streamer& read_stream(*bam_streams[current.sample_no]);
        const bam_record& read(*(read_stream.get_record_ptr()));

        locusFinder.update(read,current.sample_no);
    }

    // debug output
    locusFinder.getLocusSet().save(opt.outputFilename.c_str());
}



void
EstimateSVLoci::
runInternal(int argc, char* argv[]) const
{

    ESLOptions opt;

    parseESLOptions(*this,argc,argv,opt);
    runESL(opt);
}
