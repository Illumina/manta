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

#include "SVLocusSetFinder.hh"
#include "starling_common/align_path_bam_util.hh"

#include "boost/foreach.hpp"



SVLocusSetFinder::
SVLocusSetFinder(const ESLOptions& opt)
{
    // pull in insert stats:
    _rss.read(opt.statsFilename.c_str());

    // cache the insert stats we'll be looking up most often:
    BOOST_FOREACH(std::string file, opt.alignmentFilename) {
        const boost::optional<unsigned> index(_rss.getGroupIndex(file));
        assert(index);
        const ReadGroupStats rgs(_rss.getStats(*index));

        _stats.resize(_stats.size()+1);
        _stats.back().min=rgs.fragSize.quantile(opt.breakendEdgeTrimProb);
        _stats.back().max=rgs.fragSize.quantile((1-opt.breakendEdgeTrimProb));
    }
}



void
SVLocusSetFinder::
update(const bam_record& read,
       const unsigned defaultReadGroupIndex)
{
    // shortcut to speed things up:
    if(read.is_proper_pair()) return;

    // start out looking for chimeric reads only:
    if(read.is_chimeric())
    {
        ALIGNPATH::path_t apath;
        bam_cigar_to_apath(read.raw_cigar(),read.n_cigar(),apath);

        const unsigned read_size(apath_read_length(apath));

        // expected breakpoint range is from the end of the read alignment to the (probabalistic) end of the fragment:
        //if(read.is_fwd_strand())
        if(read_size==0) {
            exit(0);
        }
    }
}


