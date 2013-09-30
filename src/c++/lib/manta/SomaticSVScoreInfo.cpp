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
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#include "manta/SomaticSVScoreInfo.hh"
#include "blt_util/log.hh"

#include "boost/foreach.hpp"

#include <iostream>

SVAlignmentInfo::
SVAlignmentInfo(
			const SVCandidate& sv,
			const SVCandidateAssemblyData& assemblyData)
{
	// consider 2-locus events first
	// TODO: to add local assembly later

	// for imprecise SVs, split-read evidence won't be assigned
	if ((assemblyData.isSpanning) &&
		(!sv.isImprecise()))
	{
		log_os << "bestAlignmentIndex=" << assemblyData.bestAlignmentIndex <<"\n";
		if (assemblyData.contigs.size() == 0)
		{
			log_os << "no contigs available in the assembly data.\n";
			return;
		}

		contigSeq = assemblyData.contigs[assemblyData.bestAlignmentIndex].seq;
		log_os << "contigSeq=" << contigSeq <<"\n";
		const JumpAlignmentResult<int>& alignment = assemblyData.spanningAlignments[assemblyData.bestAlignmentIndex];

		// get offsets of breakpoints in the contig
		const unsigned align1Size(apath_read_length(alignment.align1.apath));
		const unsigned insertSize(alignment.jumpInsertSize);
		bp1ContigOffset = align1Size - 1;
		bp2ContigOffset = align1Size + insertSize;
		log_os << "bp1ContigOffset=" << bp1ContigOffset << "\n";
		log_os << "bp2ContigOffset=" << bp2ContigOffset << "\n";
		bp1ContigReversed = false;
		bp2ContigReversed = false;
		if (sv.bp1.state == sv.bp2.state)
		{
			if (sv.bp1.state == SVBreakendState::RIGHT_OPEN)
				bp2ContigReversed = true;
			else
				bp1ContigReversed = true;
		}

		// get reference regions
		const reference_contig_segment& bp1Ref = assemblyData.bp1ref;
		const reference_contig_segment& bp2Ref = assemblyData.bp2ref;
		bp1RefSeq = bp1Ref.seq();
		bp2RefSeq = bp2Ref.seq();
		log_os << "bp1RefSeq=" << bp1RefSeq << "\n";
		log_os << "bp2RefSeq=" << bp2RefSeq << "\n";
		// get offsets of breakpoints in the reference regions
		bp1RefOffset = sv.bp1.interval.range.begin_pos() - bp1Ref.get_offset();
		bp2RefOffset = sv.bp2.interval.range.begin_pos() - bp2Ref.get_offset();
	}
}


std::ostream&
operator<<(std::ostream& os, const SVSampleInfo& si)
{
    os << "SVSampleInfo bp1SpanReads=" << si.bp1SpanReads << " bp2SpanReads=" << si.bp2SpanReads << " spanPairs=" << si.spanPairs << std::endl;
    return os;
}


std::ostream&
operator<<(std::ostream& os, const SomaticSVScoreInfo& ssi)
{
    os << "SomaticSVScoreInfo bp1MaxDepth=" << ssi.bp1MaxDepth << " bp2MaxDepth=" << ssi.bp2MaxDepth << " somaticScore=" << ssi.somaticScore << std::endl;
    os << "Tumor sample info " << ssi.tumor;
    os << "Normal sample info " << ssi.normal;
    BOOST_FOREACH(const std::string& filter, ssi.filters)
    {
        os << " " << filter;
    }
    return os;
}


