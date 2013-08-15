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

#pragma once

#include "alignment/Alignment.hh"
#include "alignment/GlobalJumpAligner.hh"

#include <cmath>

#include <boost/utility.hpp>


const double invlog10(1./(std::log(10.)));
const int phredScoreOffset = 33;

inline
unsigned
alignEnd(const Alignment& align)
{
    return (align.alignStart + ALIGNPATH::apath_ref_length(align.apath));
}

/// tests if prefix of aligned sequence matches target, returns length of alignment (zero if no match)
unsigned
hasAlignedPrefix(const Alignment& al, const unsigned minAlignContext = 0);


/// tests if suffix of aligned sequence matches target, returns length of alignment (zero if no match)
unsigned
hasAlignedSuffix(const Alignment& al, const unsigned minAlignContext = 0);


bool
bothEndsAligned(const Alignment& al, const unsigned minAlignContext = 0);


/// check a jump alignment for consistency (only one end aligning)
/// FIXME: not used, need to think what makes an alignment consistent
/// (how about : total number of matches shouldn't exceed sequence length?)
bool
isConsistentAlignment(const JumpAlignmentResult<int>& res, const unsigned minAlignContext = 0);


int
estimateBreakPointPos(const Alignment& al, const unsigned refOffset);



struct ReadScorer : private boost::noncopyable {

	/** Instance getter
	 *
	*/
	static const ReadScorer& get() {
		static const ReadScorer rs;
		return rs;
	}


	double
	getSemiAlignedMetric(const std::string& matchDescriptor, const std::string& qualities);

private:
	explicit
	ReadScorer();
	~ReadScorer() {}

	enum { MAX_Q = 128 };
	const int _qmin;
	double _logpcorrectratio[MAX_Q];
};
