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
/// \author Xiaoyu Chen
///

#pragma once

#include <string>
#include <iostream>

struct SRAlignmentInfo
{
    SRAlignmentInfo():
        _leftSize(0),
        _rightSize(0),
        _leftMismatches(0),
        _rightMismatches(0),
        _alignScore(0)
    {}

    unsigned get_leftSize() const
    {
        return _leftSize;
    }

    unsigned get_rightSize() const
    {
        return _rightSize;
    }

    unsigned get_leftMismatches() const
    {
        return _leftMismatches;
    }

    unsigned get_rightMismatches() const
    {
        return _rightMismatches;
    }

    unsigned get_alignScore()const
    {
        return _alignScore;
    }

    void set_sizes(unsigned lfs, unsigned rfs)
    {
        _leftSize = lfs;
        _rightSize = rfs;
    }

    void set_mismatches(unsigned lmm, unsigned rmm)
    {
        _leftMismatches = lmm;
        _rightMismatches = rmm;
    }

    void set_score(unsigned asc)
    {
        _alignScore = asc;
    }


    void update_alignInfo(
        const SRAlignmentInfo& info)
    {
        _leftSize = info._leftSize;
        _rightSize = info._rightSize;
        _leftMismatches = info._leftMismatches;
        _rightMismatches = info._rightMismatches;
        _alignScore = info._alignScore;
    }

private:
    unsigned _leftSize;
    unsigned _rightSize;
    unsigned _leftMismatches;
    unsigned _rightMismatches;
    unsigned _alignScore;

};


struct splitReadAlignment
{
    splitReadAlignment():
        _hasEvidence(false),
        _evidence(0)
    {}

    float get_evidence() const
    {
        return _evidence;
    }

    bool has_evidence() const
    {
        return _hasEvidence;
    }

    const SRAlignmentInfo& get_alignment() const
    {
        return _alignment;
    }

    unsigned
    calculateAlignScore(
        const std::string& querySeq,
        const std::string& scanWindow);

    void
    align(const std::string& querySeq,
          const std::string& targetSeq,
          const unsigned bpOffset);


private:

    void set_evidence();

    SRAlignmentInfo _alignment;
    bool _hasEvidence;
    float _evidence;

};

std::ostream&
operator<<(std::ostream& os, const SRAlignmentInfo& info);

std::ostream&
operator<<(std::ostream& os, const splitReadAlignment& srAlign);
