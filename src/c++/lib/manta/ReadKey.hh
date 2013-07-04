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

#pragma once

#include "blt_util/bam_record.hh"

#include <iosfwd>
#include <string>


/// information required to uniquely identify a bam read:
///
struct ReadKey
{

    ReadKey(const bam_record& br)
      : _qname(br.qname())
      , _readNo(br.read_no())
    {}

    ReadKey(const std::string& qname,
            const int read_no)
      : _qname(qname)
      , _readNo(read_no)
    {}

    int
    readNo() const { return _readNo; }

    const char*
    qname() const { return _qname.c_str(); }

    bool operator<(const ReadKey& rhs) const {
        if (readNo() < rhs.readNo())
            return true;
        if (readNo() == rhs.readNo()) {
            return (strcmp(qname(), rhs.qname()) < 0);
        }
        return false;
    }

    bool operator==(const ReadKey& rhs) const {
        return ((readNo() == rhs.readNo()) and ((0 == strcmp(qname(), rhs.qname()))));
    }

private:
    std::string _qname;
    int _readNo;
};


std::ostream& operator<<(std::ostream& , const ReadKey& );
