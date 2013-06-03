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

/// \author Chris Saunders
///

#pragma once

#include "boost/program_options.hh"

#include <string>
#include <vector>


namespace bpo = boost::program_options;


namespace manta {


/// base-class to help standardize option parsing
///
struct options_parser
{
    options_parser();

    virtual
    ~options_parser() {}

protected:

    void
    add_description(const bpo::options_description& desc) {
        _options.add(desc);
    }

    void
    parse_base(int argc, const char* argv[]);

    const std::vector<std::string>&
    get_unrecognized() {
        return _unrecognized_args;
    }

    const bpo::variables_map&
    get_vmap() const {
        return _vm;
    }

    void
    usage() const;

    virtual
    void
    usage_header(std::ostream& os) const {}

    virtual
    void
    usage_footer(std::ostream& os) const {}

private:
    bpo::options_description
    get_full_desc() const;

    bpo::options_description _options;
    bpo::options_description _help_options;

    //output:
    bpo::variables_map _vm;
    std::vector<std::string> _unrecognized_args;
};


}
