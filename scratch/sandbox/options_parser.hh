//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2018 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
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
