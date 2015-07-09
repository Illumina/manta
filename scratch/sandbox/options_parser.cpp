// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2015 Illumina, Inc.
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

#include "manta_common/program_options.hh"


namespace manta {


options_parser::
options_parser()
     : _options("Allowed Options"),
       _help_options("Help")

{
    _help_options.add_options()("help,h", "produce help message and exit");
    _help_options.add_options()("version,v", "print program version information");
}



bpo::options_description
options_parser::
get_full_desc() const {

    bpo::options_description full(_options);
    full.add(_help_options);
    return full;
}



void
options_parser::
usage() const {
    usage_header(log_os);
    log_os << get_full_desc();
    usage_footer(log_os);
}


 void
 options_parser::
 parse_base(int argc, const char* argv[]) {

     try
     {
         bpo::parsed_options parsed(bpo::command_line_parser(argc,argv).options(get_full_desc).allow_unregistered().run());
         bpo::store(parsed,_vm);
         bpo::notify(_vm);

         // allow remaining options to be parsed using starling command-line parser:
         _unrecognized_args = bpo::collect_unrecognized(parsed.options,bpo::include_positional);


         if (_vm.count("help"))
         {
             usage();
             exit(2);
         }
     }
     catch (const boost::program_options::multiple_values &e)
     {
         log_os << usage() << '\n';
         log_os << "Failed to parse the options: " << e.what() << ": " << e.get_option_name() << std::endl;
         exit(2);
     }
     catch (const boost::program_options::multiple_occurrences &e)
     {
         log_os << usage() << '\n';
         log_os << "Failed to parse the options: " << e.what() << ": " << e.get_option_name() << std::endl;
         exit(2);
     }
     catch (const boost::program_options::required_option &e)
     {
         log_os << usage() << '\n';
         log_os << "Failed to parse the options: " << e.what() << ": " << e.get_option_name() << std::endl;
         exit(2);
     }
     catch (const std::exception &e)
     {
         log_os << usage() << '\n';
         log_os << "Failed to parse the options: " << e.what() << std::endl;
         exit(2);
     }
 }


}
