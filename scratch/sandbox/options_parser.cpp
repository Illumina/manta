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
