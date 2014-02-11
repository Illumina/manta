// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

/// \author Chris Saunders
///

#include "blt_util/blt_exception.hh"
#include "blt_util/log.hh"
#include "blt_util/sig_handler.hh"
#include "common/Exceptions.hh"
#include "manta/Program.hh"

#include <cstdlib>

#include <iostream>


static
void
dump_cl(int argc,
        char* argv[],
        std::ostream& os)
{

    os << "cmdline:";
    for (int i(0); i<argc; ++i)
    {
        os << ' ' << argv[i];
    }
    os << std::endl;
}



static
void
post_catch(int argc,
           char* argv[],
           std::ostream& os)
{

    os << "...caught in program.run()\n";
    dump_cl(argc,argv,log_os);
    exit(EXIT_FAILURE);
}


namespace manta
{



int
Program::
run(int argc, char* argv[]) const
{

    try
    {
        std::ios_base::sync_with_stdio(false);

        std::string cmdline;
        for (int i(0); i<argc; ++i)
        {
            if (i) cmdline += ' ';
            cmdline += argv[i];
        }

        initialize_blt_signals(name(),cmdline.c_str());

        runInternal(argc,argv);
    }
    catch (const blt_exception& e)
    {
        log_os << "FATAL_ERROR: " << name() << " EXCEPTION: " << e.what() << "\n";
        post_catch(argc,argv,log_os);
    }
    catch (const illumina::common::ExceptionData& e)
    {
        log_os << "FATAL_ERROR: " << name() << " EXCEPTION: "
               << e.getContext() << ": " << e.getMessage() << "\n";
        post_catch(argc,argv,log_os);
    }
    catch (const boost::exception& e)
    {
        log_os << "FATAL_ERROR: " << name() << " EXCEPTION: "
               << boost::diagnostic_information(e) << "\n";
        post_catch(argc,argv,log_os);
    }
    catch (const std::exception& e)
    {
        log_os << "FATAL_ERROR: EXCEPTION: " << e.what() << "\n";
        post_catch(argc,argv,log_os);
    }
    catch (...)
    {
        log_os << "FATAL_ERROR: UNKNOWN EXCEPTION\n";
        post_catch(argc,argv,log_os);
    }
    return EXIT_SUCCESS;
}

}
