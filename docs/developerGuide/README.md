Manta Developer Guide
=====================

This guide provides information for manta development, including protocols for
contirbuting new methods, debugging stability problems, suspected false or missing variant calls and some high-level internal methods documentation.

For end user documentation describing how to run Manta and interpret its output, please see the [Manta User Guide](../userGuide/README.md).

## Table of Contents
* [Developer Build Notes](#developer-build-notes)
* [Coding Guidelines](#coding-guidelines)
* [Special Topic Guides](#special-topic-guides)

## Developer Build Notes

The following section provides a supplement to the standard build
instructions including additional details of interest to methods
developers.

### Building from source repository vs. versioned code distribution:

When Manta is cloned from github, it is configured for development
rather than user distribution. In this configuration all builds are strict
such that:
* all warnings are treated as errors
* if cppcheck is found any detected cppcheck issue is converted to a build error

Note that in all build configurations, all of Manta's unit tests are run and required
to pass as part of the default build procedure.

### Source auto-documentation

If doxygen is found in the path (and optionally dot as well) during
build configuration, then c++ documentation is available as an
additional "doc" target for the makefile:

    make doc

There is no installation for the documentation outside of the build
directory, the root doxygen page after completing this target will be:

    ${MANTA_BUILD_PATH}/c++/doxygen/html/index.html

### Improving build time

#### ccache

The build system is configured to use ccache whenever this is
found in the path

#### Bundled dependencies

Note that during the configuration step, the following dependencies will be
built from source if they are not found:

* cmake 2.8.0+
* boost 1.56.0+

To avoid the extra time associated with this step, ensure that (1)
cmake 2.8.0+ is in your PATH and (2) BOOST\_ROOT is defined to point
to boost 1.56.0 or newer.

### General Debugging: Address Sanitizer

The build system offers first-class support for google address sanitizer
when a supporting compiler is detected. To use this mode, start a fresh
installation process with the additional configure option `--build-type=ASan`,
extending from the configuration example in the above build instructions, use:

    ../manta-A.B.C.release_src/src/configure --jobs=4 --prefix=/path/to/install --build-type=ASan

### Windows development support

Manta does not link or run on windows. However, the build system does
facilitate Visual Studio (VS) users. When cmake configuration is run
on windows, all linking is disabled and most third party libraries are
unpacked for header include access, but are not compiled. Cmake VS
solutions allow the c++ code to be browsed, analyzed and compiled to
the library level.  Note that unit test codes are compiled to
libraries but cannot be run.

C++11 features used by manta require at least VS2013. A Windows
installation of cmake is also required to configure and compile.
Note that the minimum cmake version for Windows is 3.1.0

## Coding Guidelines

### Source formatting

* Basic formatting restrictions on c++ code:
  * spaces instead of tabs
  * 4-space indents
  * "ANSI" bracket style
* Note the above restrictions are enforced by an astyle script which is occasionally run on the master branch (see [run_cxx_formatter.bash](../../scratch/source_check_and_format/run_cxx_formatter.bash))
* Otherwise, just follow local code conventions

### Unit tests
* Unit tests are enabled for a subset of the c++ code in manta
* All tests use the boost unit test framework
* All unit tests are required to run and pass as part of every build (including end-user builds)
* Unit tests are already enabled for every library "test" subdirectory, additional tests in these directories will be automatically detected 
  * Example [svgraph unit tests directory](../../src/c++/lib/svgraph/test)

## Special Topic Guides

The following items provide more in-depth details on a subsection of the methods/debugging protocol, etc.

* [Alignment Library](alignment.md)
* [Breakend Graph Queries](breakendGraph.md)
* [Debug Single SV](debugSingleSV.md)
* [Debug Full Manta Run](debugSingleSV.md)
* [Manta VCF ID field](ID.md)

