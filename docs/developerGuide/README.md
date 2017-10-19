Manta Developer Guide
=====================

## Table of Contents

[//]: # (BEGIN automated TOC section, any edits will be overwritten on next source refresh)

* [Scope](#scope)
* [Developer Build Notes](#developer-build-notes)
  * [Building from source repository vs. versioned code distribution](#building-from-source-repository-vs-versioned-code-distribution)
  * [Static analysis](#static-analysis)
  * [Source auto-documentation](#source-auto-documentation)
  * [Improving build time](#improving-build-time)
    * [ccache](#ccache)
    * [Bundled dependencies](#bundled-dependencies)
* [General Debugging Notes](#general-debugging-notes)
  * [General Debugging: Address Sanitizer](#general-debugging-address-sanitizer)
  * [General Debugging: Inspecting temporary files](#general-debugging-inspecting-temporary-files)
  * [Windows development support](#windows-development-support)
  * [Automating Portable Binary Builds for Linux](#automating-portable-binary-builds-for-linux)
* [Coding Guidelines](#coding-guidelines)
  * [Source formatting](#source-formatting)
  * [Git conventions](#git-conventions)
    * [Commit messages](#commit-messages)
    * [Commit consolidation](#commit-consolidation)
  * [Changelog conventions](#changelog-conventions)
* [Branching and release tagging guidelines](#branching-and-release-tagging-guidelines)
  * [Error handling](#error-handling)
    * [General Policies](#general-policies)
    * [Exception Details](#exception-details)
    * [Logging](#logging)
  * [Unit tests](#unit-tests)
* [IDE support](#ide-support)
  * [Clion](#clion)
* [Special Topic Guides](#special-topic-guides)

[//]: # (END automated TOC section, any edits will be overwritten on next source refresh)


## Scope

This guide provides:
* protocols for contributing new or modified methods
* methods to debug stability or runtime issues
* methods to debug suspected false or missing variant calls
* high-level architectural documentation

Information is added as pertinent questions/discussions come up in the contributor community,
so this guide is not intended to provide complete coverage of the above topics.

For end user documentation describing how to run an analysis and interpret its output,
please see the [User Guide](../userGuide/README.md).

## Developer Build Notes

The following section provides a supplement to the standard build
instructions including additional details of interest to methods
developers.

### Building from source repository vs. versioned code distribution

When the source repository is cloned from github, it is configured for development
rather than user distribution. In this configuration all builds are strict
such that:
* all warnings are treated as errors
* if cppcheck is found any detected cppcheck issue is converted to a build error (see details below)

Note that all unit tests are always run and required to pass for the build
procedure to complete.

### Static analysis

When the build is configured for development, static analysis is run on all project c++ source using cppcheck, so long
as an appropriate cppcheck version can be found in the user's path. This static analysis step is configured as follows:
* cppcheck will be used if it can be found in the user's path and the version found is at least 1.69
* The script which runs and interprets cppcheck is [run_cppcheck.py](../../src/srcqc/run_cppcheck.py).
* Any cppcheck warning will be treated as an error. A few warnings are suppressed, depending on the cppcheck version.
* All cppcheck warnings are reformatted to follow standard gcc/clang error message style to work correctly in most IDEs.
* All project c++ code will be analyzed, but third-party/redistributed packages are ignored.

### Source auto-documentation

If doxygen is found in the path (and optionally dot as well) during
build configuration, then c++ documentation is available as an
additional "doc" target for the makefile:

    make doc

There is no installation for the documentation outside of the build
directory, the root doxygen page after completing this target will be:

    ${MANTA_BUILD_PATH}/src/c++/doxygen/html/index.html

### Improving build time

#### ccache

The build system is configured to use ccache whenever this is
found in the path

#### Bundled dependencies

Note that during the configuration step, the following dependencies will be
built from source if they are not found:

* cmake 2.8.12+
* boost 1.58.0+

To avoid the extra time associated with this step, ensure that (1)
cmake 2.8.12+ is in your PATH and (2) BOOST\_ROOT is defined to point
to boost 1.58.0 or newer.

## General Debugging Notes

### General Debugging: Address Sanitizer

The build system offers first-class support for google address sanitizer
when a supporting compiler is detected. To use this mode, start a fresh
installation process with the additional configure option `--build-type=ASan`,
for example:

    ../configure --jobs=4 --prefix=/path/to/install --build-type=ASan

### General Debugging: Inspecting temporary files

Manta's configuration step includes an extended option to keep all temporary
files which would normally be deleted by the workflow as it runs. Keeping these
files may be helpful in various debugging scenarios. To turn on this option, add
`--retainTempFiles` as a configuration argument:

    configManta.py [other_options...] --retainTempFiles

### Windows development support

Manta does not link or run on windows. However, the build system does
facilitate Visual Studio (VS) users. When cmake configuration is run
on windows, all linking is disabled and most third party libraries are
unpacked for header include access, but are not compiled. Cmake VS
solutions allow the c++ code to be browsed, analyzed and compiled to
the library level.  Note that unit test codes are compiled to
libraries but cannot be run.

C++11 features in use require at least VS2013. A Windows
installation of cmake is also required to configure and compile.
Note that the minimum cmake version is 3.1.0 for Windows.

### Automating Portable Binary Builds for Linux

A script is provided to enable a dockerized build process which issues
binary tarballs for a variety of OS/compiler configurations, with the
goal of creating a reasonably portable linux binary build by using a
relatively old OS/glibc version. To use this script, ensure you
have permission to `docker run` on the current system and execute the
following script:

```
${SOURCE_PATH}/scratch/docker/deployment/dockerBuildBinaryTarball.bash ${SOURCE_PATH2} ${BINARY_BUILD_PREFIX}
```

The term `${SOURCE_PATH2}` can point to the current git repository (ie. `${SOURCE_PATH}`),
or to an extracted source release tarball previously created using the script:

```
${SOURCE_PATH}/scratch/make_release_tarball.bash
```

The choice of virtualized build environment is hard-coded in the deploy script for the time being,
see the `builderImage` variable.

## Coding Guidelines

### Source formatting

* Basic formatting restrictions on c++ code:
  * spaces instead of tabs
  * 4-space indents
  * "ANSI" bracket style
* Note the above restrictions are enforced by an astyle script which is occasionally run on the master branch (see [run_cxx_formatter.bash](../../scratch/source_check_and_format/run_cxx_formatter.bash))
* Name formatting for all newly introduced code:
  * Make variable/type names self-documenting whereever this is practical, e.g. sampleCount, breakpointRegion, etc.
  * Lowercase camelCase variable names
  * Uppercase CamelCase type names
* Otherwise, follow local code conventions

### Git conventions

#### Commit messages

All git commit messages should be prepended with either the associated JIRA or github issue id. For example:

```
MANTA-123 Improve insertion genotype accuracy

Improve assembly and realignemnt of large insertions to reduce hom->het undercall.
```

Very minor updates may be made without an associated ticket. In these cases providing a category prefix for
the minor change is still encouraged to provide context, for instance:

```
docs: Fix paths in user guide examples
```

```
build: Close new warnings from clang-4.0.0
```

All git commit messages should conform to practices outlined here:

http://chris.beams.io/posts/git-commit/


#### Commit consolidation

On any single-developer research branch, history editing is encouraged within the branch to collapse bugs and
build a more clear feature-by-feature story for other other developers to follow.

In all other situations history editing is discouraged and a conventional merge-based workflow is preferred.

### Changelog conventions

The primary function of the changelog is to help end-users weigh the benefits and risks or updating to a newer version.

To this end:
- Changelog entries should be made for any major branch merge, bug fix, or generally something that would change a user's interaction at the configuration step, change the format/interpretation of the output or change the variant calling performance.
- The project changelog follows many of the conventions suggested in keepachangelog. Consistent with this formatting, changelog entries should be more descriptive and end-user focused than git commit entries. JIRA or github ticket IDs should be listed at the end of the Changelog description to help developers link issues, without making this the emphasis of each changelog entry.
- Each major branch with an end-user impact should add a changelog entry to the Unreleased section of the changelog prior to merging the branch.

## Branching and release tagging guidelines

All features and bugfixes are developed on separate branches. Branch names should contain the corresponding JIRA ticket
id or contain the key "github${issueNumber}' to refer to the corresponding issue on github.com. After code
review and testing, all branches intended for release are merged to the 'develop' branch. Releases are tagged
on the development branch, and the corresponding release tag can be pushed to the 'master' branch. Thus the master
branch always represents the most recent stable release. As required version specific development branches are
created to support continued bugfix iterations on older releases. These branches replace the patch number in the
release version with an "x", for instance the "v2.4.x" branch would be used to provide bugfix updates "v2.4.2",
"v2.4.3", etc.

### Error handling

#### General Policies

* Exceptions with informative contextual details are encouraged whenever possible.
* To quickly express invariants it is acceptable to add `assert()`'s first, and transition to exceptions as code stabilizes.
* Note that the build process will never define `NDEBUG` to compile out assert statements, even in release code.
* Exceptions are never thrown with the intent to recover -- this is not a web browser. The goal is to:
  * Fail at the first sign of trouble.
  * Provide as much helpful contextual information as possible, including context from multiple layers of the stack.
* Warnings are discouraged. If considering a warning you should probably just fail per the above policy.

#### Exception Details

* Preferred exception pattern is to use an internal class derived from `boost::exception`:

```c++

#include "common/Exceptions.hh"

#include <sstream>

void
foo(const char* name)
{
    using namespace illumina::common;

    std::ostringstream oss;
    oss << "ERROR: unrecognized variant scoring model name: '" << name << "'\n";
    BOOST_THROW_EXCEPTION(LogicException(oss.str()));
}
```

* Context at the original throw site is often supplemented by a 'catch and release' block to add
information at a few critical points on the stack. Typically this is information which
is unavailable at the throw site. Example code is:

```c++
try
{
    realign_and_score_read(_opt,_dopt,sif.sample_opt,_ref,realign_buffer_range,rseg,sif.getIndelBuffer());
}
catch (...)
{
    log_os << "ERROR: Exception caught in align_pos() while realigning segment: "
	   << static_cast<int>(r.second) << " of read: " << (*r.first) << "\n";
    throw;
}
```

#### Logging

* At the workflow (python) layer, please write all logging messages through pyflow's logging interface as follows:
```python
self.flowLog("Initiating Starling workflow version: %s" % (__version__)
```

* At the binary (c++) layer, there is no logger at present. Direct all error messaging to `std::cerr`.

### Unit tests

* Unit tests are enabled for a subset of the c++ code
* All tests use the boost unit test framework
* All unit tests are required to run and pass as part of every build (including end-user builds)
* Unit tests are already enabled for every library "test" subdirectory, additional tests in these directories will be automatically detected
  * Example [svgraph unit tests directory](../../src/c++/lib/svgraph/test)

## IDE support

Little support for any specific IDE is provided, except as made available by cmake generators. IDE-specific configuration files maintained in the project are described below.

### Clion

A subset of code formatting settings which can be imported into Clion are available in the configuration file

`${STRELKA_REPO_PATH}/scratch/ideConfiguration/CLion/referenceCodeStyleSettings.xml`

..note that the automated `astyle` formatting settings still define the project defaults, the above configuration simply provides a starting point for CLion which is closer to the project's formatting norms.


## Special Topic Guides

The following items provide more in-depth details on a subsection of the methods/debugging protocol, etc.

* [Alignment Library](alignment.md)
* [Breakend Graph Queries](breakendGraph.md)
* [Debug Single SV](debugSingleSV.md)
* [Debug Full Manta Run](debugSingleSV.md)
* [Manta VCF ID field](ID.md)
* [Test assembler directly from BAM input](testAssembler.md)

