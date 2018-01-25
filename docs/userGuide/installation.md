Manta User Guide - Installation
===============================

[User Guide Home](README.md)

## Table of Contents

[//]: # (BEGIN automated TOC section, any edits will be overwritten on next source refresh)

* [Prerequisites to build from source](#prerequisites-to-build-from-source)
* [Runtime prerequisites](#runtime-prerequisites)
* [Operating System Guidelines](#operating-system-guidelines)
    * [Linux](#linux)
    * [OS X](#os-x)
    * [Windows](#windows)
* [Linux Package Additions](#linux-package-additions)
    * [Ubuntu 14.04 and 16.04](#ubuntu-1404-and-1604)
    * [Ubuntu 12.04](#ubuntu-1204)
    * [CentOS 7](#centos-7)
    * [CentOS 6](#centos-6)
* [Build procedure](#build-procedure)
    * [Workflow relocation](#workflow-relocation)
* [Demo](#demo)

[//]: # (END automated TOC section, any edits will be overwritten on next source refresh)


It is recommended to start from one of the [binary distributions on
the releases page][releases] if a suitable binary is available
(note that the CentOS 6 binary is expected to support a
large variety of linux systems). If building from source, then start from
the release distributions of the source code, also provided on the
[releases page][releases]. Cloning/archiving the source
directly from git could result in missing version number entries,
undesirably stringent build requirements or an unstable development
version between releases. Additional build notes for methods developers can
be found in the [developer guide][developerGuide].

[releases]:https://github.com/Illumina/manta/releases
[DeveloperGuide]:../developerGuide/README.md


### Prerequisites to build from source

[![Build Status][tcistatus]][tcihome]

[tcistatus]:https://travis-ci.org/Illumina/manta.svg?branch=master
[tcihome]:https://travis-ci.org/Illumina/manta

A C++11 compiler is required. These are the current minimum compiler versions enforced
by the build system:

* python 2.6+
* gcc 4.8+ OR clang 3.2+ (OR Visual Studio 2013+, see windows note below)
* libz (including headers)

### Runtime prerequisites

* python 2.6+

### Operating System Guidelines

##### Linux

Manta is known to build and run on the following linux distributions
(with additional packages as described below):

- Ubuntu 12.04, 14.04, 16.04
- CentOS 6, 7

##### OS X

Manta builds and passes basic tests on OS X 10.9, but full WGS analyses
are not tested for this platform.

##### Windows

Manta does not build or run on windows. Library-level compilation is
possible for Visual Studio users. See the [developer guide][DeveloperGuide] for details.

### Linux Package Additions

##### Ubuntu 14.04 and 16.04

    apt-get update -qq
    apt-get install -qq bzip2 gcc g++ make python zlib1g-dev

##### Ubuntu 12.04

    apt-get update -qq
    apt-get install -qq bzip2 gcc g++ make python python-software-properties zlib1g-dev
    # add gcc 4.8 from ubuntu ppa:
    add-apt-repository -y ppa:ubuntu-toolchain-r/test
    apt-get update -qq
    apt-get install -qq gcc-4.8 g++-4.8

    # Prior to build configuration, set CC/CXX to gcc 4.8:
    export CC=/usr/bin/gcc-4.8
    export CXX=/usr/bin/g++-4.8

##### CentOS 7

    yum install -y tar bzip2 make gcc gcc-c++ libstdc++-static zlib-devel

##### CentOS 6

    yum install -y tar wget bzip2 make gcc gcc-c++ zlib-devel

    # Add gcc 4.9 from developer tools v3:
    yum install -y centos-release-scl
    yum install -y devtoolset-3-gcc devtoolset-3-gcc-c++

    # Prior to build configuration, set CC/CXX to gcc 4.9:
    export CC=/opt/rh/devtoolset-3/root/usr/bin/gcc
    export CXX=/opt/rh/devtoolset-3/root/usr/bin/g++

### Build procedure

After acquiring a release distribution of the source code, the build
procedure is:

* Unpack source code
* Create and move to a separate `build` directory (out-of-source build is required.)
* Configure build
* Compile & Install

Example (building on 4 cores):

    wget https://github.com/Illumina/manta/releases/download/v${MANTA_VERSION}/manta-${MANTA_VERSION}.release_src.tar.bz2
    tar -xjf manta-${MANTA_VERSION}.release_src.tar.bz2
    mkdir build && cd build
    # Ensure that CC and CXX are updated to target compiler if needed, e.g.:
    #     export CC=/path/to/cc
    #     export CXX=/path/to/c++
    ../manta-${MANTA_VERSION}.release_src/configure --jobs=4 --prefix=/path/to/install
    make -j4 install

Note that there are two other dependencies: cmake and boost. These are different than the requirements discussed
above, in that they can optionally be provided by the user. They will automatically be built from source if not
detected. The minimum required versions of these tools for users planning to provide them to the build process are

* cmake 2.8.12+
* boost 1.58.0+ (must include static libraries)

...the build process will find an existing cmake version on the user's `PATH` and an existing boost installation
indicated by the environment variable `BOOST_ROOT`.

If not detected, then versions of cmake and boost will be built from source and installed to temporary locations under
the build directory automatically. This step can make installation more convenient, but does increase the time
requiered for configuration. To accelerate this process the configuration step can be parallelized over multiple cores,
as demonstrated in the example above with the`--jobs=4` argument to configure.

To see more configure options, run:

    ${MANTA_SRC_PATH}/configure --help

##### Workflow relocation

After Manta is built the installation directory can be relocated to
another directory.  All internal paths used in the workflow are
relative.

### Demo

To help verify a successful installation, Manta includes a small demo
data set and test script. After completing the installation steps
above, the demo can be run as follows:

    python ${MANTA_INSTALL_PATH}/bin/runMantaWorkflowDemo.py

This script creates a `MantaDemoAnalysis` directory under the current
working directory, runs Manta on a small demo dataset, and compares the
somatic structural variant output to an expected result.

See [the demo README](../../src/demo/README.md) for additional information
on the test script and data.

