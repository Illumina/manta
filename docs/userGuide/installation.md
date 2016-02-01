Manta User Guide - Installation
===============================

[User Guide Home](README.md)

It is recommended to start from one of the [binary distributions on
the Manta releases page] [releases] if a suitable version is available
(note that the CentOS 5 binary distribution is expected to support a
large variety of linux systems).  If building from source start from
the release distributions of the source code, also provided on the
[Manta releases page] [releases]. Cloning/archiving the source
directly from git could result in missing version number entries,
undesirably stringent build requirements or an unstable development
version between releases. Additional build notes for Manta developers can
be found in the [manta developer guide] [developerGuide].

[releases]:https://github.com/Illumina/manta/releases
[DeveloperGuide]:../developerGuide/README.md


### Prerequisites to build from source

[![Build Status] [tcistatus]] [tcihome]

[tcistatus]:https://travis-ci.org/Illumina/manta.svg?branch=master
[tcihome]:https://travis-ci.org/Illumina/manta

Manta requires a compiler supporting most of the C++11 standard. These
are the current minimum versions enforced by the build system:

* python 2.4+
* gcc 4.7+ OR clang 3.1+ (OR Visual Studio 2013+, see windows note below)

### Runtime prerequisites

* python 2.4+

### Operating System Guidelines

##### Linux 

Manta is known to build and run on the following linux distributions
(with additional packages as described below):

- Ubuntu 12.04,14.04
- CentOS 5,6,7

##### OS X

Manta builds and passes basic tests on OS X 10.9, but full WGS analyses
are not tested for this platform.

##### Windows

Manta does not build or run on windows. Library-level compilation is
possible for Visual Studio users. See the the [manta developer guide] [DeveloperGuide] for details.

### Linux Package Additions

##### Ubuntu 14.04

    apt-get update -qq
    apt-get install -qq gcc g++ make python

##### Ubuntu 12.04

    apt-get update -qq
    apt-get install -qq bzip2 gcc g++ make python python-software-properties
    # add gcc 4.8 from ubuntu ppa:
    add-apt-repository -y ppa:ubuntu-toolchain-r/test
    apt-get update -qq
    apt-get install -qq gcc-4.8 g++-4.8

    # Prior to build configuration, set CC/CXX to gcc 4.8:
    export CC=/usr/bin/gcc-4.8
    export CXX=/usr/bin/g++-4.8

##### CentOS 7

    yum install -y tar bzip2 make gcc gcc-c++

##### CentOS 5 and 6

    yum install -y tar wget bzip2 make gcc gcc-c++
    # add gcc 4.8 from developer tools v2:
    wget http://people.centos.org/tru/devtools-2/devtools-2.repo -O /etc/yum.repos.d/devtools-2.repo
    yum install -y devtoolset-2-gcc devtoolset-2-gcc-c++ devtoolset-2-binutils

    # Prior to build configuration, set CC/CXX to gcc 4.8:
    export CC=/opt/rh/devtoolset-2/root/usr/bin/gcc
    export CXX=/opt/rh/devtoolset-2/root/usr/bin/g++

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

Note that during the configuration step, the following dependencies
will be built from source if they are not found:

* cmake 2.8.0+
* boost 1.56.0+

To accelerate this process the configuration step can be parallelized
over multiple cores, as demonstrated in the example above with the
`--jobs=4` argument to configure.

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

