Manta Structural Variant Caller
===============================

Manta calls structural variants (SVs), small indels and large insertions from paired-end
sequencing reads. The method is designed for very rapid/low-cost analysis of individuals
or tumor-normal sample pairs. Manta combines paired and split-read evidence during SV
discovery and scoring to improve accuracy, but does not require split-reads
or successful breakpoint assemblies to report a variant in cases where there is
strong evidence of an imprecise variant. It provides scoring models
for germline variants in individual diploid samples and somatic variants in matched
tumor-normal sample pairs. Manta can detect all classes
of structural variants and indels which can be identified in the absence of copy number
analysis and large-scale assembly.
See the [user guide] [UserGuide] for a full description of
capabilities and limitations.

[UserGuide]:src/markdown/mantaUserGuide.md

Build instructions
------------------

[![Build Status](https://travis-ci.org/sequencing/manta.svg?branch=master)](https://travis-ci.org/sequencing/manta)

It is recommended to start from one of the [binary distributions
on the Manta releases page] [releases] if a suitable version is available (note
that the CentOS 5 binary distribution is expected to support a large variety
of linux systems). 
If building from source start from the release distributions of the source code,
also provided on the [Manta releases page] [releases]. Cloning/archiving the
source directly from git could result in missing version number entries, undesirably
stringent build requirements or an unstable development version between
releases. Additional build notes for contributors can be found below.

Note that this README is _NOT_ part of a tagged source-code release.

[releases]:https://github.com/sequencing/manta/releases

### Build prerequisites:

Manta requires a compiler supporting most of the C++11 standard. These are the
current minimum versions enforced by the build system:

* python 2.4+
* gcc 4.7+ OR clang 3.2+ (OR Visual Studio 2013+, see windev note below)
* libz (including headers)

### Runtime prerequisites

* python 2.4+

### Operating System Guidelines

##### Linux 

Manta is known to build and run on the following linux distributions
(with additional packages as described below):

- Ubuntu 12.04,14.04
- CentOS 5,6,7

##### OS X

Build and run testing on OS X is ad-hoc. Last successful build was on
v0.26.1-15-gfe9a48d using OS X 10.9.5 / Xcode 6.2

##### Windows

Manta does not build or run on windows. Library-level compilation is
possible for Visual Studio users. See Developer section below.

### Linux Package Additions

##### Ubuntu 14.04

    sudo apt-get update -qq
    sudo apt-get install -qq gcc g++ make zlib1g-dev python

##### Ubuntu 12.04

    apt-get update -qq
    apt-get install -qq bzip2 gcc g++ make zlib1g-dev python python-software-properties
    # add gcc 4.8 from ubuntu ppa:
    add-apt-repository -y ppa:ubuntu-toolchain-r/test
    apt-get update -qq
    apt-get install -qq gcc-4.8 g++-4.8

    # Prior to manta configuration, set CC/CXX to gcc 4.8:
    export CC=/usr/bin/gcc-4.8
    export CXX=/usr/bin/g++-4.8

##### CentOS 7

    yum install -y tar bzip2 make gcc gcc-c++ zlib-devel

##### CentOS 5 and 6

    yum install -y tar wget bzip2 make gcc gcc-c++ zlib-devel
    # add gcc 4.8 from developer tools v2:
    wget http://people.centos.org/tru/devtools-2/devtools-2.repo -O /etc/yum.repos.d/devtools-2.repo
    yum install -y devtoolset-2-gcc devtoolset-2-gcc-c++ devtoolset-2-binutils

    # Prior to manta configuration, set CC/CXX to gcc 4.8:
    export CC=/opt/rh/devtoolset-2/root/usr/bin/gcc
    export CXX=/opt/rh/devtoolset-2/root/usr/bin/g++

### Build procedure

After acquiring a release distribution of the source code, the build procedure is:

* Unpack source code
* Create and move to a separate build directory (out-of-source build is required.)
* Configure build
* Compile & Install

Example (building on 4 cores):

    wget https://github.com/sequencing/manta/releases/download/vA.B.C/manta-A.B.C.release_src.tar.bz2
    tar -xjf manta-A.B.C.release_src.tar.bz2
    mkdir build && cd build
    # Ensure that CC and CXX are updated to target compiler if needed 
    ../manta-A.B.C/src/configure --jobs=4 --prefix=/path/to/install
    make -j4 install

Note that during the configuration step, the following dependencies will be
built from source if they are not found:

* cmake 2.8.0+
* boost 1.53.0

To accelerate this process the configuration step can be parallelized over
multiple cores, as demonstrated in the example above with the `--jobs=4`
argument to configure.

To see more configure options, run:

    ${MANTA_SRC_PATH}/configure --help

##### Workflow relocation

After Manta is built the installation directory can be relocated to another directory.
All internal paths used in the workflow are relative.

### Demo / Build Verification 

To help verify a successful installation, Manta includes a small demo data
set and test script. After completing the installation steps above, run the
demo as follows:

    ${MANTA_INSTALL_PATH}/bin/runMantaWorkflowDemo.bash

The demo data includes aligned reads near the breakends of an HCC1954 somatic translocation
([COSMIC] id: COST16011). The demo script runs Manta on these regions and verifies that the
somatic translocation is called as expected. Data, expected results and additional details 
can be found in the Manta installation demo directory:

    ${MANTA_INSTALL_PATH}/share/demo/manta

[COSMIC]:http://cancer.sanger.ac.uk/cosmic

Data analysis and Interpretation
--------------------------------

After completing installation, see the [Manta user guide] [UserGuide] for instructions on
how to run Manta, interpret results, and a high-level overview of the method.

The user guide is also available within any Manta binary distribution here:

    ${MANTA_INSTALL_PATH}/doc/html/mantaUserGuide.html


Contributor build configuration
-------------------------------

When the Manta source is cloned from github, it is configured for development
rather than user distribution. As such, builds are strict: all warnings are
treated as errors and if cppcheck is found any detected issue is converted
to a build error.

#### Source documentation

If doxygen is found in the path (and optionally dot as well) during build
configuration, then c++ documentation is available as an additional "doc"
target for the makefile:

    make doc

There is no installation for the documentation outside of the build directory,
the root doxygen page after completing this target will be:

    ${MANTA_BUILD_DIR}/c++/doxygen/html/index.html

#### Improving build time

##### ccache

Manta's build system is configured to use ccache whenever this is found in the path

##### Bundled dependencies

Note that during the configuration step, the following compilation
dependencies will be built from source if they are not found:

* cmake 2.8.0+
* boost 1.53.0

To avoid the extra time associated with this step, ensure that (1) cmake 2.8.0+
is in your PATH and (2) BOOST\_ROOT is defined to point to boost 1.53.0 (the
boost version is required to be an exact match).

#### Windows development support

Manta does not link or run on windows. However, the build system does
facilitate Visual Studio (VS) users. When cmake configuration is run on
windows, all linking is disabled and third party libraries are unpacked
for header include access, but are not compiled. Cmake VS solutions allow
the c++ code to be browsed, analyzed and compiled to the library level.
Note that unit test codes are compiled to libraries but cannot be run.

C++11 features used by manta require at least VS2013. Windows installations
of cmake and zlib are also required to configure and compile. Windows zlib
is provided by the 
[gnuwin32 package](http://gnuwin32.sourceforge.net/packages/zlib.htm).
among others.

