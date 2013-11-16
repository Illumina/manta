Manta
=====

Version: NOT RELEASED

Manta is a software package to call structural variants and indels from short
paired-end sequencing reads. It combines paired-end and split read evidence to
provide genotype and quality score for variants in single diploid samples. It
will also call somatic variants when a matched tumor sample is provided. Manta
can detect all classes of structural variants which are identifiable in the
absence of copy number analysis. See the user guide for a full description of
capabilities and limitations.

_Note: This software under development and not supported for general use._


End-user build instructions
---------------------------

Manta has been built and tested on linux systems only. It is currently
maintained for Centos5,6 and Ubuntu 12.04.

For Manta end users it is strongly recommended to start from one of the
release distributions of the source code provided on the Manta [releases]
page. Acquiring the source via a git clone or archive could result in missing
version number entries, undesirably stringent build requirements, or an unstable
development intermediate between releases.

Note that this is _NOT_ an end-user release distribution.

[releases]:https://github.com/StructuralVariants/manta/releases

### Build prerequisites

* python 2.4+
* gcc 4.1+ OR clang 3.2+
* libz (including headers)

### Runtime prerequisites

* python 2.4+

### Prerequisite package names (RHEL/Centos)

* g++
* make
* zlib-devel

### Build procedure

After acquiring a release distribution of the source code, the build procedure is:

* Unpack the source code
* Create a separate `build` directory (note an out-of-source build is
  required.)
* Configure the build
* Compile
* Install

Example:

    wget https://github.com/StructuralVariants/manta/releases/download/vA.B.C/manta-A.B.C.tar.bz2
    tar -xjf manta-A.B.C.tar.bz2
    mkdir build
    cd build
    ../manta-A.B.C/src/configure --prefix=/path/to/install
    make
    make install

Note that during the configuration step, Manta will build the following
compilation dependencies if these are missing:

* cmake 2.8.0+
* boost 1.49.0

To avoid this extra step, ensure that (1) cmake is in your PATH and (2)
BOOST\_ROOT is defined to point to boost 1.49.0 (the boost version is required to
be an exact match). If either of these dependencies are not found, they will be
built during the configuration step, To accelerate this process it may be
desirable to parallelize the configure step over multiple cores. To do so
provide the `--jobs` argument to the configuration script. For example:

    ${MANTA_SRC_PATH}/configure --prefix=/path/to/install --jobs=4

Compiling Manta itself can also be parallelized by the standard make procedure, e.g.

    make -j4
    make install

To see more configure options, run:

    ${MANTA_SRC_PATH}/configure --help


Data analysis and Interpretation
--------------------------------

After completing the installation, please see the Manta User Guide for
instructions on how to use Manta, interpret results, and a high-level overview
of the method. The user guide can be found in:

    ${MANTA_INSTALL_PATH}/doc/html/mantaUserGuide.html


Developer build configuration
-----------------------------

When the Manta source is cloned from github, it is configured for development
rather than end-user distribution. As such, all builds include -Werror. If
cppcheck is found any detected issue is converted to a build error.

