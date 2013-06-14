# Manta

This is a pre-alpha project.


## Minimal Build instructions:

Manta must be built on a *nix flavor os. It is currently developed on Centos5 and
6 only.

An out-of-source build is recommended. Example:

> "
git clone $MANTA_GIT_URL manta
mkdir build
cd build
../manta/src/configure --prefix=/path/to/install
make
make install
"

All configure/make steps can be parallelized as follows:

> "
../manta/src/configure --prefix=/path/to/install --jobs=4
make -j4
make -j4 install
"

To see more options:
> "
../manta/src/configure --prefix=/path/to/install --jobs=4
make -j4
make -j4 install
"


### Known RHEL package dependencies:
* g++
* make
* zlib-devel


### Included build dependencies
Manta requires cmake 2.8+ and boost 1.49.0 Both of these packages are included in
the distribution and will be built from source if not present.


### Build configuration:
The build is configured for development rather than distribution. All unit tests are on by default and all builds include -Werror. If cppcheck is found any detected issue is converted to a build error.

