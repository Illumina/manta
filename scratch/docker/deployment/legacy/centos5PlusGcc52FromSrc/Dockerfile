#
# This is a simple image used to assist with deploying portable
# binaries, by allowing us to build in a virtual centos 5
# environment.
#
# At present we simply add all of the projects build requirements onto
# the centos5 base image.
#
# This version is extended to build gcc 5.2 from source. The purpose is
# to investigate a small performance gap in the deployed binary compared
# to our non-docker based deployment:
#

FROM astj/centos5-vault

# add standard centos5 packages && swap in newer cmake as default:
RUN yum install -y bzip2 make gcc gcc-c++ tar wget zlib-devel git && \
    yum install -y epel-release && \
    yum install -y cmake28 && cd /usr/bin && ln -s cmake28 cmake

## update binutils
## 
## this was a test to get lto working in centos5
#
#RUN mkdir -p /download/binutils && cd /download/binutils && wget ftp://ftp.gnu.org/gnu/binutils/binutils-2.25.tar.bz2 && \
#    tar -xjf binutils-2.25.tar.bz2 && cd binutils-2.25 && ./configure --prefix=/opt/binutils-2.25 --enable-lto && \
#    make -j2 && make install && cd / && rm -rf /download
#
## other lto stuff for gcc configure below:
#        --with-ld=/opt/binutils-2.25/bin/ld \
#        --with-as=/opt/binutils-2.25/bin/as \

# build gcc 5.2 from source
RUN GCCVER=5.2.0 && mkdir -p /download/gcc-${GCCVER} && cd /download/gcc-${GCCVER} && wget ftp://ftp.gnu.org/gnu/gcc/gcc-${GCCVER}/gcc-${GCCVER}.tar.bz2 && \
    tar -xjf gcc-${GCCVER}.tar.bz2 && cd gcc-${GCCVER} && ./contrib/download_prerequisites && \
    cd .. && mkdir -p build && cd build && \
    ../gcc-${GCCVER}/configure \
        --prefix=/opt/gcc-${GCCVER} \
        --disable-multilib \
        --disable-bootstrap \
        --enable-lto \
        --with-system-zlib \
        --enable-languages=c,c++ && \
    make -j2 && make install && cd / && rm -rf /download && \
    GCC_PATH=/opt/gcc-${GCCVER} && SPECS_PATH=${GCC_PATH}/lib/gcc/x86_64-unknown-linux-gnu/${GCCVER} && SPECS_FILE=${SPECS_PATH}/specs && \
    ${GCC_PATH}/bin/gcc -dumpspecs > ${SPECS_FILE} && echo '*link:'$'\n'+ -R ${GCC_PATH}/lib64$'\n'  >> ${SPECS_FILE}

# Prior to build configuration, set CC/CXX to new gcc version: 
ENV CC /opt/gcc-5.2.0/bin/gcc
ENV CXX /opt/gcc-5.2.0/bin/g++
