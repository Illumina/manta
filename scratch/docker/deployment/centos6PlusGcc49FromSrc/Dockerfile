#
# This is a simple image used to assist with deploying portable
# binaries, by allowing us to build in a virtual centos 6
# environment.
#
# At present we simply add all of the project's build requirements onto
# the centos6 base image.
#
# This version is extended to build gcc 4.9 from source. The purpose is
# to investigate a small performance gap in the deployed binary compared
# to our non-docker based deployment:
#

FROM centos:6

# add standard centos6 packages
# adding cmake here for extra speed even though starka will bootstrap this
RUN yum install -y bzip2 gcc gcc-c++ tar wget zlib-devel cmake git

# build gcc 4.9 from source
RUN mkdir -p /download/gcc-4.9.3 && cd /download/gcc-4.9.3 && wget ftp://ftp.gnu.org/gnu/gcc/gcc-4.9.3/gcc-4.9.3.tar.bz2 && \
    tar -xjf gcc-4.9.3.tar.bz2 && cd gcc-4.9.3 && ./contrib/download_prerequisites && \
    cd .. && mkdir -p build && cd build && \
    ../gcc-4.9.3/configure \
        --prefix=/opt/gcc-4.9.3 \
        --disable-multilib \
        --disable-bootstrap \
        --enable-cloog-backend=isl \
        --enable-lto \
        --with-system-zlib \
        --enable-languages=c,c++ && \
    make -j2 && make install && cd / && rm -rf /download && \
    GCC_PATH=/opt/gcc-4.9.3 && SPECS_PATH=${GCC_PATH}/lib/gcc/x86_64-unknown-linux-gnu/4.9.3 && SPECS_FILE=${SPECS_PATH}/specs && \
    ${GCC_PATH}/bin/gcc -dumpspecs > ${SPECS_FILE} && echo '*link:'$'\n'+ -R ${GCC_PATH}/lib64$'\n'  >> ${SPECS_FILE}

# Prior to build configuration, set CC/CXX to new gcc:
ENV CC /opt/gcc-4.9.3/bin/gcc
ENV CXX /opt/gcc-4.9.3/bin/g++
