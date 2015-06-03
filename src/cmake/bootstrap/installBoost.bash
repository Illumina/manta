#!/usr/bin/env bash
#
# Manta
# Copyright (c) 2013-2014 Illumina, Inc.
#
# This software is provided under the terms and conditions of the
# Illumina Open Source Software License 1.
#
# You should have received a copy of the Illumina Open Source
# Software License 1 along with this program. If not, see
# <https://github.com/sequencing/licenses/>
#

################################################################################
##
## Script to install boost
##
## author Come Raczy
##
################################################################################

set -o nounset
set -o pipefail


REDIST_DIR=$1
INSTALL_DIR=$2
if [[ $# -ge 3 ]] ; then PARALLEL=$3 ; else PARALLEL=1 ; fi

# test that these values are defined:
test=${THIS_BOOST_VERSION}
test=${THIS_BOOST_BUILD_COMPONENTS}

script_dir="$(dirname "$0")"
source $script_dir/common.bash

BUILD_DIR=${INSTALL_DIR}/build
BIN_DIR=${INSTALL_DIR}/bin
LIB_DIR=${INSTALL_DIR}/lib
INCLUDE_DIR=${INSTALL_DIR}/include

SCRIPT=`basename "$0"`
VERSION=`echo ${THIS_BOOST_VERSION} | sed "s/\./_/g"`
SOURCE_PREFIX=boost_${VERSION}
SOURCE_TARBALL=${REDIST_DIR}/${SOURCE_PREFIX}.tar.bz2
TARBALL_COMPRESSION=j
SOURCE_DIR=${BUILD_DIR}/${SOURCE_PREFIX}

BOOST_LIBRARY_LIST=$(echo ${THIS_BOOST_BUILD_COMPONENTS} | sed "s/;/,/g")

if [ -z "${BOOTSTRAP_OPTIONS+xxx}" ]; then BOOTSTRAP_OPTIONS=""; fi
if [ -z "${BJAM_OPTIONS+xxx}" ]; then BJAM_OPTIONS=""; fi


common_options $@

if [[ $CLEAN ]] ; then
    ilog removing $SOURCE_DIR
    rm -rf $SOURCE_DIR ${INCLUDE_DIR}/boost ${LIB_DIR}/libboost_*.{a,so}
    exit 0
fi

common_create_source

if [ $CMAKE_CXX_COMPILER_ID == "GNU" ]; then
    echo "using gcc : : \"$CXX\" ;" >> ${SOURCE_DIR}/tools/build/v2/user-config.jam
fi

cd ${SOURCE_DIR} \
    && ./bootstrap.sh ${BOOTSTRAP_OPTIONS} --prefix=${INSTALL_DIR} --with-libraries=$BOOST_LIBRARY_LIST \
    && ./bjam -j$PARALLEL ${BJAM_OPTIONS} --libdir=${INSTALL_DIR}/lib --layout=system link=static threading=single install \
    && touch ${INSTALL_DIR}/boost_install_complete

if [ $? != 0 ] ; then ilog "$SCRIPT: build failed: Terminating..."; exit 1 ; fi

ilog "Cleaning up ${SOURCE_DIR}"
rm -rf ${SOURCE_DIR}

ilog "boost-$VERSION installed successfully"
