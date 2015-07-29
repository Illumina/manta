#
# This is a simple image used to assist with deploying portable
# binaries, by allowing us to build in a virtual centos 5
# environment.
#
# At present we simply add all of the projects build requirements onto
# the centos5 base image.
#

FROM centos:5

MAINTAINER Chris Saunders <csaunders@illumina.com>

# add standard centos5 packages:
RUN yum install -y bzip2 make gcc gcc-c++ tar wget zlib-devel && \
    yum install -y epel-release && \
    yum install -y cmake28
RUN cd /usr/bin && ln -s cmake28 cmake

# add gcc 4.8 from developer tools v2:
RUN wget http://people.centos.org/tru/devtools-2/devtools-2.repo -O /etc/yum.repos.d/devtools-2.repo && \
    yum install -y devtoolset-2-gcc devtoolset-2-gcc-c++ devtoolset-2-binutils

# Prior to build configuration, set CC/CXX to gcc 4.8:
ENV CC /opt/rh/devtoolset-2/root/usr/bin/gcc
ENV CXX /opt/rh/devtoolset-2/root/usr/bin/g++
