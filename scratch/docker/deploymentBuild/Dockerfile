#
# This is a simple image used to assist with deploying portable
# manta binaries, by allowing us to build in a virtual centos 5
# environment.
#
# At present we simply add all of manta's build requirements onto 
# the centos5 base image.
#

FROM centos:5

MAINTAINER Chris Saunders <csaunders@illumina.com>

# add standard centos5 packages:
RUN yum install -y tar wget bzip2 make gcc gcc-c++ zlib-devel

# add cmake 2.8
RUN yum install -y epel-release
RUN yum install -y cmake28
CMD cd /usr/bin && ln -s cmake28 cmake

# add gcc 4.8 from developer tools v2:
RUN wget http://people.centos.org/tru/devtools-2/devtools-2.repo -O /etc/yum.repos.d/devtools-2.repo
RUN yum install -y devtoolset-2-gcc devtoolset-2-gcc-c++ devtoolset-2-binutils

# Prior to build configuration, set CC/CXX to gcc 4.8:
ENV CC /opt/rh/devtoolset-2/root/usr/bin/gcc
ENV CXX /opt/rh/devtoolset-2/root/usr/bin/g++

