#
# This is a simple image used to assist with deploying portable
# binaries, by allowing us to build in a virtual centos 5
# environment.
#
# At present we simply add all of the projects build requirements onto
# the centos5 base image.
#

FROM astj/centos5-vault

# add standard centos5 packages:
RUN yum install -y bzip2 make gcc gcc-c++ tar wget zlib-devel git && \
    yum install -y epel-release && \
    yum install -y cmake28 && \
    cd /usr/bin && ln -s cmake28 cmake

# add gcc 4.8 from developer tools v2:
RUN wget http://people.centos.org/tru/devtools-2/devtools-2.repo -O /etc/yum.repos.d/devtools-2.repo && \
    yum install -y devtoolset-2-gcc devtoolset-2-gcc-c++ devtoolset-2-binutils

# Prior to build configuration, set CC/CXX to new gcc version:
ENV CC /opt/rh/devtoolset-2/root/usr/bin/gcc
ENV CXX /opt/rh/devtoolset-2/root/usr/bin/g++
