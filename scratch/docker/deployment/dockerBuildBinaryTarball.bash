#!/usr/bin/env bash

set -o nounset
set -o xtrace

buildLabel=
if ! [ $# != 3 ]; then
    cat <<END
Create binary distro tarball using docker.

usage: $0 distro_root_directory binary_distro_prefix

  distro_root_directory - root directory of the package git repository
  binary_distro_prefix - filename prefix of tarball and directory name
                         used when tarball is unpacked
END
    exit 2
fi

rootDir="$1"
buildLabel="$2"


rel2abs() {
  cd $1 && pwd -P
}

builderImage=centos5PlusGcc49FromSrc
scriptDir=$(rel2abs $(dirname $0))
#rootDir=$(rel2abs $scriptDir/../../..)
rootDir=$(rel2abs $rootDir)

# check that rootDir conatins expected files:
if ! [ -f $rootDir/configure ]; then
    echo "Can't find package configure script. Expected location is '$rootDir/configure'" 2>&1
    exit 1
fi


# in dockerfile directory:
tag="deployment:$builderImage"
sudo docker build -t $tag $scriptDir/$builderImage

# in scratch
#unpack src tarball and cd into tarball root

dmount=/builder

installScriptFilename=buildBinaryTarball.bash
installScriptPath=$rootDir/$installScriptFilename

trap "{ rm -f $installScriptPath; }" EXIT

cat << ENDE >| $installScriptPath
set -o errexit
set -o nounset

# build
mkdir -p build
cd build
$dmount/configure --prefix=$dmount/install --jobs=2
make -j2 install

# make tarball
cd $dmount
mv install $buildLabel
tar -cj $buildLabel -f $buildLabel.tar.bz2
rm -rf $buildLabel
chmod 777 $buildLabel.tar.bz2
ENDE

sudo docker run -v $rootDir:$dmount -t $tag bash $dmount/$installScriptFilename
mv $rootDir/$buildLabel.tar.bz2 .

