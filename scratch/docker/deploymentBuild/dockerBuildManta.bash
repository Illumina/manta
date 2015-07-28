#!/usr/bin/env bash

set -o nounset
set -o xtrace

buildLabel="manta-unknown-version"
if ! [ $# == 1 ]; then
    cat <<END
usage: $0 [binary_distro_prefix]

Create centos5 binary distro tarball using docker
END
    exit 2
fi

buildLabel="$1"


rel2abs() {
  cd $1 && pwd -P
}

scriptDir=$(rel2abs $(dirname $0))

echo scriptDir $scriptDir

# in dockerfile directory:
tag="deployment:MantaBuilderImage"
sudo docker build -t $tag $scriptDir

# in scratch
#unpack src tarball and cd into tarball root

dmount=/builder

installScript=installManta.bash

cat << ENDE >| $installScript 
set -o errexit
set -o nounset

# build
mkdir -p build
cd build
$dmount/src/configure --prefix=$dmount/install --jobs=2
make -j2 install

# make tarball
cd $dmount
mv install $buildLabel
tar -cj $buildLabel -f $buildLabel.tar.bz2
rm -rf $buildLabel
chmod 777 $buildLabel.tar.bz2
ENDE

sudo docker run -v $(pwd):$dmount -t $tag bash $dmount/$installScript

