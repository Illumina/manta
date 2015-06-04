#!/usr/bin/env bash

set -o errexit
set -o nounset
set -o pipefail
set -o xtrace


package_name=manta

pname_root=""
if [ $# -gt 1 ]; then
    cat <<END
usage: $0 [tarball_rootname]

This script makes the $package_name source tarball

- script assumes that it is located in the release checkout version
- the tarball is written to the caller's working directory
END
    exit 2
elif [ $# == 1 ]; then
    pname_root=$1
fi

rel2abs() {
    cd $1; pwd -P
}

thisdir=$(rel2abs $(dirname $0))
outdir=$(pwd)

cd $thisdir
# don't bother passing --dirty to git describe, the
# archive step below will remove any changes from
# HEAD anyway:
gitversion=$(git describe --match "v[0-9]*" | sed "s/^v//")
if [ $? != 0 ]; then
    echo "ERROR: 'git describe' failed" 1>&2
    exit 1
fi

if [ "$pname_root" == "" ]; then
    pname_root=${package_name}-$gitversion
fi

pname=$outdir/$pname_root

if [ -d $pname ]; then
    echo "ERROR: directory already exists: $pname"
    exit 1
fi
mkdir -p $pname

cd ..
git archive --prefix=$pname_root/ HEAD: | tar -x -C $outdir

# substitute certain build values for distribution:
tmp_file=$(mktemp)
cml=$pname/src/CMakeLists.txt
awk '
{
    if (/^set *\(DEVELOPER_MODE/) printf "set (DEVELOPER_MODE false)\n";
    else print;
}' $cml >| $tmp_file
mv $tmp_file $cml

cml=$pname/src/cmake/getBuildTimeConfigInfo.cmake
awk -v gver=$gitversion '
{
    if      (/^set *\(WORKFLOW_VERSION/) printf "set (WORKFLOW_VERSION \"%s\")\n",gver;
    else if (/Detected workflow version/) printf "";
    else print;
}' $cml >| $tmp_file
mv $tmp_file $cml

rme=$pname/README.md
awk -v gver=$gitversion '
{
    if      (!ver && !NF) { printf "\nVersion: %s\n\n",gver; ver=1}
    else if ($0~/_NOT_ part of/) a=1;
    else if ($0~/\[Build Status\]/) b=1;
    else if (b && !NF) b=0;
    else print;
}' $rme >| $tmp_file
mv $tmp_file $rme

# tar it up:
(
cd $outdir
tar -f $pname_root.release_src.tar.bz2 -cj $pname_root
)

rm -rf $pname
