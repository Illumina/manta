#!/usr/bin/env bash

set -o nounset


rel2abs() {
    (cd $1 && pwd -P)
}

docdir=$(rel2abs $(dirname $0))
builddir=$docdir/build

mkdir -p $builddir

mname=methods

do_latex_cmds() {
  file=$1
  latex $file
  bibtex $file
  latex $file
  latex $file
  dvipdf $file
}

for mm in $mname; do
(
cd $builddir
cp ../packages/* .
ln -sf $docdir/$mm.tex
ln -sf $docdir/$mname.bib
ln -sf $docdir/figures
do_latex_cmds $mm
mv $mm.pdf $docdir 
mv $mm.log $docdir
rm -f $mm.*
)
done


