#!/usr/bin/env bash

set -o nounset
set -o errexit


rel2abs() {
    (cd $1 && pwd -P)
}

docdir=$(rel2abs $(dirname $0))
builddir=$docdir/build

mkdir -p $builddir

mname=methods

latexCmd() {
  latex -halt-on-error -interaction=nonstopmode $1
}

do_latex_cmds() {
  file=$1
  latexCmd $file
  bibtex $file
  latexCmd $file
  latexCmd $file
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


