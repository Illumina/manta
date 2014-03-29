#!/usr/bin/env bash

reltoabs() {
  (cd $1; pwd)
}

script_dir=$(reltoabs $(dirname $0))

awk '/^#/ || /PASS/' |\
$script_dir/inversionFilter.py |\
$script_dir/largeIntrachromFilter.py --maxSize 100000 |\
$script_dir/overlapFilter.py 


