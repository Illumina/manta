#!/usr/bin/env bash

reltoabs() {
  (cd $1; pwd)
}

script_dir=$(reltoabs $(dirname $0))

$script_dir/reFilterVcf.py --minSS 10 |\
$script_dir/minSVSizeFilter.py --minSize 100 |\
awk '/^#/ || /PASS/' |\
$script_dir/inversionFilter.py |\
$script_dir/largeIntrachromFilter.py --maxSize 100000 |\
$script_dir/overlapFilter.py |\
$script_dir/pairSupportFilter.py


