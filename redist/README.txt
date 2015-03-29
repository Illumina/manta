3rd party package modification notes:

boost has been modified to remove some files according to
$ROOT/scratch/make_boost_subset.bash

samtools and htslib hava been modified to remove the test/
directories, in addition to all test and curses requirements
form the Makefiles.

cmake-modules-c99fd3 modified to show git describe --dirty 
