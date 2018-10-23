# External packages

This directory contains external packages for building manta.

## Package modification notes

### boost

To reduce size, certain files have been removed from the
full boost source distribution acccording to:

    ${ROOT_PATH}/scratch/make_boost_subset.bash

### htslib/samtools

To reduce size, both packages have been modified to remove the test
directories and test references in the makefiles. The copy of htslib
in samtools has been removed.

### cmake-modules

cmake-modules-c99fd3 modified to show git describe --dirty

