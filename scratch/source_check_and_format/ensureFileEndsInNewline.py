#!/usr/bin/env python
"""
Treat all arguments as input files and append a newline to the end of each input file if one does not exist.

Usage:
    $0 file1 [file2 [file3]]...

"""


import os.path,sys

if len(sys.argv) < 2 :
    sys.stderr.write(__doc__)
    sys.exit(2)

for arg in sys.argv[1:] :
    assert(os.path.isfile(arg))
    fp=open(arg,"rb+")
    fp.seek(-1,2)
    if fp.read() != "\n" : fp.write("\n")
    fp.close()

