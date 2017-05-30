#!/usr/bin/env python
"""
apply new header to a cpp file
"""

import os
import sys

script_name=os.path.basename(__file__)


if len(sys.argv) != 2 :
    sys.stderr.write("usage: %s new_header < file\n" % (script_name))
    sys.exit(2)
   

header_file = sys.argv[1] 

newstr = "//\n"

for line in open(header_file) :
    newstr += "//"
    if line != "\n" : newstr += " "
    newstr += line

newstr += "//\n"

is_first = True
count = 0

infp=sys.stdin
outfp=sys.stdout

for line in infp :
    if is_first :
        if line.startswith("//") :
            count += 1
            continue
        outfp.write(newstr)
        is_first=False
    outfp.write(line)

if is_first :
    raise Exception("No input!")

