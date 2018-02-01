#!/usr/bin/env python2
#
# Manta - Structural Variant and Indel Caller
# Copyright (c) 2013-2018 Illumina, Inc.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#

"""
Very simple TOC from markdown solution.

The script looks for the first occurance of /^## Table of Contents/ and inserts a
table of contents immediately below that marker, replacing anything occuring between
the 'table of contents' line and the first empty line that follows.
"""

import re,sys


def getOptions() :

    from optparse import OptionParser

    usage = "usage: %prog [options] < in.md > out.md"
    description="""
reset SOMATICSCORE filter to a new value
"""
    parser = OptionParser(usage=usage)

    parser.add_option("--depth", type="int", default=3,
                      help="TOC depth (relative range from lowest to highest header level)")

    (opt,args) = parser.parse_args()

    if len(args) != 0 :
        parser.print_help()
        sys.exit(2)

    # validate input:
    if opt.depth < 1:
        raise Exception("Invalid depth value: %i" % (opt.minSS))

    return (opt,args)


def getHeaderLevel(line):
    index=0
    while line[index] == "#" :
        index += 1
    if line[index] == " " : return index
    return 0


class Constants:
    tocPattern = re.compile("^#+[ \t]*[Tt]able[ \t]*[Oo]f[ \t]*[Cc]ontents")

    # characters which must be removed to create an anchor in GFM:
    rmLinkChars = "[^0-9a-z -]"

    # this works as a semi-inline comment on github, but may not work on other md parsers:
    warnString = "\n[//]: # (%s automated TOC section, any edits will be overwritten on next source refresh)\n\n"


def main() :
    
    (opt,args) = getOptions()
    
    infp=sys.stdin
    outfp=sys.stdout
    
    # true when current line is part of a previous table of contents or body text
    isBuilding=False

    # true when the body text of the document is being analyzed to build the new table of contents
    isScanning=False
    isValidHeaderZone=True
    lineBuffer=[]
    tocInfo=[]
    for line in infp :
        if not isBuilding :
            outfp.write(line)
            if Constants.tocPattern.match(line) :
                isBuilding=True
        elif not isScanning:
            sline = line.strip()
            if not ((sline == "") or (sline[0] == "*") or (sline[0] == "[")) : isScanning=True

        if isBuilding and isScanning :
            if line.startswith("```") :
                isValidHeaderZone = not isValidHeaderZone

            if isValidHeaderZone :
                level = getHeaderLevel(line)
                if level > 0 :
                    tocInfo.append([level,line.strip().split(None, 1)[1]])
 
            lineBuffer.append(line)
    
    # build TOC:
    lowDepth=None
    for tocEntry in tocInfo :
        if (lowDepth is not None) and (tocEntry[0] >= lowDepth) : continue
        lowDepth = tocEntry[0]

    if isBuilding :
        outfp.write(Constants.warnString % ("BEGIN"))

    linkTags=set()

    for tocEntry in tocInfo :
        relDepth=tocEntry[0]-lowDepth
        if relDepth >= opt.depth : continue
        for _ in range(relDepth) :
            outfp.write("  ")

        # Replicate the github anchor generation method. Summary below copied from @TomOnTime gist comments:
        #
        # The code that creates the anchors is here:
        #   https://github.com/jch/html-pipeline/blob/master/lib/html/pipeline/toc_filter.rb
        #
        # 1) downcase the string
        # 2) remove anything that is not a letter, number, space or hyphen (see the source for how Unicode is handled)
        # 3) change any space to a hyphen.
        # 4) If that is not unique, add "-1", "-2", "-3",... to make it unique
        #
        linkTagRoot=tocEntry[1].lower()
        linkTagRoot=re.sub(Constants.rmLinkChars,"",linkTagRoot).replace(' ', '-')

        linkTag=linkTagRoot
        linkIndex=1
        while linkTag in linkTags :
            linkTag = linkTagRoot + "-" + str(linkIndex)
            linkIndex += 1
        linkTags.add(linkTag)

        outfp.write("* [%s](#%s)\n" % (tocEntry[1],linkTag))
        
    if isBuilding :
        outfp.write(Constants.warnString % ("END"))
        outfp.write("\n")

    for line in lineBuffer :
        outfp.write(line)


main()
