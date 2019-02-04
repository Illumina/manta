# Manta Developer Guide - Debugging a single SV in Manta

[Developer Guide Home](README.md)

## Table of Contents

[//]: # (BEGIN automated TOC section, any edits will be overwritten on next source refresh)

* [Summary](#summary)
* [Scenario 0: Debug SV call in both graph creation and SV candidate generation steps](#scenario-0-debug-sv-call-in-both-graph-creation-and-sv-candidate-generation-steps)
* [Scenario 1 : Debug gold-standard SV call which is already covered by an edge in the SVLocus graph](#scenario-1--debug-gold-standard-sv-call-which-is-already-covered-by-an-edge-in-the-svlocus-graph)
  * [Step 1: Identify the graph edge corresponding to the SV](#step-1-identify-the-graph-edge-corresponding-to-the-sv)
    * [S1 - Step 1A : Query Node in the SV Locus graph](#s1---step-1a--query-node-in-the-sv-locus-graph)
    * [S1 - Step 1B : Query Node in the SV Locus graph](#s1---step-1b--query-node-in-the-sv-locus-graph)
    * [S1 - Step 1C : Query Locus in the SV Locus graph](#s1---step-1c--query-locus-in-the-sv-locus-graph)
    * [S1 - Step 2 : Run candidate generation on specific SV locus or specific SV locus edge](#s1---step-2--run-candidate-generation-on-specific-sv-locus-or-specific-sv-locus-edge)
* [Debugging Infrastructure TODO](#debugging-infrastructure-todo)

[//]: # (END automated TOC section, any edits will be overwritten on next source refresh)


## Summary

Manta tries to reduce intermediate I/O as much as possible, which is helpful in production, but not so for development/debug. As debuging cases come up, speciallized debug/localized running modes have been added to certain Manta tools, and more will certainly be added. This page documents the workflows that exist – especially those that facilitate debugging a single known FN or FP case. We also note workflows which would be useful and should be added in the future. For related debug discussions involving properties of an entire run, see the related page: [Debugging a full Manta run](debugFullRun.md)

## Scenario 0: Debug SV call in both graph creation and SV candidate generation steps

Manta workflow has the ability to run one to many sub-segments of the genome, which can accelerate debugging for by setting the region(s) to cover the breakends of any SVs of interest.

To enable a regional build, simply specify the additional flag "--region" at configuration-time and provide a region in samtools format. This region flag can be repeated multiple times to (for instance) cover both breakend regions of a translocation. Documentation for this and other hidden debug/development options appears when you provide the '--allHelp' flag as follows:

    ${MANTA_INSTALL_ROOT}/bin/configManta.py --allHelp

A full example of a rapid debug region run is:

    ${MANTA_INSTALL_ROOT}/bin/configManta.py --normalBam myBam.bam --region chr2:19000-20000

An example debug build for a translocation is:

    ${MANTA_INSTALL_ROOT}/bin/configManta.py --normalBam myBam.bam --tumorBam myTumorBam.bam --region chr2:19000-20000 --region chr20:1000-2000

## Scenario 1 : Debug gold-standard SV call which is already covered by an edge in the SVLocus graph

If a known clean SV is either missing from the output, or we would like to repeatedly/quickly run just this one SV during development of a new feature, there is limited support to run a specific SV (actually a specific disjoint subgraph of the SVgraph – often for gold-standard calls, the call is the only member of a disjoint sub-graph).

### Step 1: Identify the graph edge corresponding to the SV

Every edge in the breakend graph has a locus index (each locus is a connected sub-graph), a node index (nodes are numbered from 0 within each locus) and a second node index (this case equal the first node index for a self-edge). It is possible that there is no edge corresponding to your SV, in which case the FN problem extends all the way back to graph candidate generation, you can determine if this is the case in step 1B below.

The following subsections of step1 describe how to extract these graph indices for an SV which is already being printed out at least to the candidate vcf file (Step 1A – easier case), or for the general case (Step 1B, a bit more involved).

#### S1 - Step 1A : Query Node in the SV Locus graph

The following record from a Manta candidate VCF shows the locus-index<sup>**A**</sup>, node1-index<sup>**B**</sup> and node2-index<sup>**C**</sup> , per the corresponding superscripts applied to each field:

> chrX    70129275        MantaINV:572680<sup>**A**</sup>:1<sup>**B**</sup>:3<sup>**C**</sup>:0:0 A       <INV>   .       .       END=70140033;SVTYPE=INV;SVLEN=-10758;UPSTREAM_PAIR_COUNT=29;DOWNSTREAM_PAIR_COUNT=29;PAIR_COUNT=29;CIPOS=0,7;CIEND=-7,0;HOMLEN=7;HOMSEQ=ACATGGA

If these are available for your SV of interest, then you can skip to Step 2 below.

#### S1 - Step 1B : Query Node in the SV Locus graph

Hypothesis generation starts from a graph of connected genomic regions called the SV locus graph. After manta has been run, a binary file containing the graph can be found in:

    ${MANTA_RUN_DIR}/workspace/svLocusGraph.bin

The most useful way to examine the content of the graph is by dumping all regions/edges which overlap a specific genomic region, for example:

    ${MANTA_INSTALL_ROOT}/libexec/DumpSVLoci --graph-file svLocusGraph.bin --region chr2:2000000-2000100

Making such a query over the region where either of a gold-standard SV's breakends are expected should reveal if the region was even detected as SV associated.

For example, to check for the first breakend of COSMIC deletion COST17172, the following query can be run on the HCC1187 graph file:

    ${MANTA_INSTALL_ROOT}/libexec/DumpSVLoci --graph-file svLocusGraph.bin --region chrX:154205937-154205937

..which yields:

```
SVNode LocusIndex:NodeIndex : 13716:0
LocusNode: count: 50 GenomeInterval: 23:[154205488,154206145) n_edges: 1 out_count: 50 evidence: [154205542,154205948)
        EdgeTo: 1 out_count: 50
```

This shows that the locus index is **13716**

#### S1 - Step 1C : Query Locus in the SV Locus graph

In the output above (end of Step 1), we gain access to the locus index number of the disjoint subgraph which covers (at least) the COST17172 breakends This locus index is 13716 from the first line of the output above). Using this index number we can query the graph file for this disjoint subgraph to see all components connected to the initial region:

    ${MANTA_INSTALL_ROOT}/libexec/DumpSVLoci --graph-file svLocusGraph.bin --locus-index 13716

Result:

```
LOCUS BEGIN INDEX 13716
NodeIndex: 0 LocusNode: count: 50 GenomeInterval: 23:[154205488,154206145) n_edges: 1 out_count: 50 in_count: 55 evidence: [154205542,154205948)
        EdgeTo: 1 out_count: 50 in_count: 55
NodeIndex: 1 LocusNode: count: 55 GenomeInterval: 23:[154226371,154226808) n_edges: 1 out_count: 55 in_count: 50 evidence: [154226586,154226908)
        EdgeTo: 0 out_count: 55 in_count: 50
LOCUS END INDEX 13716
```

This shows us that the deletion is supported by one edge in the SVLocus graph. We can also see that the edge connects node index 0 to node index 1. Together with the locus index, we have all graph indices required to run only the edge corresponding to our SV of interest.


#### S1 - Step 2 : Run candidate generation on specific SV locus or specific SV locus edge

We can also run Candidate SV generation for an entire locus, all edges which connect to one node in a locus or only one edge in a locus. To do this we extract the GenerateSVCandidates command from a manta run and add the `--locus-index ARG` flag. ARG can be:

option | example | description
------ | ------- | -----------
--locus-index LocusIndex | --locus-index 13716 | Run all edges in the specified locus only
--locus-index LocusIndex:NodeIndex | --locus-index 13716:0 | Run all edges which connect to the specified node in the specified locus
--locus-index LocusIndex:NodeIndex1:NodeIndex2 | --locus-index 13716:0:1 | Run the single edge connecting the two specified nodes (or the self-edge if NodeIndex1==NodeIndex2)

An example full command-line is:

> ${MANTA_INSTALL_ROOT}/libexec/GenerateSVCandidates --align-stats /tmp/manta_test/assemble_test/testAssm3/workspace/alignmentStats.xml --graph-file /tmp/manta_test/assemble_test/testAssm3/workspace/svLocusGraph.bin --ref genome.fa --candidate-output-file /tmp/manta_test/assemble_test/testAssm3/workspace/svHyGen/candidateSV.0103.vcf --somatic-output-file /tmp/manta_test/assemble_test/testAssm3/workspace/svHyGen/somaticSV.0103.vcf --chrom-depth /tmp/manta_test/assemble_test/testAssm3/workspace/chromDepth.txt --align-file sorted.bam --tumor-align-file tsorted.bam **--locus-index 13716:0:1** **--verbose**

The additional `--locus-index ARG` command is highlighted, together with the new `--verbose` option. In the example, SV generation runs for the specified edge "13716:0:1" only. This makes it easier to run modifications of the SV generator with various types of verbose debugging outputs, etc...  To get started in this direction, the example includes the --verbose option to provide some quick high level logging without recompiling – for many problems more specific/noising debug output will have to be compiled in as part of a follow-up step.

## Debugging Infrastructure TODO

The above process could be more streamlined, especially for cases where an SV is part of a large disjoint subgraph. New features:
* GenerateSVCandidates should accept a --region1 and --region2 argument, and only call SVGraph edges connecting those two regions.


