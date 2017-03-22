# Manta Developer Guide - Breakend Graph

[Developer Guide Home](README.md)

## Table of Contents

[//]: # (BEGIN automated TOC section, any edits will be overwritten on next source refresh)

* [Summary](#summary)
* [Querying the graph](#querying-the-graph)

[//]: # (END automated TOC section, any edits will be overwritten on next source refresh)


## Summary

The breakend graph is a central intermediate  in the manta workflow. It contains the following inforation:
* "Nodes" – this are contiguous regions of the genome. These are associated with one or more breakends.
  * Evidence range: Each node has an additional chromosomal range where it's read evidence was originally aligned. (if the node is built from indirect evidence only, then evidence range should be the same as the node range.
  * Edges between nodes – Every node with an edge has a return pointer in a properly formatted SV graph, but the evidence count on an out-edge represents that evidence for the edge was observed at the node the edge is coming from.
* "Loci" – A locus is a disjoint subgraph

## Querying the graph
The graph can be queried as follows:

* Given a genomic region, write out all nodes which intersect that region
* List a specific locus id (such as one identified when querying a region)
* List the whole graph

Example:

    ${MANTA_INSTALL_ROOT}/libexec/DumpSVLoci --graph-file foo --region chr10:1000000-1001000

You can also get summary metrics from the graph:
* List of node count, edge count, observation count, etc. for every locus
* Summary of total reads used and total reads cleaned out as noise edges.
