# Manta Developer Guide - VCF ID field

[Developer Guide Home](README.md)

Manta includes an identifier (VCF `ID` field) for every VCF record in its output. This ID is guaranteed to be unique within any VCF output by manta, as requierd by the VCF spec.

The `ID` field provides information linking the variant call to the SV association graph. It can be useful for certain debugging procedures.

The identifier has two minor variants, one used for translocation breakends and one variant used for all other VCF records. An example ID for a translocation record is:

```
MantaBND:5862:0:1:0:0:0:1
```

An example for other record types is:

```
MantaDEL:47029:3:9:0:0:0
```

We describe the first example ID by describing each sub-field when the ID is split on the colon character. Note the [breakend graph description](breakendGraph.md) may be helpful for reference here.

index | ID component name | Value from above example | Description
----- | ----------------- | ------------------------ | -----------
1 | Label | MantaBND | Simple text label appending "Manta" to the SVType. The ID is still unique if this component is removed.
2 | LocusID | 5862 | Index of the SV breakend graph **locus**. Each locus is a disjoint subgraph of the full breakend graph.
3 | Node1ID | 0 | Index of the first SV breakend graph **node** forming the graph edge used to discover this variant.
4 | Node2ID | 1 | Index of the second SV breakend graph **node** forming the graph edge used to discover this variant. If Node1ID==Node2ID this is a self-edge.
5 | CandidateID | 0 | Each graph edge is analyzed for evidnece of specific SV or indel candidates. This index provides the index of the source candidate among all candidates associated with this edge.
6 | AssemblyID | 0 | For each candidate multiple contigs/paths are extracted from the assembly graph, this describes the index of the path used to generate this candidate or "0" for an IMPRECISE variant
7 | SegmentID | 0 | Multiple small variants can be extracted from each assembly contig/path alignment, this describes the alingment segment index used to produce the candidate variant. This index can only be non-zero for small indels.
8 | BNDID | 1 | This sub-field is only used in BND records. BND records are different than other VCF record types in that one variant is represented by two breakends. This component is set to either 0 or 1 to indicate the breakend number of the variant.
