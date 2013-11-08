<link rel='stylesheet' href='userguide.css' />

Manta User Guide
================

Version: @MANTA_VERSION@

## Introduction

Manta is a structural variant caller for short sequencing reads. It is capable
of discovering structural variants of any size and scoring these using both a
diploid genotype model and a somatic model (when separate tumor and normal
samples are specified). Structural variant discovery and scoring incorporate
both paired-read fragment spanning and split read evidence to improve
specificity.

## Capabilities

Primary application is WGS DNA-Seq. Requires mapped short sequencing reads
provided in the BAM file format. At least one BAM file must be submitted for
each sample.

Limited testing of WES. Configure with the `--exome` flag.

Limited testing RNA-Seq. Configure with the `--rna` flag.

### Detected variant classes

Manta is able to detect all classes of variations which can be explained as
novel DNA adjacencies in the genome. Simple insertion/deletion events can be
detected down to a configurable minimum size cutoff (defaulting to 51). All DNA
adjacencies are classified into the following categories based on the breakend
pattern:

* Deletions
* Insertions
* Inversions
* Tandem Duplications
* Interchromosomal Translocations

#### Known Limitations

Manta should not be able to detect the following variant types:

* Non-tandem repeats/amplifications
* Large insertions
    *    The maximum detectable size should correspond to approximately the
  read-pair fragment size, but note that detection power should fall off to
  impractical levels well before this size
* Small inversions
    *   The limiting size is not tested, but in theory detection falls off below
  approx 200bases. So-called micro-inversions might be detected indirectly as
  combined insertion/deletion variants.

More general repeat-based limitations exist for all variant types:

* Power to assemble variants to breakend resolution falls to zero as breakend
  repeat length approaches the read size.
* Power to detect any breakend falls to (nearly) zero as the breakend repeat
  length approaches the fragment size.
* The method cannot detect non-tandem repeats

Note that while manta classifies novel DNA-adjacencies, it does not infer the
higher level constructs implied by the classification. For instance, a variant
marked as a deletion by manta indicates an intrachromosomal translocation of high
confidence with a breakend pattern which is consistent with a deletion, however
there is no higher-order testing of depth and all intersecting adjacencies to
directly infer the SV type.

### Input requirements

Manta requires input sequencing reads to be mapped by an external tool and
provided as input in BAM format.

## Outputs

### SV predictions

All variants are reported in the vcf using  symbolic alleles unless both of these conditions are met:
*   The variant can be expressed as a combination of inserted and deleted sequence.
*   The deletion or insertion length is not 1000 or greater

### Statistics

## Run configuration and Execution

