<link rel='stylesheet' href='userguide.css' />

Manta User Guide
================

Version: @MANTA_VERSION@

## Introduction

Manta is a structural variant caller for short sequencing reads. It discovers
structural variants of any size and scores these according to a diploid genotype
model. It will also score structural variants according to a somatic model when
separate tumor and normal reads are specified. Both discovery and scoring
incorporate paired-read fragment spanning and split read evidence.

## Capabilities

### Detected variant classes

Manta should be able to detect simple indels.

#### Known Limitations

Manta should not be able to detect the following variant types:

* Non-tandem repeats/amplifications
* Large insertions
    *   The maximum detectable size should correspond to approximately the
  read-pair fragment size, but note that detection power should fall off to
  impractical levels well before this size
* Small inversions
    *  The limiting size is not tested, but in theory detection falls off below
  approx 1kb-500bases. So-called micro-inversions might be detected, but as
  combined insertion/deletion variants.

More general repeat-based limitations exist for all variant types:

* Power to assemble variants to breakend resolution falls to zero as breakend
  repeat length approaches the read size.
* Power to detect any breakend falls to (nearly) zero as the breakend repeat
  length approaches the fragment size.
* The method cannot detect non-tandem repeats

Additionally note that manta does not offer high-level classification of SVs
beyond finding simple DNA-adjacencies.

### Input requirements

Manta requires input sequencing reads to be mapped by an external tool and
provided as input in BAM format.

## SV Output

## Run configuration and Execution

