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

## Method Overview

Manta works by dividing the structural variant (SV) discovery process into two
primary steps: (1) scanning the genome to find SV associated regions and (2)
analysis, scoring and output of SVs found in such regions.

1. **Build SV association graph** In this step the entire genome is scanned to
discover evidence of possible SVs and large indels. This evidence is enumerated
into a graph with edges connecting all regions of the genome which have a
possible SV association. Edges may connect two different regions of the genome
to represent evidence of a long-range association, or an edge may connect a
region to itself to capture a local indel/small SV association. Note that these
associations are more general than a specific SV hypothesis, in that zero to
many SV candidates may be found on one edge, although we would typically expect
only one or two candidates per edge.

2. **Analyze graph edges to find SVs** The second step is to analyze individual
graph edges or groups of highly connected edges to discover and score any SVs
associated with the edge(s). This process entails several substeps, including
inference of SV candidates associated with the edge, attempted assembly of the
SVs breakends, scoring and filtration of the SV under various biological models
(currently diploid germline and somatic), and finally, output to VCF.

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

At configuration time, at least one bam file must be provided for the normal
sample and optionally can also read one or more bam files for a tumor sample.
If multiple bams are provided for the normal or tumor sample, these are treated
as a merged single normal or tumor.

The following limitations exist on the input BAMs provided to Manta:

* Alignments cannot contain the "=" character in the SEQ field.
* Alignments cannot use sequence the match/mismatch ("="/"X") CIGAR notation
* RG (read group) tags in the BAMs are ignored -- each BAM must represent one
  sample.
* Alignments with basecall quality values greater than 70 are rejected (these
  are not supported on the assumption that this indicates an offset error)

## Outputs

### SV predictions

The primary manta outputs are a set of [VCF 4.1][1] files, found in
`${RUNFOLDER}/results/variants`. Currently there are at least two vcf files
created for any manta run, and a third somatic vcf is produced when tumor input
is provided. These files are:

* __candidateSV.vcf.gz__ Unscored SV and indel candidates. Only a minimal
  amount of supporting evidence is required for a SV to entered as a candidate.
  An SV or indel must be a candidate to be considered for scoring, therefor if
  an SV cannot appear in the other VCF outputs if it is not present here. Note
  that this file includes indels down to a very small size (10 by default)
  intended to be passed on to a small variant caller without scoring by manta
  itself (by default manta scoring starts at size 51).
* __diploidSV.vcf.gz__ SVs and indels scored and genotyped under a diploid
  model for the normal sample. The scores in this file do not reflect any
  information in the tumour bams
* __somaticSV.vcf.gz__ SVs and indels scored under a somatic variant model.
  This file will only be produced if at least one tumour bam argument is
  supplied during configuration

All variants are reported in the vcf using  symbolic alleles unless both of these conditions are met:

* The variant can be expressed as a combination of inserted and deleted sequence.
* The deletion or insertion length is not 1000 or greater

### Statistics

Additional secondary output is provided in `${RUNFOLDER}/results/stats`

* __alignmentStatsSummary.txt__ fragment length quantiles for each input bam
* __svLocusGraphStats.tsv__ statistics pertaining to the SV locus graph

## Run configuration and Execution

Manta is run in a two step procedure: (1) configuration and (2) workflow
execution. The configuration step is used to specify the input data and any
options pertaining to the variant calling methods themselves  The execution
step is used to specify any parameters pertaining to _how_ manta is executed
(such as the total number of cores or SGE nodes over which the jobs should be
parallelized). The second execution step can also be interrupted and restarted
without changing the final result of the workflow.

### Configuration

Example:

```
${MANTA_INSTALL_DIR}/bin/configManta.py \
--normalBam NA12878_S1.bam \
--referenceFasta hg19.fa \
--runDir ${ANALYSIS_RUN_DIR}
```

### Execution

Example execution on an SGE cluster:

```
${ANALYSIS_RUN_DIR}/runWorkflow.py -m sge -j 36
```


[1]: http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
