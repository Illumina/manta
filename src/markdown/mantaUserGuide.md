<link rel='stylesheet' href='userGuide.css' />

Manta User Guide
================

Version: @MANTA_VERSION@

<script src="tableOfContents.js"></script>

## Introduction

Manta is a structural variant caller for short sequencing reads. It is capable
of discovering structural variants of any size and scoring these using both a
diploid genotype model and a somatic model (when separate tumor and normal
samples are specified). Structural variant discovery and scoring incorporate
both paired-read fragment spanning and split read evidence.

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
associations are more general than a specific SV hypothesis, in that many SV
candidates may be found on one edge, although typically only one or two
candidates are found per edge.

2. **Analyze graph edges to find SVs** The second step is to analyze individual
graph edges or groups of highly connected edges to discover and score SVs
associated with the edge(s). This substeps of this process include
inference of SV candidates associated with the edge, attempted assembly of the
SVs breakends, scoring and filtration of the SV under various biological models
(currently diploid germline and somatic), and finally, output to VCF.

## Capabilities

Manta is capable of detecting all structural variant types which are
identifiable in the absence of copy number analysis and large-scale de-novo
assembly. Detectable types are enumerated further below.

For each structural variant and indel, Manta attempts to align the breakends to
basepair resolution and report the left-shifted breakend coordinate (per the [VCF 4.1][1]
SV reporting guidelines), together with the any breakend homology sequence
and/or inserted sequence between the breakends. It is often the case that the
assembly will fail to provide a confident explanation of the data -- in such
cases the variant will be reported as `IMPRECISE`, and scored according the
paired-end read evidence alone.

The sequencing reads provided as input to Manta are expected to be from a
paired-end sequencing assay which results in an "innie" orientation between the
two reads of each sequence fragment, each presenting a read from the outer edge of
the fragment insert inward.

Manta is primarily tested for whole genome DNA-Seq experiments of single
diploid samples or subtractive analysis of a matched tumor/normal sample pair.
There has been limited testing in support of other cases:

* For exome or other targeted sequencing experiments, the workflow can be
configured with the `--exome` flag to set filtration levels more appropriate for
this case.
* RNA-Seq analysis can be configured with the `--rna` flag to also adjust filtration
levels and take other RNA-specific filtration and intron handling steps.


### Detected variant classes

Manta is able to detect all variation classes which can be explained as
novel DNA adjacencies in the genome. Simple insertion/deletion events can be
detected down to a configurable minimum size cutoff (defaulting to 8). All DNA
adjacencies are classified into the following categories based on the breakend
pattern:

* Deletions
* Insertions
    * Fully-assembled insertions
    * Partially-assembled 'inferred' insertions
* Inversions
* Tandem Duplications
* Interchromosomal Translocations

### Known Limitations

Manta should not be able to detect the following variant types:

* Dispersed duplications
* Most expansion/contraction variants of a reference tandem repeat
* Small inversions
    * The limiting size is not tested, but in theory detection falls off
  below ~200bases. So-called micro-inversions might be detected indirectly as
  combined insertion/deletion variants.
* Fully-assembled large insertions
    * The maximum fully-assembled insertions size should correspond to
  approximately twice the read-pair fragment size, but note that power to fully
  assemble the insertion should fall off to impractical levels before this
  size
    * Note that manta does detect and report very large insertions when the breakend
  signature of such an event is found, even though the inserted sequence cannot
  be fully assembled.  

More general repeat-based limitations exist for all variant types:

* Power to assemble variants to breakend resolution falls to zero as breakend
  repeat length approaches the read size.
* Power to detect any breakend falls to (nearly) zero as the breakend repeat
  length approaches the fragment size.

Note that while manta classifies novel DNA-adjacencies, it does not infer the
higher level constructs implied by the classification. For instance, a variant
marked as a deletion by manta indicates an intrachromosomal translocation with a
deletion-like breakend pattern, however there is no test of depth, b-allele
frequency or intersecting adjacencies to directly infer the SV type.

## Input requirements

The sequencing reads provided as input to Manta are expected to be from a
paired-end sequencing assay with an "innie" orientation between the two reads of
each DNA fragment, each presenting a read from the outer edge of the fragment
insert inward.

Manta can tolerate non-paired reads in the input, so long as sufficient paired-end
reads exist to estimate the paired fragment size distribution. Non-paired reads
will still be used in discovery, assembly and split-read scoring if their alignment
(or SA tag split alignments) support a large indel or SV, or mismatch/clipping 
suggests a possible breakend location.

Manta requires input sequencing reads to be mapped by an external tool and
provided as input in BAM format.

At configuration time, at least one bam file must be provided for the normal
sample. A matched tumor sample can optionally be provided as well. If multiple
bams are provided for the normal or tumor sample, these are merged and treated
as a single normal or tumor sample.

The following limitations exist on the input BAMs provided to Manta:

* Alignments cannot contain the "=" character in the SEQ field.
* Alignments cannot use the sequence match/mismatch ("="/"X") CIGAR notation
* RG (read group) tags in the BAMs are ignored -- each BAM must represent one
  sample.
* Alignments with basecall quality values greater than 70 are rejected (these
  are not supported on the assumption that this indicates an offset error)

## Outputs

### Structural Variant predictions

The primary manta outputs are a set of [VCF 4.1][1] files, found in
`${RUNFOLDER}/results/variants`. Currently there are at least two vcf files
created for any manta run, and a third somatic vcf is produced when tumor input
is provided. These files are:

* __diploidSV.vcf.gz__
    * SVs and indels scored and genotyped under a diploid model for the normal
  sample. The scores in this file do not reflect any information in the tumor
  bams
* __somaticSV.vcf.gz__
    * SVs and indels scored under a somatic variant model. This file
  will only be produced if at least one tumor bam argument is supplied during
  configuration
* __candidateSV.vcf.gz__
    * Unscored SV and indel candidates. Only a minimal amount of supporting
  evidence is required for an SV to be entered as a candidate. An SV or indel
  must be a candidate to be considered for scoring, therefore an SV cannot
  appear in the other VCF outputs if it is not present in this file. Note that
  by default this file includes indels down to a very small size (>= 8 bases).
  These are intended to be passed on to a small variant caller without scoring
  by manta itself (by default manta scoring starts at size 51).
* __candidateSmallIndels.vcf.gz__
    * Subset of the candidateSV.vcf.gz file containing only simple insertion and
  deletion variants of size 50 or less. Passing this file to a small variant caller
  like strelka or starling (Isaac Variant Caller) will provide continuous
  coverage over all indel sizes when the small variant caller and manta outputs are
  evaluated together. Alternate small indel candidate sets can be parted out of the
  candidateSV.vcf.gz file if this candidate set is not appropriate.

All variants are reported in the vcf using symbolic alleles unless they are classified 
as a small indel, in which case full sequences are provided for the vcf `REF` and `ALT`
allele fields. A variant is classified as a small indel if all of these criteria are met:

* The variant can be entirely expressed as a combination of inserted and deleted sequence.
* The deletion or insertion length is not 1000 or greater.
* The variant breakends and/or the inserted sequence are not imprecise.

When vcf records are printed in the small indel format, they will also include
the `CIGAR` INFO tag describing the combined insertion and deletion event.

### Statistics

Additional secondary output is provided in `${RUNFOLDER}/results/stats`

* __alignmentStatsSummary.txt__
    * fragment length quantiles for each input bam
* __svLocusGraphStats.tsv__
    * statistics pertaining to the SV locus graph

## Run configuration and Execution

Manta is run in a two step procedure: (1) configuration and (2) workflow
execution. The configuration step is used to specify the input data and any
options pertaining to the variant calling methods themselves. The execution
step is used to specify any parameters pertaining to _how_ manta is executed
(such as the total number of cores or SGE nodes over which the jobs should be
parallelized). The second execution step can also be interrupted and restarted
without changing the final result of the workflow.

### Configuration

The workflow is configured with the script: `${MANTA_INSTALL_DIR}/bin/configManta.py`
. Running this script with no arguments will display all standard configuration
options to specify input BAM files, the reference sequence and the output run folder.
Note that all input BAMs and reference sequence must contain the same chromosome names
in the same order. Manta's default settings assume a whole genome DNA-Seq analysis, but there
are configuration options for exome/targeted sequencing analysis in addition to RNA-Seq.

Single Sample Analysis -- Example Configuration:

```
${MANTA_INSTALL_DIR}/bin/configManta.py \
--normalBam NA12878_S1.bam \
--referenceFasta hg19.fa \
--runDir ${ANALYSIS_RUN_DIR}
```

Tumor Normal Analysis -- Example Configuration:

```
${MANTA_INSTALL_DIR}/bin/configManta.py \
--normalBam HCC1187BL.bam \
--tumorBam HCC1187C.bam \
--referenceFasta hg19.fa \
--runDir ${ANALYSIS_RUN_DIR}

```

On completion, the configuration script will create the workflow run script `${ANALYSIS_RUN_DIR}/runWorkflow.py`
. This can be used to run the workflow in various parallel compute modes per the
instructions in the [Execution] section below.

#### Advanced configuration options

There are two sources of advanced configuration options:

* Options listed in the file: `${MANTA_INSTALL_DIR}/bin/configManta.py.ini`
    * These parameters are not expected to change frequently. Changing the file
  listed above will re-configure all manta runs for the installation. To change
  parameters for a single run, copy the configManta.py.ini file to another location,
  change the desired parameter values and supply the new file using the configuration
  script's `--config FILE` option.
* Advanced options listed in: `${MANTA_INSTALL_DIR}/bin/configManta.py --allHelp`
    * These options are indented primarily for workflow development and
  debugging, but could be useful for runtime optimization in some specialized
  cases.

### Execution

The configuration step creates a new workflow run script in the requested run directory:

`{ANALYSIS_RUN_DIR}/runWorkflow.py`

This script is used to control parallel execution of Manta via the [pyFlow][2]
task engine. It can be used to parallelize structural variant analysis via one
of two modes:

1. Parallelized across multiple cores on a single node.
2. Parallelized across multiple nodes on an SGE cluster.

A running workflow can be interrupted at any time and resumed where it left
off. If desired, the resumed analysis can use a different running mode or total
core count.

For a full list of execution options, see:

`{ANALYSIS_RUN_DIR}/runWorkflow.py -h`

Example execution on a single node:

```
${ANALYSIS_RUN_DIR}/runWorkflow.py -m local -j 8
```

Example execution on an SGE cluster:

```
${ANALYSIS_RUN_DIR}/runWorkflow.py -m sge -j 36
```

#### Advanced execution options

These options are useful for Manta development and debugging:

* Stderr logging can be disabled with `--quiet` argument. Note this log is
  replicated to `${ANALYSIS_RUN_DIR}/workspace/pyflow.data/logs/pyflow_log.txt`
  so there is no loss of log information.
* The `--rescore` option can be provided to force the workflow to re-execute
  candidates discovery and scoring, but not the initial graph generation steps.

[1]: http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
[2]: http://ctsa.github.io/pyflow/
