## v1.6.0 - 2019-06-25

### Added
- Add configuration option to turn off the evidence signal filter during candidate generation (DRAGEN-1873)
  - This enables very high sensitivity during high depth tumor-only calling 

### Changed
- Change the SV candidate discovery and genotyping phase from a multi-process to a multi-thread design for better CPU utilization (MANTA-1521)
  - As a result, runtime is faster and less variable than before.
  - Runtime improvements vary by workload and server configuration, for typical WGS workloads on a modern server an improvement of ~5-10% may be expected, but improvements of up to 50% have been observed for cases where work previously was poorly distributed across processes.
  - SGE support is removed with this change.
- Update htslib/samtools to 1.9 (MANTA-1483)

## v1.5.1 - 2019-02-15
This is a minor update from v1.5.0.

### Changed
- Make existingAlignmentStats as a fallback option (MANTA-1487)
  - It prioritizes stats estimated from the alignment file and allows static stats file when stats estimation fails due to insufficient high-confidence read pairs.
- Minor updates on Methods writeup (MANTA-1522)

### Fixed
- Fix a bug in selecting large insertion candidates (MANTA-1496)
- Fix a minor bug in QC small SV alignment (MANTA-1498)

## v1.5.0 - 2018-11-12
This is a major update from v1.4.0, featuring improved precision and stability, a new configurable option of overlapping pair reads, and a few bug fixes. VCF representation is improved by introducing a couple of new filters and representing inversions as two breakend records.

### Added
- Add a configurable option to allow overlapping pairs to be used as evidence (MANTA-1398)
  - The option is available in the configure file configureManta.py.ini

### Changed
- Change SV candidate contig aligners to improve precision (MANTA-1396)
  - Change contig aligners such that variant occurrences are more heavily penalized.
- Fix multi-junction nomination (MANTA-1430)
  - Complex events with more than two junctions are no longer nominated as a group
  - Fix the problem of duplicate detection of the same SV candidate
- Add index to ensure uniqueness of evidence bam filenames (MANTA-1431)
  - It solves the potential problem of name conflicts for evidence bams if the input bam files have the same name while located in different directories.
- Change filters for easy interpretation of multi-sample germline variant vcf (MANTA-1343)
  - Add record-level filter 'SampleFT' when no sample passes all sample level filters
  - Add sample-level filter 'HomRef' for homogyzous reference calls
  - No more sample-level filter will be applied at the record level even if it applies to all samples
- Change representation of inversions in the VCF output (MANTA-1385)
  - Intrachromosomal translocations with inverted breakpoints are now reported as two breakend (BND) records.
  - Previously they were reported in the VCF using the inversion (INV) allele type.

### Fixed
- Fix the bug of stats generation with short reference sequences (MANTA-1459/[#143])
- Fix the evidence significance test in the multi-sample calling mode (MANTA-1294)
  - This issue previously caused spurious false negatives during the multi-sample calling mode. The incidence rate of the problem tended to increase with sample count.

## v1.4.0 - 2018-04-25

This is a major bugfix update from v1.3.2, featuring improved precision and vcf representation, in addition to minor user friendly improvements.

### Changed
- Refine SV candidate filter to improve precision (MANTA-1310)
  - Change assessment of assembled SV candidate contigs such that indel occurrences near breakends are more heavily penalized and indel extension penalties are removed.
- Improve the accuracy of SV breakend position for precise calls (MANTA-1178)
  - Expand the reference sequence for the rare cases where a detected breakend is close to the reference end.
- Improve sensitivity for RNA calling (MANTA-1330)
  - Reduce the length of reference sequence to which an assembled SV candidate contig is aligned.
- Add strict checks and improve error message for BED regions of size less than one (STREL-865)
- Improve user guide to clarify meaning of the INV3/INV5 vcf INFO tags (MANTA-1305)
- Remove BND from vcf ALT field to comply with vcf spec (MANTA-1314)

### Fixed
- Fix the off-by-1 bug for vcf HOMSEQ tag in certain variant types (MANTA-1311)
  - Due to position adjustment for certain breakends of inversions and duplications, HOMSEQ was off-by-1 for those variants.
- Fix a bug in breakend homology logic for large variants (MANTA-1313)
  - The breakend insertion was not considered previously when identifying breakend homology sequence.

## v1.3.2 - 2018-03-02

This is a bugfix update from v1.3.1.

### Changed
- Move remote read retrieval for insertions to a configuration file option (MANTA-1351)
  - This feature was previously hard-coded in the workflow, the default behavior (off for cancer workflows, otherwise on) has not changed.
- Turn off the complexity filter for SV breakend graph loci containing more than two nodes (MANTA-1346)
  - High complexity elements of the SV breakend graph associated with centromeres have always been filtered out to control runtime. With this
    change the filtration has been slightly relaxed to ensure that true variants will not be filtered out by transitive association with a high
    complexity graph node.
- Turn on automated task resubmission for all workflow run modes (MANTA-1354)
  - Failed tasks have always been automatically resubmitted in SGE mode, this is now extended to localhost mode as well.
  - This change is intended to work around sporadic I/O issues on network filesystems.

### Fixed
- Standardize germline FORMAT/GQ VCF tag to Integer type (MANTA-1349)
- Update to pyflow v1.1.20 to close infrequent race condition in task resolution (STREL-853)
  - Under the race condition error, a non-failing task would be logged as failing with the message "TASKNAME has stopped without a traceable cause"

## v1.3.1 - 2018-02-19

This is a bugfix update from v1.3.0, notably providing an htslib update to address issues with running manta from alignments in CRAM format.

### Changed
- Change default minimum scored SV size from 51 to 50 (MANTA-1321)
  - This change is intended to better align with GIAB SV size range conventions.

### Fixed
- Update htslib to incorporate CRAM file query fix (MANTA-1336/[#109])
  - This is expected to resolve all known issues with  manta from alignments in CRAM format.
- Fix RNA-seq split read count regression introduced by MANTA-1332 in the v1.3.0 release (MANTA-1342)
- Stop forcing static linking for zlib (STREL-842)
  - This change brings zlib policy in-line with strelka, such that options (i.e. LD\_PRELOAD) are available for improved compression performance.

## v1.3.0 - 2018-02-01

This is a major update from v1.2.2. It features improvements to candidate SV precision and genotyping accuracy, in addition to minor improvements in stability, runtime and error diagnostics.

### Changed
- Change depth estimation read filter to better match the filter used for variant calling (MANTA-1296)
  - Expected depth per chromosome and local depth per locus/variant are now computed after removing filtered, pcr-duplicate, and secondary reads.
- Lower default memory requirements for scatter phase tasks (MANTA-1307)
  - Reduce from 2Gb to 1.5Gb to enable all cores by default on AWS c4.8xlarge/other c\* servers.
  - Added new `--callMemMb` option to override this value for cases of extreme depth/chimera rate, etc
- Update htslib/samtools to 1.6 (MANTA-1331)
  - Updated from older 1.2 version to include improved checks on corrupted data.
  - Stopped vendoring zlib as part of this update, so zlib (w/ headers) is now a build requirement.
- Update the genotyping model (MANTA-1205)
  - Allow minor evidence of reference allele for hom-alt calls to tolerate a small number of noisy reads.
  - Filter out split reads with poor alignments to both alleles.
- Change the eligibility filter of split reads to only consider the most likely allele (MANTA-1332)
  - A split read is allowed to contribute to the split-read evidence if its alignment to the most like allele passes the eligibility filter.

### Fixed
- Provide clear error message when attempting to configure/run with python3 (MANTA-1285)
- Improve error message/docs for alignment records with unknown sequence (SEQ='\*') (MANTA-1295/[#111])
- Improve error message when two alignments use the same QNAME/read-number (MANTA-1293)
  - Message changed to help end-users track down the issue in the alignment file more easily - now includes chromosome name instead of contig index and 1-indexed alignment position.
- Fix CRC error from python gzip lib when generating evidence bam (MANTA-1270)
  - Remove the filter on evidence reads for SV candidates that are not in the file candidateSV.vcf.gz, because the file now contains all SV candidates without removing duplicates.
- Stop automatically clearing python environment variables (MANTA-1316)
  - This should allow python from certain module systems to be used, but may (rarely) cause instability due to conflicting content in a user's PYTHONPATH.
- Improve estimation of chimeric fragment rate (ie. fraction of reads which are split or in anomalous pairs) (MANTA-1261/[#103])
  - This fraction is used to set signal/noise thresholds important for somatic calling.
  - The secondary/supplementary segments of each split read are no longer counted as separate observations.
  - The method now accounts for many reads being classified as both anomalous and split.

## v1.2.2 - 2017-11-10

This is a minor update from v1.2.1.

### Added
- Add a configurable option `graphNodeMaxEdgeCount` (MANTA-1247)
  - During SV candidate generation, if both nodes of an edge have an edge count higher than the configured value, that edge will be skipped for evaluation.

### Changed
- Test for unsupported BAM SEQ format (MANTA-1265)
  - Test the input bam reads for use of the = symbol in the SEQ field, and provide a clear error message if this is found.
- Verify run directory has not already been configured (MANTA-1252/STREL-734/[#102])
- Update minumum boost version to 1.58 (MANTA-1250)
- Update minimum supported linux OS from Centos 5 to 6 (MANTA-1249)

### Fixed
- Fix imprecise SV filtering when CIEND is a subset of CIPOS (MANTA-1146)
- Ensure consistent BND pairs for translocations or RNA fusions are selected during vcf merging (MANTA-1243)
- Fix manta to tolerate deserialization differences in boost above/below v1.58 (MANTA-1262)
  - Impact of issue was an (infrequent) assertion using boost v1.58+: `Assertion 'size() == rhs.size()' failed`
  - Change in boost policy for certain vectors causes vector append in some instances that were previously overwritten.
- Fix target region retrieval (MANTA-1264)
  - The bug was in retrieving overlapping subregions specified by "--callRegion" and "--region" when the argument of "--region" has chromosome name only.

## v1.2.1 - 2017-10-06

This is a minor bugfix release from v1.2.0.

### Added
- Use the BAM mate CIGAR (MC) tag, when present, to improve the accuracy of accessing if a read has extended into adapter sequence (MANTA-1097)
- Add sanity check of specified target regions during configuration (MANTA-1218)
  - A configuration error, instead of a runtime error, will be generated if the chromosome portion of the region does not match a chromosome name from the input reference sequence.

### Changed
- Change candidateSV.vcf to include all candidates used to generate final calls (MANTA-1039)
  - The candidates may include redundanta entries.
  - Redundant candidates are still reduced to a single best call in the final call output.
- Move changelog to markdown format (MANTA-1245)

### Fixed
- Fix an issue in the contig output feature introduced in v1.2.0/MANTA-1207, such that contigs are provided for all SV types (MANTA-1236)

## v1.2.0 - 2017-09-27

This is a major feature update from v1.1.1

### Added
- MANTA-696 Improve the iterative assembler and set it to the default assembler
- MANTA-1207 Add an optional feature to output contig in vcf
- MANTA-1215 Add an option for specification of an existing alignment stats file
- MANTA-1043 Add more read statistics when 'too few read pair observations' exception
- MANTA-886 Add verification of the extension of alignment input file
- MANTA-1118 Add QC check of read length
- MANTA-1044 Enable build on more architectures
- MANTA-1231 Enable clang5.0 build

### Changed
- MANTA-445 RNA: Select contigs for refinement based on alignment score and read count
- MANTA-443 Improve the check of running into adapter for data with small insert size
- MANTA-1177 Improve CRAM reference handling
- MANTA-1041 Improve exception message context to include associated read and SV candidate info
- MANTA-1116 Sort duplicate vcf entries together to allow merging

### Fixed
- MANTA-1211/[#99] Fix issue with empty regions input to EstimateSVLoci
- MANTA-1038 Fix infrequent graph edge lookup failure
- MANTA-1175 Fix a bug of duplicate RG & PG in evidence bams
- MANTA-865 Fix sortVcf.py to handle vcf files without a header

## v1.1.1
- MANTA-691 Add an option to specify call regions in a bed file
- MANTA-977/[#83] Improve SA tag error reporting
- Manta-894 Add INTRON_OPEN penality to off_edge splicing
- Manta-586 Add RNA filters to vcf
- MANTA-962/[#79] Add a user error for missing ini file
- MANTA-554 Disable multiJunction scoring logic in RNA mode
- MANTA-625 Add RNA-mode specific rnaMinCandidateSize option
    - With default set to disable small variant calling
- MANTA-449 Move RNA-mode output to RnaSV.vcf.gz;
    - remove GT and related fields and filters
    - no longer output unpaired breakends

## v1.1.0
- Manta-677 Stats generation
    - Add a sampling buffer to skip regions with a large number of abnormal reads
    - Add a compiler flag to output fragment-size cdf in stats summary
- Sync shared library/build changes with strelka
- Scratch: Merge src formatting changes from strelka
- Docs: Update comments in automated TOC generation; Fix nonstandard markdown links
- MANTA-583 Document multi-junction handling
- MANTA-518 Add sample level check of genotyping for multi-junction events
- MANTA-486/[#57] clarify rna mode configuration/improve errors
- MANTA-436 Classify indels with unknown size insertions as insertions
- MANTA-304 Merge RNAMaster branch
    - RNA: Include some spliced reads in insert-size estimation
    - RNA: Include split-reads as evidence for imprecise translocations
    - RNA: Use proper-pair BAM flag for abnormal pair detection

## v1.0.3
- MANTA-302 fix the bug of duplicate variant ID
- MANTA-279 fix the bug of reporting variants beyond chromosome end position

## v1.0.2
- SPW-316 update denovo scoring script to write valide vcf when the input vcf contains more than three samples
- Remove old/unused reference utility functions

## v1.0.1
- MANTA-299 documentation and method cleanup for handling split read evidence with a couple of bug fixes
- MANTA-296 Test assembler directly from bam input
- [#43] fix freed string pointer use in sample name

## v1.0.0
- [#36] support setting SGE task runtime limit
- MANTA-290 add debug option to generate candidate SV evidence BAMs

## v0.29.7
- Format change of de novo calls, moving DQ from INFO to FORMAT
- [#40] improve robustness to user locale settings
- close infrequent memory leak in hts_streamer
- A number of changes in build and docs

## v0.29.6
- MANTA-287 add the unit test of consistency of supporting reads
- MANTA-287 fix a bug in the small assembler by checking consistency when adding a new supporting read of contig

## v0.29.5
- [#32] Preserve file path softlinks so that sidecar index files can be found
- Fix field name in denovo variant search script

## v0.29.4
- Add a python script to identify de novo calls
- MANTA-285 add option to keep all temp files to support workflow debug
- MANTA-276 mv configure to top-level, mv guides to docs dir, add methods docs
- MANTA-284 improve windows shell support by shortening very long cmdlines

## v0.29.3
- MANTA-280 filter supplementary reads without SA tags and read pairs with unmatched mate information
- STARKA-306 fix rare chunk size boundary defect in RangeMap

## v0.29.2
- MANTA-277 fix invalid genome region requested during insert size estimation
- Update to pyflow 1.1.12 to improve SGE filesystem delay handling
and fix issue between SGE and recent bash shellshock fix

## v0.29.1
- [#22] Add new manta developer guide to source docs
- [#21] improve fragment size estimation for very short fragments
- RNA: Improve fusion detection sensitivity
- MANTA-261 Transfer stable components from Manta windows port
- MANTA-273 fix support for "csi"-style BAM indices
- MANTA-273/[#14] allow bam index filenames in single-extension (Picard) style

## v0.29.0
- MANTA-267/[#12] Support contig names with colons (for HLA contigs in 1kg hg38)
- MANTA-252 Complete support for CRAM input
- MANTA-264 Remove samtools from manta dependencies
- MANTA-252 Change default chrom depth to median estimate from alignments
- MANTA-263 Improve performance/stability for references with
large numbers of small contigs
- MANTA-261 Transfer stable components from Manta windows port

## v0.28.0
- MANTA-259 Support joint analysis of multiple diploid samples
- MANTA-260 Add per-sample filtration to separate QUAL and GQ filters for diploid case
- MANTA-252 Add fast chrom median depth estimator (partally enables CRAM)
- MANTA-258 Add PL values to diploid output

## v0.27.2
- MANTA-257 Fix rare failure condition for graph merge
- MANTA-255 include zlib in build, simplify win64 development
- MANTA-254 Fix handling of off-edge splicing in the RNA Jump Intron Aligner
- MANTA-253 improve alignment corner cases and debugging features

## v0.27.1
- [#6] Fix assertion caused by filtered graph edges on bin boundaries.
- [#5] Improve robustness to filesystem delay (update ot pyflow v1.1.7)

## v0.27.0
- MANTA-188/[#4] fix off-by-one position issues in some precise duplication
and inversion breakends
- MANTA-229 Add initial support for tumor-only analysis

## v0.26.5
- Update pyflow to v1.1.6: fixes multithread bug introduced in v1.1.5.
    - Manta should be isolated from this issue in theory.

## v0.26.4
- Update license to GPLv3
- Update to relicensed pyflow v1.1.5
- MANTA-244 Handle unstranded RNA data
- MANTA-239 Use RNA bam alignments for ref read scoring
- Fix core/memory auto-detect for OSX

## v0.26.3
- Fix OSX build and demo run (req. update to boost 1.56)
- Update travis CI OSX build, static analyzer

## v0.26.2
- MANTA-242 cleanup code portability and documentation:
    - Updated all build/installation and contributor guidelines
    - Updated Manta user guide
    - Added Travis CI configuration file for clang/gcc build and demo run
    - Minor code edits to clean compile on OS X 10.9, CentOS 5,6,7, Ubuntu 12.04 and 14.04
    - Updated demo: new dataset covers COSMIC HCC1954 variants, added test to verify expected output from demo run.

## v0.26.1
- Remove python reflection from run configuration process, fixes rare config issue
- VCF output formatting: FileData corrected to FileDate

## v0.26.0
- MANTA-224 improve short-fragment handling for RNA
- MANTA-235 kmer reference mask to accelerate RNA contig alignments
- MANTA-232 filter large SVs with no read pair support
- MANTA-236 expand conditions for large insertion search to normalize BWA-mem/Isaac performance
- MANTA-231 expand scoring phase split read search around breakends to find soft clipped ref and alt support.
- MANTA-234 remove discovery pair counts from scored output files
- MANTA-222 support bwa-mem '-M' split read format
- MANTA-218 filter spanning candidates without significant signal/noise
- MANTA-155 handle N's and lowqual bases during assembly
- MANTA-219 Build system improvements for auditing, visual studio support
- MANTA-217 Improve runtime for large min candidate size settings, creating a high speed large-event mode.

## v0.25.0
- RNA - parameter adjustments, additional vcf output, orientation fixes
- MANTA-187 treat split reads symmetrically to improve RNA fusion detection
    - detect split reads directly rather than via associated soft-clip sequence
    - exclude split reads from contributing to local assembly evidence

## v0.24.0
- MANTA-213 FFPE runtime optimization:
    - Recognize new BAM format for fragments shorter than read length, prevent these reads from triggering assembly.
    - Improve insert size distribution estimation by including fragments shorter than read length.
    - Add new runtime instrumentation report for candidate generation
    - Add runtime summary to existing graph report
    - Reduce SW edit matrix size with short k-mer match bounds on ref seq
    - Improve graph noise filtration by testing specific region for evidence signal threshold
    - Use adaptive noise rates in hypoth gen step: background read anomaly rate is used to determine if signal is significant. Currently applies to small assembly candidates only.
    - Collapse redundant assembly candidates for small indels to single copy

## v0.23.1
- MANTA-206 Disable remote read search in T/N analysis
- MANTA-195 Improve filtration of short fragments in FFPE samples

## v0.23.0
- MANTA-200 Add filter for overlapping diploid calls which can't be explained as two haplotypes.
- MANTA-199 Fix low-frequency fragment/breakpoint mismatch (primarily
an issue when large numbers of short ref contigs were used)
- Add CRAM support for individual tools - still need a quick chrom depth estimator for full workflow CRAM support
- c++11 update
- MANTA-185 remove transloc calls with neg position (from circular genome)
- [mantadev] #1 Fix SA split read breakpoint position
- MANTA-183 handle paired/single read mixture in input alignments
- MANTA-182 submit config on cmdline

## v0.22.0
- MANTA-181 fix assembler path and coverage issues
- MANTA-177 filter redundant partial insertions
- MANTA-142 improve contig alignment for large events
- MANTA-160 add pseudo-coloring to assembly, and improve multi-pass read/contig association
- MANTA-170 improve pair allele support accuracy
- MANTA-139 add shadow and chimera reads into pair counts

## v0.21.0
- MANTA-167 semi-mapped som pair correction
- MANTA-156 filter out assm poison reads
- MANTA-161 make assembly robust to seed k-mer selection
- MANTA-164 batch retrieve assembly mate reads
- MANTA-158 improve small indel contig alignment specificity
- MANTA-146 use MAPQ0 mate reads in assembly
- MANTA-48 use shadow reads for split scoring
- MANTA-157 improve shadow read filter
- MANTA-153 fix diploid prior
- MANTA-150 fix SV scoring size cutoff
- MANTA-148 check bam records for region errors
- RNA: rna-scoring

## v0.20.1
- Lower default min candidate size to 8

## v0.20.0
- MANTA-136 turn on conservative large insertion calling
- MANTA-126 multi-junction SV scoring
- MANTA-131 improve large somatic sv specificity with expanded
supporting evidence search
- RNA: track candidate orientation

## v0.19.0
- MANTA-127 RG based insert stats (default off)
- MANTA-128 Improved pair orientation estimation and error checks
- RNA: Improve fusion specificity
- MANTA-125 add experimental large insertion calling (default off)
- MANTA-125 add tier2 permissive split score to reduce small somatic FP deletions
- MANTA-125 add tier2 chimera rate to reduce somatic FP calls

##  v0.18.0
- MANTA-125 modify pair weight for small SVs
- MANTA-120 Improve stability of SV scoring as a function of read length

##  v0.17.0
- Filter SA split read segments by MAPQ value
- MANTA-116 better handle BWA-mem SA split-reads for inversions
- MANTA-118 static libstdc++ for gcc 4.5+

##  v0.16.0
- MANTA-117 add somatic quality score
- fix SA tag parsing
- MANTA-27 accept bam/fasta filenames with spaces

##  v0.15.0
- MANTA-108 combine clip/semi-aligned evidence detection, don't detect overlapping reads as assembly evidence
- MANTA-98 make fewer bam scans during scoring
- MANTA-106 add high depth limit to candgen and assembler
- MANTA-75 Better match reads to SV candidates to improve runtime and lower repeat observations (part 2)
- MANTA-105 filter poorly supported candidates to reduce per-edge compute time
- MANTA-103 fix issue in RNA and WES modes introduced by MANTA-99

##  v0.14.0
- MANTA-102 filter calls with high MQ0 fraction
- MANTA-99 add high-depth graph filter to improve FFPE runtime
- MANTA-100 allow for neighboring variants during assembly
- MANTA-83 sort vcfs in bam chrom order
- MANTA-96 Keep matching read pairs after candidate generation read buffer fills
- MANTA-89 Use semi-mapped read pairs to improve germline/somatic classification.
- MANTA-92 Add edge runtime performance log
- MANTA-75 Better match reads to SV candidates to improve runtime and lower repeat observations
- MANTA-85 Increase uniformity of tags in vcf output

##  v0.13.2
- First complete pass at installation and user guide

##  v0.13.1
- MANTA-81 Fix small indel somatic false negatives introduced in MANTA-63
- MANTA-80 Additional workflow options: run subsections of the genome, finer task parallelization control and merge multiple input BAMs per sample.

##  v0.13.0
- MANTA-63 Incorporate read-pair evidence into small SVs/indels
- MANTA-77 Fix assertion for rna-seq test
- MANTA-17 Include semi-aligned reads in discovery and scoring
- MANTA-69 Update score/write filter to account for CIGAR and SA-read candidates, and new uniform candidate scheme for self-edges.
- MANTA-70 Correct filters to allow for small inversion and tandem dup detection
- MANTA-68 SVLEN not set correctly for non-deletions
- MANTA-64 Improve candidate generation for small regions
- MANTA-43 allow manta installation to be relocated
- MANTA-55 compile python code as part of build/install

##  v0.12.1
- MANTA-58 fix issue with breakends near contig boundaries
- MANTA-61 add markdown-based user guide to build
- MANTA-30 initial integration of known variant tracing framework

##  v0.12
- MANTA-20 incorporate split-reads into quality score
- MANTA-42 SV finder mismatches various read pair / sv-candidate combinations
- MANTA-53 Enable --rescore option in runWorkflow.py
- MANTA-40 Don't call splicing-events in RNA-seq as deletions
- MANTA-20 include split read counts for short reads
- MANTA-44 Fix Rhodobacter analysis

##  v0.11
- Adjust all vcf output to pass vcf-validator
- MANTA-20 fix split read breakpoint location

##  v0.10.1
- Fix low-frequency assertion due to unexpected alignment pattern

##  v0.10
- MANTA-20 Limit split read counts to those uniquely supporting each allele, where P(allele|read)>0.999
- MANTA-20 Add likelihood based QUAL,GQ scores to diploid output, adjust thresholds of somatic output to incorporate ref pairs and split reads.
- MANTA-41 Fails when chrom name not in [a-zA-z0-9_-]+
- MANTA-25 Partial support for BWA-MEM SA split reads
- MANTA-36 Segfault on RNA-Seq BAM input
- MANTA-20 Combined reference spanning read and split read evidence per variant
- MANTA-20 Diploid vcf output for non-tumor sample, diploid genotype inference score still todo
- MANTA-39 prevent crash on large CIGAR deletions
- MANTA-20 split read evidence counts for all large spanning SVs

##  v0.9
- MANTA-20 preliminary work on this branch allows assembly skip and control of min indel candidate size and min indel score size
- MANTA-33 reduce SV graph ram requirement to ~1/3 of its previous value, increase all post-merge task memory requests.
- MANTA-17 merged shadow reads into assembly and adjusted assembly parameters. Large (50+ base) insertion sensitivity improves by ~.35-.4 as a result.
- Improvements to vcf output and cmake build.

##  v0.8
- MANTA-28 Add prototype discovery/local-assembly of small events down to 10 bases
- MANTA-24 Better handle very high depth and chimeric noise complexity based
on BWA-mem FFPE examples
- MANTA-26 Extend fragment stats to provide estimate of full fragment size
distribution
- Large event assembly fixes
- MANTA-23 enable use of pre-existing depth and stats files (for sparse bams)

##  v0.7
- Add assembly of large-event breakends and basepair resolution SV reporting
- MANTA-19 Correctly parse large deletion reads from Isaac and incorporate this into discovery

##  v0.6
- Fix sensitivity problems caused by unexpected proper pair bit settings, fix several self-edge issues. Detect intrachrom variants down to ~2kb.

##  v0.5
- Expand POC calls to include intrachromosomal variants down to ~5kb.
- Minor modifications to method based on FFPE testing.

##  v0.4
- POC somatic transloc output

##  v0.3
- POC translation of graph into candidate transloc vcf

##  v0.2
- working proof of concept denoised sv locus graph

##  v0.1
- initial prototype code tag
