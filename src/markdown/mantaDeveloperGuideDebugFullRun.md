
# Manta Developer Guide - Debugging a full Manta run

[Return to Manta Developer Guide](mantaDeveloperGuide.md)

## Summary

This page describes debug/analysis capabilities which are especially useful to full Manta runs -- for instance for scenarios when iterating on a general improvement to the methods. When debugging a single SV, see the related page on [Debugging a single SV in Manta](mantaDeveloperGuideDebugSingleSV.md)

## Rerun SVCandidateGeneration without recreating the SVLocus graph:

This is useful if you're working on a component of candidate generation/scoring which doesn't impact graph creation, and frequently rerunning a test. To use this option provide the '–rescore" option to runWorkflow.py. When this is provided candidate generation and scoring will always be re-run, but the graph will only be created if it doesn't already exist. Example

    ${RUN_DIR}/runWorkflow.py -m sge -j 24 --rescore


## Comparing VCF output between runs:
To assist in evaluating the quality of predictions from a full run compared to a stable benchmark (master, etc), there is a manta utility script to suppress some of the noise expected from a simple 'diff' of two manta vcfs. It takes two vcf files as arguments. These can be gzipped or uncompressed. A usage example is:

    ${MANTA_GIT_CLONE_DIR}/scratch/util/compareMantaVcfs.bash ../m67_test/m67_12_redo_control/results/variants/diploidSV.vcf.gz m63/results/variants/diploidSV.vcf.gz | grep -v MaxDepth | less


## Options to accelerate a small test case:

If not running an analysis on single small genome segment (see [Debugging a single SV in Manta](mantaDeveloperGuideDebugSingleSV.md)) there are various options to make small bam subsegments run a bit faster:

At config time you can reduce/increase the total number of tasks by making each job do more/less :

```
    --scanSizeMb=scanSizeMb
                        Maximum sequence region size (in Mb) scanned by each
                        task during SV locus graph generation. (default: 12)

    --candidateBins=candidateBins
                        Provide the total number of tasks which candidate
                        generation  will be sub-divided into. (default: 256)
```

At run time you can shut down stderr logging, this log is replicated to $runDir/workflow/pyflow.data/logs/pyflow_log.txt so there is no loss of information. To do so, provide 'runWorkflow.py' with the "--quiet" option.
