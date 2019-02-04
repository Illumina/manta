# Manta Developer Guide - Debugging a full Manta run

[Developer Guide Home](README.md)

## Table of Contents

[//]: # (BEGIN automated TOC section, any edits will be overwritten on next source refresh)

* [Summary](#summary)
* [Rerun SVCandidateGeneration without recreating the SVLocus graph](#rerun-svcandidategeneration-without-recreating-the-svlocus-graph)
* [Comparing VCF output between runs](#comparing-vcf-output-between-runs)
* [Options to accelerate a small test case](#options-to-accelerate-a-small-test-case)

[//]: # (END automated TOC section, any edits will be overwritten on next source refresh)


## Summary

This page describes debug/analysis capabilities which are especially useful to full Manta runs -- for instance for scenarios when iterating on a general improvement to the methods. When debugging a single SV, see the related page on [Debugging a single SV in Manta](debugSingleSV.md)

## Rerun SVCandidateGeneration without recreating the SVLocus graph

This is useful if you're working on a component of candidate generation/scoring which doesn't impact graph creation, and frequently rerunning a test. To use this option provide the '–rescore" option to runWorkflow.py. When this is provided candidate generation and scoring will always be re-run, but the graph will only be created if it doesn't already exist. Example

    ${RUN_DIR}/runWorkflow.py -j 24 --rescore


## Comparing VCF output between runs
To assist in evaluating the quality of predictions from a full run compared to a stable benchmark (master, etc), there is a manta utility script to suppress some of the noise expected from a simple 'diff' of two manta vcfs. It takes two vcf files as arguments. These can be gzipped or uncompressed. A usage example is:

    ${MANTA_GIT_CLONE_DIR}/scratch/util/compareMantaVcfs.bash ../m67_test/m67_12_redo_control/results/variants/diploidSV.vcf.gz m63/results/variants/diploidSV.vcf.gz | grep -v MaxDepth | less


## Options to accelerate a small test case

If not running an analysis on single small genome segment (see [Debugging a single SV in Manta](debugSingleSV.md)) there are various options to make small bam subsegments run a bit faster:

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
