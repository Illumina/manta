# Manta User Guide

Version: @MANTA_VERSION@

## Contents
[TOC]

## Introduction

Manta is a general purpose structural variant caller for optimized for short NGS reads.
It is designed to use reads aligned in BAM format as input. Its output is a set of
structural variant predictions scored according to either diploid genotype
probability in a single sample or somatic probability given a tumor and normal sample.
The structural variant scores incorporate evidence from both fragment spanning and
split reads.

## Capabilities

### Detected variant classes

Manta should be able to detect simple indels.

Manta should not be able to detect the following variant types:

* Power to assemble variants to breakend resolution falls to zero as breakend repeat length approaches the read size.
* Power to detect any breakend falls to (nearly) zero as the breakend repeat length approaches the fragment size.
* The method cannot detect non-tandem repeats

Additionally note that manta does not offer high-level classification of SVs beyond finding simple DNA-adjacentcies.

### Input requirements 

## SV Output

## Run configuration and Execution

