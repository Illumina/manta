# Manta User Guide

Version: @MANTA_VERSION@

## Contents
[TOC]

## Introduction

Manta is a general purpose structural variant caller for short sequencing reads.

## Capabilitiies

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

