# Manta Developer Guide - Assembler Test

[Developer Guide Home](README.md)

## Table of Contents

[//]: # (BEGIN automated TOC section, any edits will be overwritten on next source refresh)

* [Introduction](#introduction)
* [Operation](#operation)

[//]: # (END automated TOC section, any edits will be overwritten on next source refresh)


## Introduction

Manta's assembler is internally accessed through the assembler's API during SV analysis. For certain debugging scenarios it is helpful to be able to invoke the assembler on a specified set of input reads. The application `TestAssembler` provides such a capability, reading in read input from one or more bam files and producing a set of contigs in fasta format as output.

## Operation

Example command-line

    ${MANTA_INSTALL_PATH}/libexec/TestAssembler --align-file foo.bam > contigs.fa

Only use this with very small BAMs -- assumes everything is input.

Manta itself has complex selection and orientation logic. In this routine, everything in the bam is selected as input to the assembler. Read orientation is changed only for unmapped reads with mapped read pairs, in which case the unmapped read will be given the opposite strand orientation of its mapped partner.

