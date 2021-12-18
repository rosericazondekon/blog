---
title: Integrated RNA-seq Analysis Application
type: page
---

The Integrated RNA-seq Analysis Application (IRNAA) is a shiny application built for the purpose of facilitating Differential Gene Expression (DGE) analysis. IRNAA integrates command line tools such as `salmon`, `fastQC`, `multiQC` with DGE `R` packages such as `DESeq2`, `edgeR`, and `limma-voom`, and hides all the abstractions and coding hurdles from the end user. 

IRNAA makes it possible for the user to either preprocess fastq files and import read counts into the application from the fastq preprocessed files, or directly import their read counts into the application, or import RNA-seq datasets directly from the The Cancer Genome Atlas (TCGA). 

Finally, IRNAA also integrates other tools such as `gProfiler` and `WebGestalt` for Gene Pathway Analysis. So far, IRNAA preprocesses both **paired-end** and **single-end** fastq files. IRNAA is a work in progress and might, in the future, integrate other tools such as `STAR`, `samR`, etc.

* Release date: April, 2019
* [Source code](https://github.com/rosericazondekon/irnaa)
