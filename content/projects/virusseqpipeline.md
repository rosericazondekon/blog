---
title: VirusSeq Pipeline
type: page
---

The VirusSeq software is a next generation sequencing software of human cancer tissue. Its goal is to identify viruses and their integration sites in human cancer cells. VirusSeq was originally developed by [Chen and collaborators](https://academic.oup.com/bioinformatics/article/29/2/266/202055 "VirusSeq: software to identify viruses and their integration sites using next-generation sequencing of human cancer tissue") and reported in a 2012 paper published in Bioinformatics. 

Although the software is open source and freely available, its installation is not at the ease of the novice user. VirusSeq has not been designed to effectively and efficiently process multiple samples. In its present shape, it requires the user to write a shell script to every single sample to process. 

Furthermore, the user has to engage in the tedious and error prone tasks of path modification. All these limitations of the VirusSeq software may tremendously limit its use in the scientific community. 

In this project, we aim at addressing the aforementioned limitations by designing an easy to install, easy to implement pipeline for the VirusSeq software. Our approach does not require the user to have knowledge of any shell command.

* Release date: February, 2018
* [Source code](https://github.com/rosericazondekon/virusSeqPipeline)
