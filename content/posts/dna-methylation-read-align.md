---
title: "DNA Methylation Read Sequence Alignment"
date: 2022-01-08T13:08:00+08:00
description: "A tutorial to Bisulfite-seq DNA Methylation read alignment with Bismark Aligner."
tags: ["dna", "methylation", "cpg-islands", "bisulfite-seq", "bam", "sam", "fastq", "alignment"]
type: post
weight: 15
---

**Disclaimer:** This tutorial was originally written on May 30, 2019.

## Introduction
In this tutorial, we show you how to download raw Bisulfite-seq DNA methylation sequence data from the European instance of the SRA, which can be accessed via https://www.ebi.ac.uk/ena. On the European Nucleotide Archive (ENA) website, the sequencing reads are directly available in FASTQ or SRA formats, which will be explained below.

For this tutorial, we need [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [`multiQC`](https://multiqc.info/), the [`SRA toolkit`](https://hpc.nih.gov/apps/sratoolkit.html), a powerful suite of tools designed to interact with SAM and BAM files called [`samtools`](https://sourceforge.net/projects/samtools/files/), and the [`Bismark`](https://github.com/FelixKrueger/Bismark) aligner to align the Bisulfite-seq reads to the reference genome.

All the above mentioned tools need to be installed and referenced in the environment variable `PATH`.

We set our working directory to the `tuto` folder created in our first tutorial.

```r
setwd('./tuto')
```

Let's first check if these requirements are met:

```shell
fastqc --version
```

```shell
multiqc --version
```

```shell
fastq-dump --version
```

```shell
samtools --version
```

```shell
bismark --version
```


If at least one of the above commands produces an error, please, check your installation of the tool and try again.

Now let's create a working directory for our bisulfite-seq DNA methylation project (let's call it `tuto`).

```shell
mkdir -p tuto && cd tuto
```

## 1. Data Download
To download a set of SRA files:

1. Go to https://www.ebi.ac.uk/ena.
2. Search for the accession number of the project, e.g., SRP041828 (should be indicated in the published paper).
3. There are several ways to start the download, here we show you how to do it through the command line interface on GNU/Linux.

    * copy the link’s address of the "SRA files" column (right mouse click), go to the command line, move to the target directory, type: `wget < link copied from the ENA website >`
    * If there are many samples as it is the case for the project referenced here (accession number: SRP041828), you can download the summary of the sample information from ENA by right-clicking on "TEXT" or "TSV" and copying the link location.

Now let's download the file from the link copied earlier.

```shell
wget -O all_samples.txt "https://www.ebi.ac.uk/ena/portal/api/filereport?\
accession=PRJNA274258&result=read_run&fields=study_accession,sample_accession,\
experiment_accession,run_accession,tax_id,scientific_name,fastq_ftp,\
submitted_ftp,sra_ftp&format=tsv&download=true&limit=0"
```

You may try to open the `all_samples.txt` file with LibreOffice or Excel to view it. For this project, we are only interested in the paired-end first 4 Bisulfite-seq samples (2 normal cells samples vs 2 breast cancer cells samples). Since the first line in `all_samples.txt` contains the header, we will generate another file containing only the first 4 lines of  `all_samples.txt` with the following command:

```shell
sed '1d' all_samples.txt > all_samples2.txt
head -4 all_samples2.txt > samples.txt
rm all_samples2.txt
```

Now, let's create a new folder for our SRA files.

```shell
mkdir -p sra_files
```

Since the decommission of the SRA Fuse/FTP site on December 1, 2019, the `prefetch` utility included in the `SRA toolkit` is recommended.

<!-- According to https://www.ncbi.nlm.nih.gov/books/NBK158899/, the FTP root to download files from NCBI is ftp://ftp-trace.ncbi.nih.gov/ and the remainder path follow the specific pattern `/sra/sra-instant/reads/ByRun/sra/{SRR|ERR|DRR}/<first 6 characters of accession>/<accession>/<accession>.sra`.-->

Notice that the accession number for the SRA files are located in the fourth column `run_accession` in all_samples.txt. We proceed to the download of the SRA files of the samples listed in samples.txt with the following code: (**Attention: The download may take a long time!**)

```shell
# The command below may take too long to download.
cut -f4 samples.txt | xargs -i sh -c \
        'run_accession={}; \
         prefetch $run_accession --output-directory sra_files'
```

## 2. Converting SRA files to FASTQ files
Now that the download is complete, let's convert the SRA files into FASTQ files with the following command: (**Attention: This may take a long time!**)

```shell
cut -f4 samples.txt | xargs -i sh -c \
    'run_accession={}; \
     fastq-dump --outdir fastq/${run_accession} \
                --gzip \
                --skip-technical \
                --split-3 sra_files/${run_accession}/${run_accession}.sra'
```


## 3. Quality Control of the FASTQ files
Up to this point, we have all our Bisulfite-seq DNA methylation FASTQ files ready for Quality Control (QC) check. This is done with the `fastqc` tools developed by the Babraham Institute. Run the following command to perform QC check for all the samples: (**This may take some time!**)

```shell
cut -f4 samples.txt | xargs -i sh -c 
    'run_accession={}; \
     mkdir -p fastqc_reports/${run_accession}; \
     fastqc fastq/${run_accession}/*fastq.gz -o fastqc_reports/${run_accession}'
```

Next, let's summarize the QC reports (for all the samples) into one unique report using `multiqc`:

```shell
multiqc fastqc_reports --dirs -o multiQC_report/
```

Let's examine the summary `multiqc` report either by double-clicking on `multiQC_report/multiqc_report.html` or by executing the following code:

```shell
xdg-open multiQC_report/multiqc_report.html
```

## 4. Read Alignment

The assignment of sequencing reads to the most likely locus of origin is called read alignment or mapping and it is a crucial step in most types of high-throughput sequencing experiments.

The general challenge of short read alignment is to map millions of reads accurately and in a reasonable time, despite the presence of sequencing errors, genomic variation and repetitive elements. The different alignment programs employ various strategies that are meant to speed up the process (e.g., by indexing the reference genome) and find a balance between mapping fidelity and error tolerance.


### 4.1. Reference genomes and annotation

Genome sequences and annotation are often generated by consortia such as (mod)ENCODE, The Mouse Genome Project, The Berkeley Drosophila Genome Project, and many more. The results of these efforts can either be downloaded from individual websites set up by the respective consortia or from more comprehensive data bases such as the one hosted by the University of California, Santa Cruz ([UCSC](https://genome.ucsc.edu/)) or the European genome resource ([Ensembl](http://www.ensembl.org/)).

Reference sequences are usually stored in plain text FASTA files that can either be compressed with the generic gzip command. 

The reference sequences file can be obtained from [NCBI](https://www.ncbi.nlm.nih.gov/genome/51), [ENSEMBL](https://www.ensembl.org/info/data/ftp/index.html) or [UCSC Genome Browser](http://hgdownload.soe.ucsc.edu/downloads.html#human).

For this DNA methylation (Bisulfite-seq) tutorial, we will align the reads against the genome (DNA) reference sequences. We obtain both the genome refernce sequences and our gene annotation files from [UCSC](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/). This is very important as we intend to perform all downstream DNA methylation analysis using the [`methylKit`](https://bioconductor.org/packages/release/bioc/html/methylKit.html) Bioconductor R package which works nominally with UCSC genome references.

```shell
# Download the latest human genome
wget -P reference http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
```

### 4.2. Aligning reads using `Bismark` aligner

#### 4.2.1. Generate genome index
**This step has to be performed only once per genome type (and alignment program)**.

```shell
bismark_genome_preparation --verbose ./reference
```

#### 4.2.2. Alignment
This step has to be performed for each individual FASTQ file.

**This step may take a long time! (may take several days to complete).**

```shell
# execute Bismark aligner
cut -f4 samples.txt | xargs -i bash -c \
  'run_accession={}; 
  mkdir -p alignment_Bismark/${run_accession}; \
  bismark --parallel 8 
          --gzip \
          --fastq \
          --output_dir alignment_Bismark/${run_accession} \
          --genome ./reference -1 fastq/${run_accession}/${run_accession}_1.fastq.gz \
                               -2 fastq/${run_accession}/${run_accession}_2.fastq.gz'
```

#### 4.2.3. Sorting BAM files and converting to SAM files

We sort sort the `BAM` files using the `samtools sort` command:

```shell
# Sorting the bam files and converting to 
for i in alignment_Bismark/*/*; do
    if [ "${i}" != "${i%pe.bam}" ];then
        samtools sort -l 0 \
                      -T $(dirname ${i})/$(basename ${i} \
                                  _1_bismark_bt2_pe.bam)_temp \
                      -O sam -@ 8 \
                      -o $(dirname ${i})/$(basename ${i} .bam).sort.sam ${i}
    fi
done
```

Either [SeqMonk](https://www.bioinformatics.babraham.ac.uk/projects/download.html#seqmonk) or the the [Integrative Genomics Viewer (IGV)](https://software.broadinstitute.org/software/igv/) can be used to visualize the resulting sorted `SAM` files.

We will later use the `methylKit` Bioconductor R package to import the methylation data into `R` from the sorted `SAM` files.


## 5. Methylation extraction using Bismark methylation extractor

With the `bismark_methylation_extractor` command, we extract the methylation call for every single Cytosine analyzed. This process takes as input the resulting `BAM` file from `Bismark` aligner. The `bismark_methylation_extractor` command writes the position of every single Cytosine to a new output file, depending on its context (CpG, CHG or CHH), whereby methylated Cytosines are labelled as forward reads (+), non-methylated Cytosines as reverse reads (-).

`SeqMonk`, a genome viewer, can be used to visualize the output files.

We store the output of the `Bismark` methylation extractor in the `methylation_data` folder.

```shell
# Extract methylation data
cut -f4 samples.txt | xargs -i bash -c \
        'run_accession={}; \
         mkdir -p bismark_methCalls/${run_accession}; \
         bismark_methylation_extractor \
              --parallel 8 \
              --gzip \
              --bedGraph \
              --buffer_size 40G \
              --merge_non_CpG \
              --comprehensive \
              --output bismark_methCalls/${run_accession} alignment_Bismark/${run_accession}/*_pe.bam'
```

In another tutorial, we will analyze DNA methylation data from the generated sorted SAM files from this tutorial using the `MethylKit` Bioconductor package in R.

Find the Original source code for this tutorial [here](https://github.com/rosericazondekon/bioinformatics-tuto/blob/main/DNA%20methylation%20analysis/bisulfite_seq_tuto.ipynb).