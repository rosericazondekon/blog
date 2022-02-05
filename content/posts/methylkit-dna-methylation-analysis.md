---
title: "DNA Methylation Analysis using methylKit"
date: 2022-02-05T11:47:00+08:00
description: "A DNA Methylation Analysis tutorial Using the methylKit R package."
tags: ["methylkit", "cpg-islands", "dmr", "promoter", "intron", "exon", "coverage", "dna", "methylation", "fastq"]
type: post
weight: 25
---

**Disclaimer:** This tutorial was originally written on June 04, 2019.

## Introduction
In our [DNA Methylation Sequence Read Alignment](/posts/bisulfite-treated-read-align) tutorial, we showed you how to download and process Bisulfite-treated sequence DNA methylation FASTQ files for read alignment on a reference sequence. In this tutorial, we show you how to run DNA methylation analysis using the Bioconductor [`methylKit`](https://bioconductor.org/packages/release/bioc/html/methylKit.html) R package.

We set our working directory to the `tuto` folder created in a [previous tutorial](/posts/bisulfite-treated-read-align).

```r
setwd('./tuto')
```

Now, let's install all the required packages for this tutorial.

```r
# Indicate package repositories to R...
repositories <- c("https://cloud.r-project.org", 
                   "https://bioconductor.org/packages/3.7/bioc",
                   "https://bioconductor.org/packages/3.7/data/annotation", 
                   "https://bioconductor.org/packages/3.7/data/experiment",
                   "https://www.stats.ox.ac.uk/pub/RWin", 
                   "http://www.omegahat.net/R", 
                   "https://R-Forge.R-project.org",
                   "https://www.rforge.net", 
                   "https://cloud.r-project.org", 
                   "http://www.bioconductor.org",
                   "http://www.stats.ox.ac.uk/pub/RWin")

# Package list to download
packages <- c("BiocGenerics","Biobase","S4Vectors","IRanges",
              "GenomicRanges","GenomeInfoDb","AnnotationDbi", 
              "genomation", "fastseg", "methylKit")

# Install and load missing packages
new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]

if(length(new.packages)){
    install.packages(new.packages, repos = repositories)
}

lapply(packages, library, character.only = TRUE)
```


## 1. Obtaining methylation percentage from sorted Bismark alignments

We read in the methylation calls directly from the `SAM` files obtained from the last tutoral. The SAM files must be sorted and be generated from Bismark aligner. For that purpose, we use the `processBismarkAln()` function from the Bioconductor [`methylKit`](https://bioconductor.org/packages/release/bioc/html/methylKit.html) R package as described below:


```r
# Find all sorted SAM files
file_loc <- file.path(getwd(),'alignment_Bismark')

file.list <- as.list(list.files(file_loc, "\\.sort.sam$",
                                full.names = TRUE, recursive = TRUE))

file.list
```

Now, let's read in the methylation data from the `SAM` files (**this has to be run once**).

The following script may take a long time to run!

```r
# Read the methylation data into an object (Has to be run once)
metyl.obj <- processBismarkAln(
                file.list, 
                sample.id = list("normal1", "normal2",
                                 "cancer1", "cancer2"),
                treatment = c(0,0,1,1), 
                assembly = "hg38",
                read.context = "CpG",
                save.folder = file.path(getwd(), "methylation_calls")
)
```

```r
head(metyl.obj,3)
```

When computed, the set of methylation call files can be imported using the `methRead()` function provided by the `methylKit` R package.

```r
# This code may be run if the methylation call files already exist
file_loc2 <- file.path(getwd(),'methylation_calls')
file.list2 <- as.list(
                list.files(file_loc2, "\\CpG.txt$",
                           full.names = TRUE,
                           recursive = TRUE)
)

file.list2
```

```r
methyl.obj2 <- methRead(
                  file.list2, 
                  sample.id = list("cancer1","cancer2",
                                   "normal1","normal2"), 
                  treatment = c(1, 1, 0, 0), 
                  assembly = "hg38"
)
```

```r
head(methyl.obj2,3)
```


## 2. Quality check and basic features of the data

Let's check the basic stats about the methylation data such as coverage and percent methylation.

```r
# Descriptive statistics on the samples
getMethylationStats(metyl.obj[[1]], plot = FALSE, both.strands = FALSE) # normal1
getMethylationStats(metyl.obj[[2]], plot = FALSE, both.strands = FALSE) # normal2
getMethylationStats(metyl.obj[[3]], plot = FALSE, both.strands = FALSE) # cancer1
getMethylationStats(metyl.obj[[4]], plot = FALSE, both.strands = FALSE) # cancer2
```

```r
# Histogram of % CpG methylation
getMethylationStats(metyl.obj[[1]], plot = TRUE, both.strands = FALSE) # normal1
getMethylationStats(metyl.obj[[2]], plot = TRUE, both.strands = FALSE) # normal2
getMethylationStats(metyl.obj[[3]], plot = TRUE, both.strands = FALSE) # cancer1
getMethylationStats(metyl.obj[[4]], plot = TRUE, both.strands = FALSE) # cancer2
```

From the histograms generated by the previous scripts, the numbers on the bars represent the percentage of locations that are contained in each bin. Percent methylation histogram are expected to have two peaks on both ends.

Alternatively, the histogram of the read coverage per base information can be plotted as well. Experiments that are highly suffering from PCR duplication bias are expected to have a secondary peak towards the right hand side of the histogram.

```r
# Histogram of CpG coverage
getCoverageStats(metyl.obj[[1]], plot = TRUE, both.strands = FALSE) # normal1
getCoverageStats(metyl.obj[[2]], plot = TRUE, both.strands = FALSE) # normal2
getCoverageStats(metyl.obj[[3]], plot = TRUE, both.strands = FALSE) # cancer1
getCoverageStats(metyl.obj[[4]], plot = TRUE, both.strands = FALSE) # cancer2
```


## 3. Filtering samples based on read coverage

It is a good practice to filter samples based on coverage, and discard bases that have coverage below 10X bases that have more than 99.9th percentile of coverage in each sample. This can be achieved with the following code:

```r
filtered.metyl.obj <- filterByCoverage(
                        metyl.obj, 
                        lo.count = 10,
                        lo.perc = NULL,
                        hi.count = NULL,
                        hi.perc = 99.9
)
```

```r
combined_data.filtered
```

Let's assess once again the basic stats about the methylation data such as coverage and percent methylation.

```r
# Histogram of % CpG methylation
getMethylationStats(filtered.metyl.obj[[1]], plot = TRUE, both.strands = FALSE) # normal1
getMethylationStats(filtered.metyl.obj[[2]], plot = TRUE, both.strands = FALSE) # normal2
getMethylationStats(filtered.metyl.obj[[3]], plot = TRUE, both.strands = FALSE) # cancer1
getMethylationStats(filtered.metyl.obj[[4]], plot = TRUE, both.strands = FALSE) # cancer2
```

```r
# Histogram of CpG coverage
getCoverageStats(filtered.metyl.obj[[1]], plot = TRUE, both.strands = FALSE) # normal1
getCoverageStats(filtered.metyl.obj[[2]], plot = TRUE, both.strands = FALSE) # normal2
getCoverageStats(filtered.metyl.obj[[3]], plot = TRUE, both.strands = FALSE) # cancer1
getCoverageStats(filtered.metyl.obj[[4]], plot = TRUE, both.strands = FALSE) # cancer2
```


## 4. Sample correlation

To conduct sample correlation, we will need to merge all samples to one object for base-pair locations that are covered in all samples.

```r
# Merging all samples
merged.obj = unite(filtered.metyl.obj, destrand = FALSE)

# Taking a glance at the data...
head(merged.obj,3)
```

The sample correlation is computed using the `getCorrelation()` function available in the `methylKit` package.

```r
# Sample correlation
getCorrelation(merged.obj, plot = TRUE)
```


## 5. Clustering samples

The `clusterSamples()` function in `methylKit` can be used to perform the hierarchical clustering of the samples based on their methylation profiles.

```r
# Clustering samples...
clusterSamples(merged.obj, dist = "correlation", method = "ward", plot = TRUE)
```

Principal Component Analysis (PCA) is another available method to cluster the samples. We perform PCA using the the `PCASamples()` function, then plot the Screeplot.

```r
# Screeplot (PCA analysis) 
PCASamples(merged.obj, screeplot = TRUE)
```

We may plot PC1 (principal component 1) and PC2 (principal component 2) axis and a scatter plot of our samples on those axes which reveals how they cluster.

```r
# Scatterplot (PCA analysis)
PCASamples(merged.obj)
```


## 6. Getting differentially methylated bases

The function `calculateDiffMeth()` is the main function to calculate differential methylation.

```r
# Finding differentially methylated bases or 
# regions (using 8 cores for faster calculations)
methyl.diff <- calculateDiffMeth(merged.obj, num.cores = 8)
```

Following bit selects the bases that have q-value < 0.01 and percent methylation difference larger than 25\%.

```r
# get hyper methylated bases
myDiff25p.hyper <- getMethylDiff(
                      methyl.diff, 
                      difference = 25, 
                      qvalue = 0.01, 
                      type = "hyper"
)

# get hypo methylated bases
myDiff25p.hypo <- getMethylDiff(
                      methyl.diff, 
                      difference = 25,
                      qvalue = 0.01, 
                      type = "hypo"
 )

# get all differentially methylated bases
myDiff25p <- getMethylDiff(
                      methyl.diff, 
                      difference = 25,
                      qvalue = 0.01
)
```

The `methylKit` package can summarize methylation information over tiling windows or over a set of predefined regions (promoters, CpG islands, introns, etc.) rather than doing base-pair resolution analysis.

```r
tiles <- tileMethylCounts(merged.obj, win.size = 1000, step.size = 1000)

head(tiles, 3)
```

## 7. Differential methylation events per chromosome

We may also visualize the distribution of hypo/hyper-methylated bases/regions per chromosome using the `diffMethPerChr()` function.

```r
# Return a list having per chromosome
# differentially methylation events will be returned
diffMethPerChr(methyl.diff,
               plot = FALSE,
               qvalue.cutoff = 0.01,
               meth.cutoff = 25
              )
```

```r
# visualize the distribution of hypo/hyper-methylated
# bases/regions per chromosome
diffMethPerChr(methyl.diff,
               plot = TRUE,
               qvalue.cutoff = 0.01,
               meth.cutoff = 25
              )
```


## 8. Annotating differential methylation events

We may annotate the differentially methylated regions/bases based on gene annotation. We therefore need to read the gene annotation from a bed file and annotate the differentially methylated regions with that information. Similar gene annotation can be created using the Bioconductor [`GenomicFeatures`](https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html) R package.

Let's first download an annotation `.bed` file in the annotation folder created in our first tutorial on DNA methylation. We get ours from [sourceforge.net](https://sourceforge.net/).

```r
# RDownload Annotation bed file
url <- "https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens"
system(paste("wget -P annotation", file.path(url, "hg38_RefSeq.bed.gz")))
```

```r
gene.obj <- readTranscriptFeatures(
              file.path(getwd(),
                        "annotation",
                        "hg38_RefSeq.bed.gz")
)
```

```r
annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj)
```

Similarly, we can read the CpG island annotation and annotate our differentially methylated bases/regions with them.

```r
# read the shores and flanking regions and name the flanks as shores 
# and CpG islands as CpGi
cpg.obj <- readFeatureFlank(file.path(getwd(), "annotation", "hg38_RefSeq.bed.gz"),
                            feature.flank.name = c("CpGi", "shores"))
```

```r
# convert methylDiff object to GRanges and annotate
diffCpGann <- annotateWithFeatureFlank(
                  as(myDiff25p,"GRanges"), 
                  cpg.obj$CpGi, 
                  cpg.obj$shores, 
                  feature.name = "CpGi", 
                  flank.name = "shores"
)
```

```r
diffCpGann
```

## 9. Regional analysis

We now summarize methylation information over a set of defined regions such as promoters or CpG islands with the `regionCounts()` function.

```r
promoters <- regionCounts(metyl.obj,gene.obj$promoters)

head(promoters[[1]])
```

## 10. Getting the distance to TSS and nearest gene name

After getting the annotation of differentially methylated regions, we can get the distance to TSS and nearest gene name using the `getAssociationWithTSS()` function from genomation package.

```r
# Annotate with gene parts
diffAnn <- annotateWithGeneParts(as(myDiff25p, "GRanges"), gene.obj)

# target.row is the row number in myDiff25p
head(getAssociationWithTSS(diffAnn), 3)
```

Getting the percentage/number of differentially methylated regions that overlap with intron/exon/promoters:

```r
genomation::getTargetAnnotationStats(
              diffAnn, 
              percentage = TRUE, 
              precedence = TRUE
)
```

## 11. Working with annotated methylation events

We may also plot the percentage of differentially methylated bases overlapping with exon/intron/promoters:

```r
genomation::plotTargetAnnotation(
              diffAnn, 
              precedence = TRUE, 
              main = "Differential methylation annotation"
)
```

It is also possible to plot the CpG island annotation showing the percentage of differentially methylated bases that are on CpG islands, CpG island shores and other regions.

```r
genomation::plotTargetAnnotation(
              diffCpGann, 
              col = c("green", "gray", "white"), 
              main = "differential methylation annotation"
)
```


Find the Original source code for this tutorial [here](https://github.com/rosericazondekon/bioinformatics-tuto/blob/main/DNA%20methylation%20analysis/dna_methylation_tuto.ipynb).