# genomecompR

This R package is provided to describe the method section of an accompanying article. It does not aim to be extended in terms of features. Please respect the packages version described below for full functionality.

## Description

*genomecompR* creates genomic compartments using ChIP-seq and ATAC-seq data as follows:

1) Active promoter: K27ac in TSS-/+1Kb.
2) Transcription initiation: TSS-1/+1kb overlapping with Ser5P peaks.
3) Transcription elongation: TSS+1kb to TES overlapping Ser2P.
4) Transcription termination: TES+50bp
5) Bivalent promoters: H3K4me3/H3K27me3 overlap TSS-/+1Kb. 
6) Active enhancer: H3K27ac/H3K4me1/ATAC-seq. They should not overlap with the combination of UCSC refGene, NCBI RefSeq, and GENCODE VM25.
7) Poised enhancer: H3K27me3/H3K4me1/PRC2. They should not overlap with the combination of UCSC refGene, NCBI RefSeq, and GENCODE VM25.
8) Polycomb domain: Suz12 and RING1B overlap.
9) Heterochromatin: H3K9me3
10) SINE: Overlap with SINE annotations of repeat maskers.
11) LINE: Overlap with LINE annotations of repeat maskers.
12) LTR: Overlap with LTR annotations of repeat maskers.

*genomecompR* also provides several plotting methods: upsetDiagram, boxplotGlcnacLevels, outputGlcPeaksCoordPerCompartment, retrieveGlcPeakVal, complexUpsetDiagram, extractCompCoordWithPeak, violinplotExpression

## Installation

This package is made to run on R-4.2.0. See installation instruction [here](https://cran.r-project.org/) or use:

```
## with mamba
mamba create -n genomecompr
mamba activate genomecompr
mamba install r-base=4.2.0
```

Here is the installation code to run from R:

```
## not tested yet
install.packages("MASS")
install.packages("remotes")
library("remotes")
install_version("ggplot2", version="3.4.1")

ERROR: dependencies ‘MASS’, ‘mgcv’ are not available for package ‘ggplot2’


install("ggplot2", version="3.4.1")
install("tibble", version="3.1.6")
install("ggupset", version="0.3.0")
install("ComplexUpset", version="1.3.3")

install("S4Vectors", version="0.36.2")
install("GenomeInfoDb", version="1.34.9")
install("rtracklayer", version="1.58.0")
install("biomaRt", version="2.54.1")
install("IRanges", version="2.32.0")
install("GenomicRanges", version="1.50.2")
install("chipenrich", version="2.22.0")



  package,
  version = NULL,
  dependencies = NA,
  upgrade = c("default", "ask", "always", "never"),
  force = FALSE,
  quiet = FALSE,
  build = FALSE,
  build_opts = c("--no-resave-data", "--no-manual", "--no-build-vignettes"),
  build_manual = FALSE,
  build_vignettes = FALSE,
  repos = getOption("repos"),
  type = "source",

  install.packages("BiocManager", repos = "https://cloud.r-project.org")
library(BiocManager)

  ...
)


```

## Usage

See the following scripts: XXX, XXXX

## Manual

See [manual.pdf](https://github.com/descostesn/genomecompR/blob/main/manual.pdf)
