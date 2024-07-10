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


## Manual

See [manual.pdf](https://github.com/descostesn/genomecompR/blob/main/manual.pdf)
