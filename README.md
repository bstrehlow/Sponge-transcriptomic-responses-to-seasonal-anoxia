# Sponge-transcriptomic-responses-to-seasonal-anoxia

These scripts and data are part of the following study: 
‘Transcriptomic responses of sponge holobionts to in situ, seasonal anoxia and hypoxia.’
Two sponge species, Eurypon sp. 2 (Es) and Hymeraphia stellifera (Hs), survive seasonal anoxia for months at a time. Each sponge species possessed a unique microbiome, but the microbiomes of each species were dominated by a species-specific Thaumarchaeon (Thaum) and a Gammaproteobacterium (Gamma). For each sponge species, a reference transcriptome was created, and metagenome assembled genomes (MAGs) were generated for their symbionts. Full methods are elaborated in the paper. 
Additional sequencing data and assemblies are available at NCBI under BioProject number PRJNA893197.

## Assemblies 
Scripts for the creation of transcriptomes for two sponge species and MAGs for their symbionts are in the folder: Assembly_annotations_and_scripts. Annotations in the form of indexed tab files for gene names, GO terms, KEGG and KOG are also in this folder using formats and scripts from https://github.com/z0on/annotatingTranscriptomes. All protein coding genes from the sponge mitochondria are in the folder: Mitochondrial_assemblies.

## Differentially expressed genes 
The R script for determining differentially expressed genes can be found in: Differential_Expression_Analyses. Differential expression analyses with performed with Desq2 (Love et al., 2014). Analyses of GO and KOG enrichments were done using scripts in said folder modified from https://github.com/z0on/GO_MWU and https://cran.r-project.org/web/packages/KOGMWU/index.html, respectively.

Lists of all differentially expressed genes for each sponge species and their microbial symbionts are present in *_report.pdf files. These files also include boxplots of expression for each differentially expressed gene as well as logFC, p values, and adjusted p values. Report files were made using the ReportingTools package in R (Huntley et al., 2013).

## Refrences
Huntley MA, Larson JL, Chaivorapol C, Becker G, Lawrence M, Hackney JA, Kaminker JS (2013). “ReportingTools: an automated result processing and presentation toolkit for high throughput genomic analyses.” Bioinformatics. doi: 10.1093/bioinformatics/btt551.

Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550. doi: 10.1186/s13059-014-0550-8.
