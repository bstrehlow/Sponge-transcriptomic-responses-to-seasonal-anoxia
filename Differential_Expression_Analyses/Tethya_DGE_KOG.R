Sys.setenv(LANG = "en")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.12")

#have to run R as an 

#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#Install packages 

#BiocManager::install("BiocStyle")
library(BiocStyle)
#install.packages('rmarkdown')
library(rmarkdown)
#BiocManager::install('geneplotter')
library(geneplotter)
#^might not have worked due to permissions???
#install.packages('ggplot2')
library(ggplot2)
#install.packages('plyr')
library(plyr)
#install.packages('LSD')
library(LSD)
#BiocManager::install("DESeq2")
#^might not have worked due to permissions??? but the package loads..
library(DESeq2)
#install.packages('gplots')
library(gplots)
#install.packages('RColorBrewer')
library(RColorBrewer)
#install.packages('stringr')
library(stringr)
#BiocManager::install('topGO')
#^might not have worked due to permissions???
library(topGO)
#BiocManager::install('genefilter')
#^might not have worked due to permissions???
library(genefilter)
#BiocManager::install('biomaRt')
#^might not have worked due to permissions???
library(biomaRt)
library(dplyr)
#BiocManager::install('EDASeq')
#^might not have worked due to permissions???
library(EDASeq)
#BiocManager::install('fdrtool')
#^permission issue again??
library(fdrtool)
#BiocManager::install('org.Mm.eg.db')
#^permission issue again??
#install.packages('ggpubr')
library(ggpubr)
#do we really need this one_
#library(org.Mm.eg.db)
#install.packages('pheatmap')
library(pheatmap)
#BiocManager::install('apeglm')
library(apeglm)
#BiocManager::install('ashr')
#^permission issue again??
library(ashr)
#BiocManager::install("ReportingTools")
#^permission issue again??
library(ReportingTools)
#install.packages('stringi')
library(stringi)
library(KOGMWU)
#install.packages('tidyverse')
library(tidyverse)

setwd("C:/Users/strehlow/OneDrive - Syddansk Universitet/SDU/Transcriptomics/DGE_lough_hyne/Final_files")

####################
####################################
#generating KOG enrichment files for T. wilhelma from Mills and Francis et al 2018 eLife
#Mills, D. B., Francis, W. R., Vargas, S., Larsen, M., Elemans, C. P. H., Canfield, D. E., & Wörheide, G. (2018). The last common ancestor of animals lacked the HIF pathway and respired in low-oxygen environments. ELife, 7, 1-17. https://doi.org/10.7554/eLife.31176
####################################

#two separate experiments
#short term exporsure to anoxia (shock)
#long term exposure to anoxia (long)

#start with short term

counts_tHost=read.table("twi_combined_curatedgenes_mito.fasta.despliced_Shock.csv",header=TRUE,row.names=1) #Reading in the table of counts per isogroup by sample - you must open R in the same folder or change the working directory

counts_t=counts_tHost[,order(names(counts_tHost))]
names(counts_t)

##conditions (metadata)
conditions=read.table("Twi_LowO2_Shock.info", header = T, sep=',')

conds_t <- conditions[ order(row.names(conditions)), ]

####load data into deseq table format
DESeq2table_t <- DESeqDataSetFromMatrix(countData = counts_t,
                                         colData = conds_t,
                                         design= ~ condition)

#make control the 'reference' level 
DESeq2table_t$condition <- relevel(DESeq2table_t$condition, "Control")

#runs differential expression anaylsis quickly
dds_t <- DESeq(DESeq2table_t)

resultsNames(dds_t) # lists the coefficients
res_t_AvN <- results(dds_t, name="condition_Shock_vs_Control")
# or to shrink log fold changes association with condition:

#significantly expressed genes based on full results
table(res_t_AvN$padj < 0.1)

#FALSE  TRUE 
#16817 10220 

#threshold at 0.01 matches Warren's data as well :).

#prep for PCA/heatmaps
#transformation 
vsd <- vst(dds_t, blind=FALSE)

df <- as.data.frame(colData(dds_t)[,c("Oxy_cat")])

#sample to sample distances
sampleDists <- dist(t(assay(vsd)))

#heatmap of distance matrices between samples 
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(vsd), vsd$Oxy_cat, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#not all anoxic samples group together. One anoxic sample and one hypoxic sample are very close to eachother 

#pca plot 
plotPCA(vsd, intgroup=c('condition'))
#looks like figure from paper :)

des2Report <- HTMLReport(shortName = 'Tethya_shock', title = 'Tetha_shock',reportDirectory = "./reports")
publish(dds_t,des2Report, pvalueCutoff=0.1, factor = colData(dds_t)$condition, reportDir="./reports")
finish(des2Report)

#Tethya control vs shock
#log fold change less than 0 -> downregulated in anoxic shock compared to control
#log fold change greater than 0 - upregulated in anoxia

#prepared csv for GO analysis using go_mwu https://github.com/z0on/GO_MWU
#df2 has all differentially expressed genes
#make log p, signed based on up or down regulation 
#need to use dataframe with only isogroup designation (dds_t2)

####################### KogMWU setup using log fold change values

#make files with log fold changes 
#anoxic shock vs control
df1 <-data.frame(res_t_AvN)
df2 <- df1 

table(df2$padj < 0.1)
#avh LFC >0 -> upregulated in anoxia
df2$lfc <- df2$log2FoldChange

#make csv with just gene names and lfc
kog <- data.frame(rownames(df2))
kog$lfc <-df2$lfc 
names(kog)[1] <- 'gene_id'

#write as csv
write.csv(kog, file='t_shock_kog_lfc_final.csv', row.names = FALSE)

#KOG mwu tethya

t_kog2gene=read.table("tethya_all_gene2kog_unique.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)

#lfc for comparison
t_shock.lfc <- read.csv('t_shock_kog_lfc_final.csv', header =TRUE)

t_shock.lth=kog.mwu(t_shock.lfc,t_kog2gene) 
t_shock.lth 


#long term experiment

counts_t2Host=read.table("twi_combined_curatedgenes_mito.fasta.despliced_LongTreat.csv",header=TRUE,row.names=1) #Reading in the table of counts per isogroup by sample - you must open R in the same folder or change the working directory

counts_t2=counts_t2Host[,order(names(counts_t2Host))]
names(counts_t2)

##conditions (metadata)
conditions=read.table("Twi_LowO2_Long.info", header = T, sep=',')

conds_t2 <- conditions[ order(row.names(conditions)), ]

####load data into deseq table format
DESeq2table_t2 <- DESeqDataSetFromMatrix(countData = counts_t2,
                                        colData = conds_t2,
                                        design= ~ condition)

#make control the 'reference' level 
DESeq2table_t2$condition <- relevel(DESeq2table_t2$condition, "Control")

#runs differential expression anaylsis quickly
dds_t2 <- DESeq(DESeq2table_t2)

resultsNames(dds_t2) # lists the coefficients
res_t2_AvN <- results(dds_t2, name="condition_LongTreat_vs_Control")
# or to shrink log fold changes association with condition:

#significantly expressed genes based on full results
table(res_t2_AvN$padj < 0.1)

#FALSE  TRUE 
#24439   400 

#threshold at 0.01 matches Warren's data as well :).

#prep for PCA/heatmaps
#transformation 
vsd <- vst(dds_t2, blind=FALSE)

df <- as.data.frame(colData(dds_t2)[,c("condition")])

#sample to sample distances
sampleDists <- dist(t(assay(vsd)))

#heatmap of distance matrices between samples 
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(vsd), vsd$Oxy_cat, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#not all anoxic samples group together. One anoxic sample and one hypoxic sample are very close to eachother 

#pca plot 
plotPCA(vsd, intgroup=c('condition'))
#looks like figure from paper :)

des2Report <- HTMLReport(shortName = 'Tethya_long', title = 'Tetha_long',reportDirectory = "./reports")
publish(dds_t2,des2Report, pvalueCutoff=0.1, factor = colData(dds_t2)$condition, reportDir="./reports")
finish(des2Report)

#Tethya control vs shock
#log fold change less than 0 -> downregulated in anoxic shock compared to control
#log fold change greater than 0 - upregulated in anoxia

#prepared csv for KOG analysis using KOG https://github.com/z0on/GO_MWU
df1 <-data.frame(res_t2_AvN)
df2 <- df1 

table(df2$padj < 0.1)
#avh LFC >0 -> upregulated in anoxia
df2$lfc <- df2$log2FoldChange

#make csv with just gene names and lfc
kog <- data.frame(rownames(df2))
kog$lfc <-df2$lfc 
names(kog)[1] <- 'gene_id'

#write as csv
write.csv(kog, file='t_long_kog_lfc_final.csv', row.names = FALSE)

#KOG mwu tethya

t_kog2gene=read.table("tethya_all_gene2kog_unique.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)

#lfc for comparison
t_long.lfc <- read.csv('t_long_kog_lfc_final.csv', header =TRUE)

t_long.lth=kog.mwu(t_long.lfc,t_kog2gene) 
t_long.lth 

#tethya KOG heatmap 

t_shock.lth
# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("Long"=t_long.lth,"Shock"=t_shock.lth))
# Making a heatmap with hierarchical clustering trees: 
pheatmap(as.matrix(ktable)) 

# exploring correlations between datasets
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor)
# p-values of these correlations in the upper panel:
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor.pval)

#save KOG tables 
#no significantly enriched KOGs in long exposure 
#Translation, ribosomal structure and biogenesis significantly enriched in short term shock

#Is the mitochondrial transcript (only one in this dataset) differentially expressed
df1 <-data.frame(res_t_AvN)

df1[c('Twilhelma_mtDNA_genome_NODE_33_length_19666_cov_12'),]

df1 <-data.frame(res_t2_AvN)

df1[c('Twilhelma_mtDNA_genome_NODE_33_length_19666_cov_12'),]

#no significant differences in mito gene expression

