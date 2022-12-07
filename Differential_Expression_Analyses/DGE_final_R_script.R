Sys.setenv(LANG = "en")

#uncomment installation scripts if needed for bioc packages

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
#install.packages('ggplot2')
library(ggplot2)
#install.packages('plyr')
library(plyr)
#install.packages('LSD')
library(LSD)
#BiocManager::install("DESeq2")
library(DESeq2)
#install.packages('gplots')
library(gplots)
#install.packages('RColorBrewer')
library(RColorBrewer)
#install.packages('stringr')
library(stringr)
#BiocManager::install('topGO')
library(topGO)
#BiocManager::install('genefilter')
library(genefilter)
#BiocManager::install('biomaRt')
library(biomaRt)
library(dplyr)
#BiocManager::install('EDASeq')
library(EDASeq)
#BiocManager::install('fdrtool')
library(fdrtool)
#BiocManager::install('org.Mm.eg.db')
#install.packages('ggpubr')
library(ggpubr)
#do we really need this one_
#install.packages('pheatmap')
library(pheatmap)
#BiocManager::install('apeglm')
library(apeglm)
#BiocManager::install('ashr')
library(ashr)
#BiocManager::install("ReportingTools")
library(ReportingTools)
#install.packages('stringi')
library(stringi)
library(KOGMWU)
#install.packages('tidyverse')
library(tidyverse)

setwd("C:/Users/strehlow/OneDrive - Syddansk Universitet/SDU/Transcriptomics/DGE_lough_hyne/Final_files")

#start with eurypon.
####################
####################################
#Eurypon sp2 (es) sponge (host)
####################################

##will need multiple comparisons too 

counts_esHost=read.table("allcounts_es_sponge_final.txt",header=TRUE,row.names=1) #Reading in the table of counts per isogroup by sample - you must open R in the same folder or change the working directory

counts_es=counts_esHost[,order(names(counts_esHost))]
names(counts_es)

##conditions (metadata)
conditions=read.table("Conditions_es.txt", header = T)

conds_es <- conditions[ order(row.names(conditions)), ]

####load data into deseq table format
DESeq2table_es <- DESeqDataSetFromMatrix(countData = counts_es,
                                         colData = conds_es,
                                         design= ~ Oxy_cat)

#make normoxic the 'reference' level 
DESeq2table_es$Oxy_cat <- relevel(DESeq2table_es$Oxy_cat, "Normoxic")

#runs differential expression anaylsis quickly
dds_es <- DESeq(DESeq2table_es)

resultsNames(dds_es) # lists the coefficients
res_es_AvN <- results(dds_es, name="Oxy_cat_Anoxic_vs_Normoxic")
# or to shrink log fold changes association with condition:

#shrunken results beased on LFC
res_es_AvN_apeglm <- lfcShrink(dds_es, coef="Oxy_cat_Anoxic_vs_Normoxic", type="apeglm")

#significantly expressed genes based on full results
table(res_es_AvN$padj < 0.1)

#FALSE  TRUE 
#27499  1571 
#1571 differentially expressed genes between normoxic and anoxic samples :) es

#significantly expressed genes based on LFC results
table(res_es_AvN_apeglm$padj < 0.1)

##1571 differentially expressed genes between normoxic and anoxic samples = same as above :)

#standard method - hypoxic vs normoxic
res_es_HvN <- results(dds_es, name="Oxy_cat_Hypoxic_vs_Normoxic")

table(res_es_HvN$padj < 0.1)
#FALSE  TRUE 
#28141   929 
#929 differentially expressed genes between hypoxic and normoxic

# or to shrink log fold changes association with condition:
res_es_HvN_apeglm <- lfcShrink(dds_es, coef="Oxy_cat_Hypoxic_vs_Normoxic", type="apeglm")

table(res_es_HvN_apeglm$padj < 0.1)
#same as above

length(res_es_AvN[,1])

#44249 genes in total for es

#comparison (contrast), between hypoxic and anoxic samples
res_es_HvA <- results(dds_es, contrast = list('Oxy_cat_Hypoxic_vs_Normoxic', 'Oxy_cat_Anoxic_vs_Normoxic'))

table(res_es_HvA$padj < 0.1)
#FALSE  TRUE 
#28791  1093

#1093 genes differentially expressed between hypoxic and anoxic

#shrunken results 
#apeglm does not have contrast ability
#need to try 'normal shrinkage, or ashr

res_es_HvA_ashr <- lfcShrink(dds_es, contrast = list("Oxy_cat_Anoxic_vs_Normoxic", "Oxy_cat_Hypoxic_vs_Normoxic"), type="ashr")

table(res_es_HvA_ashr$padj < 0.1)
#same number as above :)

#more about contrasts here: https://github.com/tavareshugo/tutorial_DESeq2_contrasts/blob/main/DESeq2_contrasts.md
#also includes instructions for combining levels using a matrix structure


#intercept does not work for normoxic vs. hypoxic... how do we get these data?!
#checked https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
#install.packages('pheatmap')

#prep for PCA/heatmaps
#transformation 
vsd <- vst(dds_es, blind=FALSE)

df <- as.data.frame(colData(dds_es)[,c("Oxy_cat")])

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
plotPCA(vsd, intgroup=c('Oxy_cat'))

#lots of separation from one normoxic sample
#some separation with anoxia and normoxia..
#maybe would be better to use continuous oxygen rather than categorical

#could be effectively the same if top left outlier is removed from the PCA plot 

#MA plots

plotMA(res_es_AvN)

plotMA(res_es_AvN_apeglm)

plotMA(res_es_HvN_apeglm)

plotMA(res_es_HvA_ashr)

##After calling plotMA, one can use the function identify to interactively detect the row number of individual genes by clicking on the plot. One can then recover the gene identifiers by saving the resulting indices
#idx <- identify(res_es_AvN_apeglm$baseMean, res_es_AvN_apeglm$log2FoldChange)
##print gene names from points you clicked on 
#rownames(res_es)[idx]

####possible removing of outliers 
boxplot(log10(assays(dds_es)[["cooks"]]), range=0, las=2)

#anoxic samples both look like outliers, but biologically are not


#plotting dispersion estimates
plotDispEsts(dds_es)


#PCA to remove outlier
PCA <- plotPCA(vsd, intgroup=c('Oxy_cat'), returnData=TRUE)
PCA

#DC33 is the outlier based in the PCA, distorting the whole view, so we can try pulling it out

###########################
##removing 33
#############################
#because it removes a whole factor, have to change the design matrix or reload the data
DESeq2table_es2 <- DESeq2table_es
#drop normoxic as level 
#DESeq2table_es2$Oxy_cat <- droplevels.factor(DESeq2table_es2$Oxy_cat, "Normoxic")

#relevel with hypoxia as reference 
#DESeq2table_es2$Oxy_cat <- relevel(DESeq2table_es2$Oxy_cat, "Hypoxic")

DESeq2table_es_no_out2 <- DESeq2table_es2[, !(colnames(DESeq2table_es2) %in% c('DC33'))]
dds_es2 <- DESeq(DESeq2table_es_no_out2)

resultsNames(dds_es2) # lists the coefficients
res_es2_AvN <- results(dds_es2, name="Oxy_cat_Anoxic_vs_Normoxic")
# or to shrink log fold changes association with condition:

#shrunken results beased on LFC
res_es2_AvN_apeglm <- lfcShrink(dds_es2, coef="Oxy_cat_Anoxic_vs_Normoxic", type="apeglm")

#significantly expressed genes based on full results
table(res_es2_AvN$padj < 0.1)

#FALSE  TRUE 
#29086  1376 
#1376 differentially expressed genes between anoxic and normoxic samples :)

#significantly expressed genes based on LFC results
table(res_es2_AvN_apeglm$padj < 0.1)
##need to run the other multiple comparisons
#standard method - hypoxic vs normoxic
res_es2_HvN <- results(dds_es2, name="Oxy_cat_Hypoxic_vs_Normoxic")

table(res_es2_HvN$padj < 0.1)
#FALSE  TRUE 
#30483   737 
#737 differentially expressed genes between hypoxic and normoxic

# or to shrink log fold changes association with condition:
res_es2_HvN_apeglm <- lfcShrink(dds_es2, coef="Oxy_cat_Hypoxic_vs_Normoxic", type="apeglm")

table(res_es2_HvN_apeglm$padj < 0.1)
#same as above

length(res_es2_AvN[,1])

#44249 genes in total for es

#comparison (contrast), between hypoxic and anoxic samples
res_es2_HvA <- results(dds_es2, contrast = list('Oxy_cat_Hypoxic_vs_Normoxic', 'Oxy_cat_Anoxic_vs_Normoxic'))

table(res_es2_HvA$padj < 0.1)
#FALSE  TRUE 
#24032  1090

res_es2_HvA_ashr <- lfcShrink(dds_es2, contrast = list("Oxy_cat_Anoxic_vs_Normoxic", "Oxy_cat_Hypoxic_vs_Normoxic"), type="ashr")

table(res_es2_HvA_ashr$padj < 0.1)

#prep for PCA/heatmaps
#transformation 
vsd_es <- vst(dds_es2, blind=FALSE)

#sample to sample distances
sampleDists <- dist(t(assay(vsd_es)))

#heatmap of distance matrices between samples 
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(vsd_es), vsd_es$Oxy_cat, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#not all anoxic samples group together. One anoxic sample and one hypoxic sample are very close to eachother 

#pca plot using EDA seq
plotPCA(vsd_es, intgroup=c('Oxy_cat'))


#more genes (top 10%)
es_pca <- plotPCA(vsd_es, intgroup=c('Oxy_cat'), ntop=4423)
#like HS, makes sep worse in anoxia 

es_pca_data <- plotPCA(vsd_es, intgroup=c('Oxy_cat'), ntop=4423, returnData=TRUE)

#get right colors, usuing ggplot
es_pca_final <- ggplot(es_pca_data,aes(x=PC1,y=PC2,col=Oxy_cat))+
  geom_point(size=4,alpha=1)+
  scale_color_manual(values = c("#56B4E9",'#D55E00',"#E69F00"))+
  xlab("PC1: 15% variance") + ylab("PC2: 18% variance") +
  theme_classic()

es_pca_final

#check with final layout
Figure_sponges <- ggarrange(es_pca_final, 
                        labels = c('A', 'B', 'C', 'D', 'E'),
                        common.legend = TRUE,
                        ncol = 2, nrow = 3, 
                        legend = 'bottom')
Figure_sponges

#fewer genes 
#plotPCA(vsd_es, intgroup=c('Oxy_cat'), ntop=100)
#looks the same

#MA plots

#save PCA as 

plotMA(res_es2_AvN)

plotMA(res_es2_AvN_apeglm)

plotMA(res_es2_HvA_ashr)

plotMA(res_es2_HvN_apeglm)
##After calling plotMA, one can use the function identify to interactively detect the row number of individual genes by clicking on the plot. One can then recover the gene identifiers by saving the resulting indices
#idx <- identify(res_es_AvN_apeglm$baseMean, res_es_AvN_apeglm$log2FoldChange)
##print gene names from points you clicked on 
#rownames(res_es)[idx]

####possible removing of outliers 
boxplot(log10(assays(dds_es2)[["cooks"]]), range=0, las=2)

#anoxic samples still look like outliers, but we need to keep them in, plus they are simlar to eachother, so likely represent the anoxic population

#plotting dispersion estimates
plotDispEsts(dds_es2)

#plot based on minimum p value
plotCounts(dds_es2, gene=which.min(res_es2_AvN$padj), intgroup=c('Oxy_cat'))

#diagnostic histogram 
hist(res_es2_AvN$pvalue[res_es2_AvN$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")

#has same freguency with high on the left :) yay
#see http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html


#########
#summary of all degs using package ReportingTools
#########

#add annotations
#work around for getting longer gene names (not just isogroup)
gn=read.table("E_sp2_sponge_ref_transcriptome_iso2geneName_final.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)

#gene_id
#gene_name

names(gn)[names(gn) == "V1"] <- "gene_id"

names(gn)[names(gn) == "V2"] <- "gene_name"

head(gn)



#add annotation to dds_es object 
#check if they match in order 
#all(rownames(dds_es2) == gn$gene_id)
#false, they do not match (expected)

#have to match them
#gn2 <- gn[match(rownames(dds_es2), gn$gene_id),]
#doesn't work because many genes are not annotated

#try with left join? 
#need to get a count table with a column name for the row data to count 
counts_esHost2 <- cbind(rownames(counts_es), counts_es)
rownames(counts_esHost2) <- NULL
colnames(counts_esHost2) <- c("gene_id", 'DC108', 'DC110', 'DC112', 'DC114', 'DC118', 'DC123', 'DC131', 'DC135', 'DC136', 'DC24', 'DC33', 'DC39', 'DC56', 'DC85', 'DC93', 'DC97')

joined <- left_join(counts_esHost2, gn, by = "gene_id")
joined$full <-paste(joined$gene_id,joined$gene_name)

#reorder so that combined column is first (also dropp extra columns not needed for count matrix)
#also got rid of names column
joined2 <-joined[, c(19, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17)]

#make first column back into rowname
joined3 <- joined2[,-1]
rownames(joined3) <- joined2[,1]

#repeat deseq anaylsis with annotated count table 
counts_es=joined3[,order(names(joined3))]
names(counts_es)

#replace all spaces in row names (gene names) with underscore
rownames(counts_es) <- gsub(" ", "_", rownames(counts_es))
#replace all special characters
rownames(counts_es) <- gsub("[()]", "_", rownames(counts_es))
rownames(counts_es) <- gsub("[:]", "_", rownames(counts_es))
rownames(counts_es) <- gsub("[/]", "_", rownames(counts_es))
rownames(counts_es) <- gsub("[.]", "_", rownames(counts_es))
rownames(counts_es) <- gsub("[,]", "_", rownames(counts_es))

#limit rownames to 100 characters
rownames(counts_es) <- strtrim(rownames(counts_es), 50)

##conditions (metadata)
#conditions=read.table("Conditions_es.txt", header = T)

#conds_es <- conditions[ order(row.names(conditions)), ]

#adding longer names
####load data into deseq table format
#includes experimental design' since it is from the quick start method
DESeq2table_es_names <- DESeqDataSetFromMatrix(countData = counts_es,
                                               colData = conds_es,
                                               design= ~ Oxy_cat)

DESeq2table_es_names$Oxy_cat <- relevel(DESeq2table_es_names$Oxy_cat, "Normoxic")


DESeq2table_es_no_out2 <- DESeq2table_es_names[, !(colnames(DESeq2table_es_names) %in% c('DC33'))]
dds_es3 <- DESeq(DESeq2table_es_no_out2)

#with releveling, of course it doesn't fit
resultsNames(dds_es3)


#hypoxia vs. anoxia
res_es3_HvA <- results(dds_es3, contrast = list('Oxy_cat_Hypoxic_vs_Normoxic', 'Oxy_cat_Anoxic_vs_Normoxic'))

table(res_es3_HvA$padj < 0.1)

#same as above once releveled with normoxia as baseline.
#FALSE  TRUE 
#24032  1090 
res_es3_HvA_ashr <- lfcShrink(dds_es3, contrast = list("Oxy_cat_Anoxic_vs_Normoxic", "Oxy_cat_Hypoxic_vs_Normoxic"), type="ashr")

table(res_es3_HvA_ashr$padj < 0.1)
#same as above :)

#Normoxia vs anoxia
res_es3_AvN <- results(dds_es3, name="Oxy_cat_Anoxic_vs_Normoxic")

table(res_es3_AvN$padj < 0.1)
#29086  1376 
#1376 differentially expressed genes between anoxic and normoxic samples :)

#shrunken results beased on LFC
res_es3_AvN_apeglm <- lfcShrink(dds_es3, coef="Oxy_cat_Anoxic_vs_Normoxic", type="apeglm")

#significantly expressed genes based on full results
table(res_es3_AvN_apeglm$padj < 0.1)

#FALSE  TRUE 
#29086  1376
#1376 differentially expressed genes between anoxic and normoxic samples :)

#significantly expressed genes based on LFC results
table(res_es3_AvN_apeglm$padj < 0.1)
##need to run the other multiple comparisons
#standard method - hypoxic vs normoxic
res_es3_HvN <- results(dds_es3, name="Oxy_cat_Hypoxic_vs_Normoxic")

table(res_es3_HvN$padj < 0.1)
#FALSE  TRUE 
#30483   737  
#737 differentially expressed genes between hypoxic and normoxic

# or to shrink log fold changes association with condition:
res_es3_HvN_apeglm <- lfcShrink(dds_es3, coef="Oxy_cat_Hypoxic_vs_Normoxic", type="apeglm")

table(res_es3_HvN_apeglm$padj < 0.1)

######################
#make report with annotations
######################

des2Report <- HTMLReport(shortName = 'ES sponge RNAseq_analysis_with_DESeq2__HvN', title = 'ES_sponge_HvN_names_no33',reportDirectory = "./reports")
publish(dds_es3,des2Report, pvalueCutoff=0.1, factor = colData(dds_es3)$Oxy_cat, reportDir="./reports")
finish(des2Report)

#es HvN
#log fold change less than 0 -> downregulated in hypoxia compared to normoxia
#log fold change greater than 0 - upregulated in hypoxia

#try with specific contrasts 
#anoxia vs normoxia
#has a built in limit of 1000
des2Report <- HTMLReport(shortName = 'ES_sponge_RNAseq_analysis_with_DESeq2_AvN', title = 'ES_sponge_AvN_names_no33',reportDirectory = "./reports")
publish(dds_es3,des2Report, pvalueCutoff=0.1, n=1500, factor = colData(dds_es3)$Oxy_cat, contrast =c("Oxy_cat", "Normoxic", "Anoxic"), reportDir="./reports")
finish(des2Report)

#es AvN
#log fold change less than 0 - upregulated in anoxia, compared to normoxia
#log fold change greater than 0 - downregulated in anoxia

#anoxia vs hypoxia
#has a built in limit of 1000
des2Report <- HTMLReport(shortName = 'ES_sponge_analysis_with_DESeq2_AvH', title = 'EES_sponge_AvH_names_no33',reportDirectory = "./reports")
publish(dds_es3,des2Report, pvalueCutoff=0.1, n=1500, factor = colData(dds_es3)$Oxy_cat, contrast =c("Oxy_cat", "Anoxic", "Hypoxic"), reportDir="./reports")
finish(des2Report)


#identify number of genes upregulated
#dataframe with differentially expressed genes
#Hypoxia vs. anoxia
df1 <-data.frame(res_es3_HvA)
df2 <- df1 %>% 
  filter(padj < 0.1)

table(df2$padj < 0.1)
#1090  true :)

#upregulated in anoxia (log fold change greater than 0)
#downregulated in anoxia (log fold change less than zero)

#up
df_up <- df2 %>% 
  filter(log2FoldChange > 0)

nrow(df_up)
#622 upregulated

#how many unannotated?
length(grep("_NA", rownames(df_up), fixed=TRUE))

#349 unannotated
nrow(df_up)-length(grep("_NA", rownames(df_up), fixed=TRUE))

#273 genes annotated and upregulated in anoxia vs hypoxia in es
#dataframe with only these genes 
df_up_an <- df_up[- grep("_NA", rownames(df_up)),]

#downregulated genes
df_down <- df2 %>% 
  filter(log2FoldChange < 0)

nrow(df_down)
#468 genes downregulated

#check with total
nrow(df_down)+nrow(df_up)
#1090 :)

#how many unannotated?
length(grep("_NA", rownames(df_down), fixed=TRUE))
#263 genes unannotated and downregulated 

#annotated, down reg
df_down_an <- df_down[- grep("_NA", rownames(df_down)),]

nrow(df_down_an)
#205 genes annotated and down regulated in anoxia vs hypoxia

#NvA
#+ = lower in anoxia
#- = higher in anoxia
df1 <-data.frame(res_es3_AvN)
df2 <- df1 %>% 
  filter(padj < 0.1)

table(df2$padj < 0.1)
#1376  true :)

#upregulated in anoxia
df_up <- df2 %>% 
  filter(log2FoldChange < 0)

nrow(df_up)
#551 upregulated

#how many unannotated?
length(grep("_NA", rownames(df_up), fixed=TRUE))

#323 unannotated
nrow(df_up)-length(grep("_NA", rownames(df_up), fixed=TRUE))

#228 genes annotated and upregulated in normoxia vs anoxia in es
#dataframe with only these genes 
df_up_an <- df_up[- grep("_NA", rownames(df_up)),]

#downregulated genes
df_down <- df2 %>% 
  filter(log2FoldChange > 0)

nrow(df_down)
#825 genes downregulated

#check with total
nrow(df_down)+nrow(df_up)
#1376 :)

#how many unannotated?
length(grep("_NA", rownames(df_down), fixed=TRUE))
#485 genes unannotated and downregulated in normoxia vs anoxia 

#annotated, down reg
df_down_an <- df_down[- grep("_NA", rownames(df_down)),]

nrow(df_down_an)
#340 genes annotated and down regulated in normoxia vs anoxia

#######NvH
#+ = up in hypoxia
#=- = Up in normoxia
df1 <-data.frame(res_es3_HvN)
df2 <- df1 %>% 
  filter(padj < 0.1)

table(df2$padj < 0.1)
#737 true :)

#upregulated in hypoxia
df_up <- df2 %>% 
  filter(log2FoldChange > 0)

nrow(df_up)
#499 upregulated

#how many unannotated?
length(grep("_NA", rownames(df_up), fixed=TRUE))

#260 unannotated
nrow(df_up)-length(grep("_NA", rownames(df_up), fixed=TRUE))

#239 genes annotated and upregulated in hypoxia vs normoxia in es
#dataframe with only these genes 
df_up_an <- df_up[- grep("_NA", rownames(df_up)),]

#downregulated genes
df_down <- df2 %>% 
  filter(log2FoldChange < 0)

nrow(df_down)
#238 genes downregulated

#check with total
nrow(df_down)+nrow(df_up)
#126 :)

#how many unannotated?
length(grep("_NA", rownames(df_down), fixed=TRUE))
#126 genes unannotated and downregulated in hypoxia vs normoxia 

#annotated, down reg
df_down_an <- df_down[- grep("_NA", rownames(df_down)),]

nrow(df_down_an)
#112 genes annotated and down regulated in normoxia vs anoxia

#prepared csv for GO analysis using go_mwu https://github.com/z0on/GO_MWU
#df2 has all differentially expressed genes
#make log p, signed based on up or down regulation 
#need to use dataframe with only isogroup designation (dds_es2)

#need to take out filter step 

df1 <-data.frame(res_es2_HvA)
df2 <- df1 

table(df2$padj < 0.1)
#avh LFC >0 -> upregulated in anoxia
df2$logP <- ifelse(df2$log2FoldChange > 0, -log(df2$pvalue), log(df2$pvalue))

#make csv with just gene names and logP
go <- data.frame(rownames(df2))
go$logP <-df2$logP 
names(go)[1] <- 'gene_id'

#write as csv
write.csv(go, file='es_go_HvA_logP_final.csv', row.names = FALSE)

#GO HvN
df1 <-data.frame(res_es2_HvN)
df2 <- df1 

table(df2$padj < 0.1)
#LFC >0 -> upregulated in hypoxia
df2$logP <- ifelse(df2$log2FoldChange > 0, -log(df2$pvalue), log(df2$pvalue))

#make csv with just gene names and logP
go <- data.frame(rownames(df2))
go$logP <-df2$logP 
names(go)[1] <- 'gene_id'

#write as csv
write.csv(go, file='es_go_HvN_logP_final.csv', row.names = FALSE)

#AvN
df1 <-data.frame(res_es2_AvN)
df2 <- df1 

table(df2$padj < 0.1)
#LFC >0 -> upregulated in normoxia (downreg in anoxia)
df2$logP <- ifelse(df2$log2FoldChange < 0, -log(df2$pvalue), log(df2$pvalue))

#make csv with just gene names and logP
go <- data.frame(rownames(df2))
go$logP <-df2$logP 
names(go)[1] <- 'gene_id'

#write as csv
write.csv(go, file='es_go_AvN_logP_final.csv', row.names = FALSE)


######################
#write csvs of results
df1 <-data.frame(res_es2_HvA)
write.csv(df1, file='es_degs_DESeq2_res_HvA_final.csv', sep=',')

df1 <-data.frame(res_es2_HvA_ashr)
write.csv(df1, file='es_degs_DESeq2_res_HvA_ashr_final.csv', sep=',')

df1 <-data.frame(res_es2_HvN)
write.csv(df1, file='es_degs_DESeq2_res_HvN_final.csv', sep=',')

df1 <-data.frame(res_es2_HvN_apeglm)
write.csv(df1, file='es_degs_DESeq2_res_HvN_apeglm_final.csv', sep=',')

df1 <-data.frame(res_es2_AvN)
write.csv(df1, file='es_degs_DESeq2_res_AvN_final.csv', sep=',')

df1 <-data.frame(res_es2_AvN_apeglm)
write.csv(df1, file='es_degs_DESeq2_res_AvN_apeglm_final.csv', sep=',')

####################### KogMWU setup - can use signed p values (e.g.Kenkel, C. D., Mocellin, V. J. L., & Bay, L. K. (2020). Global gene expression patterns in Porites white patch syndrome: Disentangling symbiont loss from the thermal stress response in reef-building coral. Molecular Ecology, 29(20), 3907-3920. https://doi.org/10.1111/mec.15608) or log fold change 
#already have signed p values (above)

#make files with log fold changes 
df1 <-data.frame(res_es2_HvA)
df2 <- df1 

table(df2$padj < 0.1)
#avh LFC >0 -> upregulated in anoxia
df2$lfc <- df2$log2FoldChange

#make csv with just gene names and lfc
kog <- data.frame(rownames(df2))
kog$lfc <-df2$lfc 
names(kog)[1] <- 'gene_id'

#write as csv
write.csv(kog, file='es_kog_HvA_lfc_final.csv', row.names = FALSE)

#GO HvN
df1 <-data.frame(res_es2_HvN)
df2 <- df1 

table(df2$padj < 0.1)
#LFC >0 -> upregulated in hypoxia
df2$lfc <- df2$log2FoldChange

#make csv with just gene names and lfc
kog <- data.frame(rownames(df2))
kog$lfc <-df2$lfc 
names(kog)[1] <- 'gene_id'

#write as csv
write.csv(kog, file='es_kog_HvN_lfc_final.csv', row.names = FALSE)

#AvN
df1 <-data.frame(res_es2_AvN)
df2 <- df1 

table(df2$padj < 0.1)
#LFC >0 -> upregulated in normoxia (downreg in anoxia)
df2$lfc <- df2$log2FoldChange*-1

#make csv with just gene names and lfc
kog <- data.frame(rownames(df2))
kog$lfc <-df2$lfc 
names(kog)[1] <- 'gene_id'

#write as csv
write.csv(kog, file='es_kog_AvN_lfc_final.csv', row.names = FALSE)


#modified from example: https://rdrr.io/cran/KOGMWU/man/kog.mwu.html
# Import KOG table 
es_kog2gene=read.table("E_sp2_sponge_ref_transcriptome_iso2kogClass_final.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)

#import signed p values for AvN
es_AvN.signedp <- read.csv('es_go_AvN_logP_final.csv', header =TRUE)
#lfc for comparison
es_AvN.lfc <- read.csv('es_kog_AvN_lfc_final.csv', header =TRUE)

es_AvN.lth=kog.mwu(es_AvN.lfc,es_kog2gene) 
es_AvN.lth 

#import signed p values for HvN
#es_HvN.signedp <- read.csv('es_go_HvN_logP_final.csv', header =TRUE)
#lfc (can switch them out to see the diffrences)
es_HvN.lfc <- read.csv('es_kog_HvN_lfc_final.csv', header =TRUE)
es_HvN.lth=kog.mwu(es_HvN.lfc,es_kog2gene)
es_HvN.lth

# es HvA
es_HvA.signedp <- read.csv('es_go_HvA_logP_final.csv', header =TRUE)
es_HvA.lfc <- read.csv('es_kog_HvA_lfc_final.csv', header =TRUE)
es_HvA.lth=kog.mwu(es_HvA.lfc,es_kog2gene)
es_HvA.lth

# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("Anoxia"=es_AvN.lth,"Hypoxia"=es_HvN.lth,"Hypoxia vs Anoxia"=es_HvA.lth))

# Making a heatmap with hierarchical clustering trees: 
pheatmap(as.matrix(ktable),clustering_distance_cols="correlation") 

# exploring correlations between datasets
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor)
# p-values of these correlations in the upper panel:
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor.pval)

# plotting individual delta-rank correlations:
corrPlot(x="Anoxia",y="Hypoxia",ktable)
corrPlot(x="Anoxia",y="Hypoxia vs Anoxia",ktable)
corrPlot(x="Hypoxia",y="Hypoxia vs Anoxia",ktable)

#what doe the various correlations mean?
#responses to anoxia and hypoxia are distinct.



#######################
#Hymeraphia stellifera, sponge host

countsHost=read.table("allcounts_hs_sponge_final.txt",header=TRUE,row.names=1) #Reading in the table of counts per isogroup by sample - you must open R in the same folder or change the working directory

counts=countsHost[,order(names(countsHost))]
names(counts)

#add annotations
gn=read.table("H_stellifera_sponge_ref_transcriptome_iso2geneName_final.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)

names(gn)[names(gn) == "V1"] <- "gene_id"

names(gn)[names(gn) == "V2"] <- "gene_name"

head(gn)

#need to get a count table with a column name for the row data to count 
counts_hsHost2 <- cbind(rownames(counts), counts)
rownames(counts_hsHost2) <- NULL
colnames(counts_hsHost2) <- c("gene_id", "DC105", "DC107", "DC116", "DC121", "DC129", "DC133", "DC54",  "DC55")

joined <- left_join(counts_hsHost2, gn, by = "gene_id")
joined$full <-paste(joined$gene_id,joined$gene_name)

#reorder so that combined column is first (also dropp extra columns not needed for count matrix)
#also got rid of names column
joined2 <-joined[, c(11, 2, 3, 4, 5, 6, 7, 8, 9)]

#make first column back into rowname
joined3 <- joined2[,-1]
rownames(joined3) <- joined2[,1]

#repeat deseq anaylsis with annotated count table 
counts_hs=joined3[,order(names(joined3))]
names(counts_hs)

#replace all spaces in row names (gene names) with underscore
rownames(counts_hs) <- gsub(" ", "_", rownames(counts_hs))
#replace all special characters
rownames(counts_hs) <- gsub("[()]", "_", rownames(counts_hs))
rownames(counts_hs) <- gsub("[:]", "_", rownames(counts_hs))
rownames(counts_hs) <- gsub("[/]", "_", rownames(counts_hs))
rownames(counts_hs) <- gsub("[.]", "_", rownames(counts_hs))
rownames(counts_hs) <- gsub("[,]", "_", rownames(counts_hs))

#limit rownames to 100 characters
rownames(counts_hs) <- strtrim(rownames(counts_hs), 50)

##conditions (metadata)
conditions=read.table("Conditions_hs.txt", header = T)

conds <- conditions[ order(row.names(conditions)), ]

####load data into deseq table format
#includes experimental design' since it is from the quick start method
DESeq2Table <- DESeqDataSetFromMatrix(countData = counts_hs,
                              colData = conds,
                              design= ~ Oxy_cat)

#make normoxic the 'reference' level 
DESeq2Table$Oxy_cat <- relevel(DESeq2Table$Oxy_cat, "Normoxic")

#runs differential expression analysis quickly
dds <- DESeq(DESeq2Table)

resultsNames(dds) # lists the coefficients
res_AvN <- results(dds, name="Oxy_cat_Anoxic_vs_Normoxic")
# or to shrink log fold changes association with condition:

#shrunken results beased on LFC
res_AvN_apeglm <- lfcShrink(dds, coef="Oxy_cat_Anoxic_vs_Normoxic", type="apeglm")

#significantly expressed genes based on full results
table(res_AvN$padj < 0.1)

#FALSE  TRUE 
#23904   284 
#284 differentially expressed genes between normoxic and anoxic samples :)

#significantly expressed genes based on LFC results
table(res_AvN_apeglm$padj < 0.1)

## differentially expressed genes between normoxic and anoxic samples = same as above :)

#standard method - hypoxic vs normoxic
res_HvN <- results(dds, name="Oxy_cat_Hypoxic_vs_Normoxic")

table(res_HvN$padj < 0.1)
#FALSE  TRUE 
#19646   523 
#523 differentially expressed genes between hypoxic and normoxic

# or to shrink log fold changes association with condition:
res_HvN_apeglm <- lfcShrink(dds, coef="Oxy_cat_Hypoxic_vs_Normoxic", type="apeglm")

table(res_HvN_apeglm$padj < 0.1)
#same as above

length(res_HvN[,1])

#41633 genes in total for HS

#comparison (contrast), between hypoxic and anoxic samples
res_HvA <- results(dds, contrast = list('Oxy_cat_Hypoxic_vs_Normoxic', 'Oxy_cat_Anoxic_vs_Normoxic'))

table(res_HvA$padj < 0.1)
#FALSE  TRUE 
#22841  2154 

#2154 genes differentially expressed between hypoxic and anoxic

#shrunken results 
#apeglm does not have contrast ability
#need to try 'normal shrinkage, or ashr

res_HvA_ashr <- lfcShrink(dds, contrast = list("Oxy_cat_Anoxic_vs_Normoxic", "Oxy_cat_Hypoxic_vs_Normoxic"), type="ashr")

table(res_HvA_ashr$padj < 0.1)
 
#same number as above :)

#more about contrasts here: https://github.com/tavareshugo/tutorial_DESeq2_contrasts/blob/main/DESeq2_contrasts.md
#also includes instructions for combining levels using a matrix structure


#intercept does not work for normoxic vs. hypoxic... how do we get these data?!
#checked https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
#install.packages('pheatmap')

#prep for PCA/heatmaps
#transformation 
vsd <- vst(dds, blind=FALSE)

df <- as.data.frame(colData(dds)[,c("Oxy_cat")])

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

#pca plot with side by side comparison to eurypon
plotPCA(vsd, intgroup=c('Oxy_cat'))


#more genes (top 10%)
hs_pca <- plotPCA(vsd, intgroup=c('Oxy_cat'), ntop=4423)
#like HS, makes sep worse in anoxia 

hs_pca_data <- plotPCA(vsd, intgroup=c('Oxy_cat'), ntop=4423, returnData=TRUE)

#get right colors, usuing ggplot
hs_pca_final <- ggplot(hs_pca_data,aes(x=PC1,y=PC2,col=Oxy_cat))+
  geom_point(size=4,alpha=1)+
  scale_color_manual(values = c("#56B4E9",'#D55E00',"#E69F00"))+
  xlab("PC1: 32% variance") + ylab("PC2: 17% variance") +
  theme_classic()

hs_pca_final

#check with final layout
Figure_sponges <- ggarrange(es_pca_final, hs_pca_final,
                            labels = c('A', 'B'),
                            common.legend = TRUE,
                            ncol = 2, nrow = 1, 
                            legend = 'bottom')
Figure_sponges

#good separation of the anoxic sample from the hypoxic and normoxic ones
#unfortunately only one normoxic sample :(
#however, given that the sponges were behaving normally under hypoxia (pumping), hypoxic could be pretty close to normal
#so hypoxic dot on the left could be an outlier...

#could be effectively the same if top left outlier is removed from the PCA plot 

#MA plots

plotMA(res_AvN)

plotMA(res_AvN_apeglm)

plotMA(res_HvN_apeglm)

plotMA(res_HvA_ashr)

##After calling plotMA, one can use the function identify to interactively detect the row number of individual genes by clicking on the plot. One can then recover the gene identifiers by saving the resulting indices
#idx <- identify(res_AvN_apeglm$baseMean, res_AvN_apeglm$log2FoldChange)
##print gene names from points you clicked on 
#rownames(res)[idx]

####possible removing of outliers 
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)

#sample DC129 has consistantly higher outliers 

#plotting dispersion estimates
plotDispEsts(dds)

#get reports
#without annotations (for now, just to check directions)


des2Report <- HTMLReport(shortName = 'HS sponge RNAseq_analysis_with_DESeq2_HvN_all_data', title = 'HS_sponge_HvN_',reportDirectory = "./reports")
publish(dds,des2Report, pvalueCutoff=0.1, factor = colData(dds)$Oxy_cat, contrast =c("Oxy_cat", "Normoxic", "Hypoxic"), reportDir="./reports")
finish(des2Report)


des2Report <- HTMLReport(shortName = 'HS_sponge_RNAseq_analysis_with_DESeq2_AvN_all_data', title = 'HS_sponge_AvN_names',reportDirectory = "./reports")
publish(dds,des2Report, pvalueCutoff=0.1, n=1500, factor = colData(dds)$Oxy_cat, contrast =c("Oxy_cat", "Normoxic", "Anoxic"), reportDir="./reports")
finish(des2Report)


des2Report <- HTMLReport(shortName = 'HS_sponge_analysis_with_DESeq2_AvH_all_data', title = 'HS_sponge_AvH_names',reportDirectory = "./reports")
publish(dds,des2Report, pvalueCutoff=0.1, n=2600, factor = colData(dds)$Oxy_cat, contrast =c("Oxy_cat", "Anoxic", "Hypoxic"), reportDir="./reports")
finish(des2Report)

#up and downregulated genes hs 
#Hypoxia vs. anoxia
df1 <-data.frame(res_HvA)
df2 <- df1 %>% 
  filter(padj < 0.1)

table(df2$padj < 0.1)
#2154  true :)

#upregulated in anoxia (log fold change greater than 0)
#downregulated in anoxia (log fold change less than zero)

#up
df_up <- df2 %>% 
  filter(log2FoldChange > 0)

nrow(df_up)
#1531 upregulated

#how many unannotated?
length(grep("_NA", rownames(df_up), fixed=TRUE))

#1093 unannotated
nrow(df_up)-length(grep("_NA", rownames(df_up), fixed=TRUE))

#438 genes annotated and upregulated in anoxia vs hypoxia in es
#dataframe with only these genes 
df_up_an <- df_up[- grep("_NA", rownames(df_up)),]

#downregulated genes
df_down <- df2 %>% 
  filter(log2FoldChange < 0)

nrow(df_down)
#623 genes downregulated

#check with total
nrow(df_down)+nrow(df_up)
#2154 :)

#how many unannotated?
length(grep("_NA", rownames(df_down), fixed=TRUE))
#337 genes unannotated and downregulated 

#annotated, down reg
df_down_an <- df_down[- grep("_NA", rownames(df_down)),]

nrow(df_down_an)
#286 genes annotated and down regulated in anoxia vs hypoxia

#NvA
#+ = lower in anoxia
#- = higher in anoxia
df1 <-data.frame(res_AvN)
df2 <- df1 %>% 
  filter(padj < 0.1)

table(df2$padj < 0.1)
#284  true :)

#upregulated in anoxia
df_up <- df2 %>% 
  filter(log2FoldChange < 0)

nrow(df_up)
#40 upregulated

#how many unannotated?
length(grep("_NA", rownames(df_up), fixed=TRUE))

#29 unannotated
nrow(df_up)-length(grep("_NA", rownames(df_up), fixed=TRUE))

#11 genes annotated and upregulated in normoxia vs anoxia in es
#dataframe with only these genes 
df_up_an <- df_up[- grep("_NA", rownames(df_up)),]

#downregulated genes
df_down <- df2 %>% 
  filter(log2FoldChange > 0)

nrow(df_down)
#224 genes downregulated

#check with total
nrow(df_down)+nrow(df_up)
#284 :)

#how many unannotated?
length(grep("_NA", rownames(df_down), fixed=TRUE))
#157 genes unannotated and downregulated in normoxia vs anoxia 

#annotated, down reg
df_down_an <- df_down[- grep("_NA", rownames(df_down)),]

nrow(df_down_an)
#87 genes annotated and down regulated in normoxia vs anoxia

#######NvH
#+ = up in hypoxia
#=- = Up in normoxia
df1 <-data.frame(res_HvN)
df2 <- df1 %>% 
  filter(padj < 0.1)

table(df2$padj < 0.1)
#523 true :)

#upregulated in hypoxia
df_up <- df2 %>% 
  filter(log2FoldChange > 0)

nrow(df_up)
#465 upregulated

#how many unannotated?
length(grep("_NA", rownames(df_up), fixed=TRUE))

#340 unannotated
nrow(df_up)-length(grep("_NA", rownames(df_up), fixed=TRUE))

#125 genes annotated and upregulated in hypoxia vs normoxia in es
#dataframe with only these genes 
df_up_an <- df_up[- grep("_NA", rownames(df_up)),]

#downregulated genes
df_down <- df2 %>% 
  filter(log2FoldChange < 0)

nrow(df_down)
#58 genes downregulated

#check with total
nrow(df_down)+nrow(df_up)
#523 :)

#how many unannotated?
length(grep("_NA", rownames(df_down), fixed=TRUE))
#49 genes unannotated and downregulated in hypoxia vs normoxia 

#annotated, down reg
df_down_an <- df_down[- grep("_NA", rownames(df_down)),]

nrow(df_down_an)
#9 genes annotated and down regulated in normoxia vs anoxia


#KOGMWU with all samples included for hs

#kog_mwu with all HS data based on log fold change
#make files with log fold changes 
df1 <-data.frame(res_HvA)
df2 <- df1 

table(df2$padj < 0.1)
#avh LFC >0 -> upregulated in anoxia
df2$lfc <- df2$log2FoldChange

#make csv with just gene names and lfc
kog <- data.frame(rownames(df2))
kog$lfc <-df2$lfc 
names(kog)[1] <- 'gene_id'

#write as csv
write.csv(kog, file='hs_kog_HvA_lfc_final.csv', row.names = FALSE)

#kog HvN
df1 <-data.frame(res_HvN)
df2 <- df1 

table(df2$padj < 0.1)
#LFC >0 -> upregulated in hypoxia
df2$lfc <- df2$log2FoldChange

#make csv with just gene names and lfc
kog <- data.frame(rownames(df2))
kog$lfc <-df2$lfc 
names(kog)[1] <- 'gene_id'

#write as csv
write.csv(kog, file='hs_kog_HvN_lfc_final.csv', row.names = FALSE)

#AvN
df1 <-data.frame(res_AvN)
df2 <- df1 

table(df2$padj < 0.1)
#LFC >0 -> upregulated in normoxia (downreg in anoxia)
df2$lfc <- df2$log2FoldChange*-1

#make csv with just gene names and lfc
kog <- data.frame(rownames(df2))
kog$lfc <-df2$lfc 
names(kog)[1] <- 'gene_id'

#write as csv
write.csv(kog, file='hs_kog_AvN_lfc_final.csv', row.names = FALSE)


#modified from example: https://rdrr.io/cran/KOGMWU/man/kog.mwu.html
#Es KOGS
es_kog2gene=read.table("E_sp2_sponge_ref_transcriptome_iso2kogClass_final.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)

#import signed p values for AvN
es_AvN.signedp <- read.csv('es_go_AvN_logP_final.csv', header =TRUE)
#lfc for comparison
es_AvN.lfc <- read.csv('es_kog_AvN_lfc_final.csv', header =TRUE)

es_AvN.lth=kog.mwu(es_AvN.lfc,es_kog2gene) 
es_AvN.lth 

#import signed p values for HvN
#es_HvN.signedp <- read.csv('es_go_HvN_logP_final.csv', header =TRUE)
#lfc (can switch them out to see the diffrences)
es_HvN.lfc <- read.csv('es_kog_HvN_lfc_final.csv', header =TRUE)
es_HvN.lth=kog.mwu(es_HvN.lfc,es_kog2gene)
es_HvN.lth

# es HvA
es_HvA.signedp <- read.csv('es_go_HvA_logP_final.csv', header =TRUE)
es_HvA.lfc <- read.csv('es_kog_HvA_lfc_final.csv', header =TRUE)
es_HvA.lth=kog.mwu(es_HvA.lfc,es_kog2gene)
es_HvA.lth


# Import KOG table Hymeraphia 

hs_kog2gene=read.table("H_stellifera_sponge_ref_transcriptome_iso2kogClass_final.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)

#lfc for comparison
hs_AvN.lfc <- read.csv('hs_kog_AvN_lfc_final.csv', header =TRUE)

hs_AvN.lth=kog.mwu(hs_AvN.lfc,hs_kog2gene) 
hs_AvN.lth 

#lfc (can switch them out to see the diffrences)
hs_HvN.lfc <- read.csv('hs_kog_HvN_lfc_final.csv', header =TRUE)
hs_HvN.lth=kog.mwu(hs_HvN.lfc,hs_kog2gene)
hs_HvN.lth

# hs HvA
hs_HvA.lfc <- read.csv('hs_kog_HvA_lfc_final.csv', header =TRUE)
hs_HvA.lth=kog.mwu(hs_HvA.lfc,hs_kog2gene)
hs_HvA.lth

#KOG mwu tethya
#long exposure to hypoxia 
t_kog2gene=read.table("tethya_all_gene2kog_unique.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)
#lfc for comparison
t_long.lfc <- read.csv('t_long_kog_lfc_final.csv', header =TRUE)
t_long.lth=kog.mwu(t_long.lfc,t_kog2gene) 
t_long.lth 

#short exposure to anoxia
#lfc for comparison
t_shock.lfc <- read.csv('t_shock_kog_lfc_final.csv', header =TRUE)
t_shock.lth=kog.mwu(t_shock.lfc,t_kog2gene) 
t_shock.lth 



# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("Es_Anox"=es_AvN.lth,"Es_Hypox"=es_HvN.lth,"Es_Hypox_vs_Anox"=es_HvA.lth,
                                "Hs_Anox"=hs_AvN.lth,"Hs_Hypox"=hs_HvN.lth,"Hs_Hypox_vs_Anox"=hs_HvA.lth,
                                "Tethya_Hypox_Long"=t_long.lth,"Tethya_Anox_short"=t_shock.lth))


#Add data from Stader et al 2016 Red fluorescence in coral larvae is associated with a diapause-like state
#dauer data
##updated rows in excel corresponding to the following kogs to make sure commas are in the right passes as follows:
#Replication, recombination and repair
#Intracellular trafficking, secretion, and vesicular transport
#Cell cycle control, cell division, chromosome partitioning
#Posttranslational modification, protein turnover, chaperones
#Secondary metabolites biosynthesis, transport and catabolism
strader <- read.csv('Strader_KOGtableDeltaRanks.csv', header =TRUE)
dauer <- data.frame(strader$X, strader$dauer)
names(dauer)[1] <- 'term'
names(dauer)[2] <- 'delta.rank'
#prepare to join ktable
dauer2 <- dauer %>% remove_rownames %>% column_to_rownames(var="term")
#name for heatmap, corrs, etc
names(dauer2)[1] <- 'Dauer'
ktable_final <- merge(ktable, dauer2, by=0)
#annoying final formatting
ktable_final2 <- ktable_final %>% remove_rownames %>% column_to_rownames(var="Row.names")

#diapuase
dia <- data.frame(strader$X, strader$diapause)
names(dia)[1] <- 'term'
names(dia)[2] <- 'delta.rank'
#prepare to join ktable
dia2 <- dia %>% remove_rownames %>% column_to_rownames(var="term")
#name for heatmap, corrs, etc
names(dia2)[1] <- 'Diapause'
#join to final table 
ktable_final3 <- merge(ktable_final2, dia2, by=0)
#annoying final formatting
ktable_final4 <- ktable_final3 %>% remove_rownames %>% column_to_rownames(var="Row.names")

#coral long term stress
coral <- data.frame(strader$X, strader$adults.lt)
names(coral)[1] <- 'term'
names(coral)[2] <- 'delta.rank'
#prepare to join ktable
coral2 <- coral %>% remove_rownames %>% column_to_rownames(var="term")
#name for heatmap, corrs, etc
names(coral2)[1] <- 'Coral_heat_tolerance'
#join to final table 
ktable_final5 <- merge(ktable_final4, coral2, by=0)
#annoying final formatting
ktable_final6 <- ktable_final5 %>% remove_rownames %>% column_to_rownames(var="Row.names")



# Making a heatmap with hierarchical clustering trees: 
pheatmap(as.matrix(ktable_final6), clustering_distance_cols="correlation") 

#turn off clustering to make figure more readable
#keep dendrogram 

pheatmap(as.matrix(ktable_final6), cluster_cols=F,clustering_distance_rows="correlation", clustering_distance_cols="correlation")

#make heatmap
#pheatmap(heat_sig_es_HvN_RRR_gn_final, color=col,cluster_cols=F,clustering_distance_rows="correlation")

# exploring correlations between datasets
pairs(ktable_final6, lower.panel = panel.smooth, upper.panel = panel.cor)
# p-values of these correlations in the upper panel:
pairs(ktable_final6, lower.panel = panel.smooth, upper.panel = panel.cor.pval)

# plotting individual delta-rank correlations where P<0.05:
#adjust for desired close up from the above


corrPlot(x="Es_Anox",y="Es_Hypox_vs_Anox",ktable_final6)
corrPlot(x="Es_Hypox",y="Hs_Hypox_vs_Anox",ktable_final6)
corrPlot(x="Hs_Anox",y="Hs_Hypox_vs_Anox",ktable_final6)
corrPlot(x="Hs_Hypox",y="Hs_Hypox_vs_Anox",ktable_final6)
corrPlot(x="Diapause",y="Es_Hypox",ktable_final6)
corrPlot(x="Hs_Hypox",y="Coral_heat_tolerance",ktable_final6)

#original correlations are fine with addition of tethya  
corrPlot(x="Es_Hypox_vs_Anox",y="Tethya_Anox_short",ktable_final6)
corrPlot(x="Es_Hypox_vs_Anox",y="Tethya_Hypox_Long",ktable_final6) #negative corr 
corrPlot(x="Es_Anox",y="Tethya_Hypox_Long",ktable_final6) #negative correlaton 

#400 pixels x 400

corrPlot(y="Diapause",x="Es_Anox",ktable_final6)
corrPlot(y="Dauer",x="Es_Anox",ktable_final6)
corrPlot(y="Diapause",x="Hs_Anox",ktable_final6)
corrPlot(y="Dauer",x="Hs_Anox",ktable_final6)

#############
#isolating genes from specific KOG and making a heatmap of those genes (only for sig (p<0.1))) 
###
#all should have an annotation because they have a KOG 
####################Now for some heatmaps - lets look at what genes annotated with interesting terms are 

#see https://support.bioconductor.org/p/133313/
vst <- vst(dds_es2, blind=FALSE)

vst <- assay(vst)

vst <- as.data.frame(vst)

#pull out list of only significant genes

#vst_sig <- vst[rownames(vst) %in% significant_gene_names,]
#heat <- t(scale(t(vst_sig)))

#sets rowmeans to zero - heatmap = zscore above or below zero 
heat <- t(scale(t(vst)))

col=color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")),bias=0.75)(30)
#lower bias number gives more blues; higher bias gives more reds
#mess with bias to get white at 0

pheatmap(head(heat),color=col,cluster_cols=F,clustering_distance_rows="correlation")

#get only significantly different genes
#for es HvN
df1 <-data.frame(res_es2_HvN)
table(df1$padj < 0.1)

df1 %>% filter(padj < 0.1) -> sig_es_HvN

length(rownames(sig_es_HvN))
#737 :)

heat_sig_es_HvN <- heat[rownames(heat) %in% rownames(sig_es_HvN),]

length(rownames(sig_es_HvN))
#737 :)

#Get genes with desired KOG 
es_kog2gene %>%
  filter(V2 == 'Replication, recombination and repair') -> es_RRR

#get genes of this kog from desired list of significant genes  
heat_sig_es_HvN_RRR <- heat_sig_es_HvN[rownames(heat_sig_es_HvN) %in% es_RRR$V1,]

length(rownames(es_RRR))
#618

length(rownames(heat_sig_es_HvN_RRR))
#16 - should be very short

#get annotations
gn_es=read.table("E_sp2_sponge_ref_transcriptome_iso2geneName_final.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)

#make isogroups into rownames 
gn_es2 <- gn_es %>% remove_rownames %>% column_to_rownames(var="V1")
#merge by isogroup
heat_sig_es_HvN_RRR_gn <- merge(heat_sig_es_HvN_RRR, gn_es2, by=0)

#get rid of spaces in gene names
#get rid of spaces at the begining of gene names since they cause problems later 
heat_sig_es_HvN_RRR_gn$V2 <- trimws(heat_sig_es_HvN_RRR_gn$V2)
heat_sig_es_HvN_RRR_gn$V2 <- gsub(" ", "_", heat_sig_es_HvN_RRR_gn$V2)

#replace all special characters? e.g.
#rownames(counts_es) <- gsub("[()]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[:]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[/]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[.]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[,]", "_", rownames(counts_es))

#make rownames the annotation rather than the isogroup
#first rows have to have individual names (two genes with the same annotation need to be made unique)
heat_sig_es_HvN_RRR_gn$V2 <- make.names(heat_sig_es_HvN_RRR_gn$V2, unique = TRUE)
#make gene names into rownames
heat_sig_es_HvN_RRR_gn %>% remove_rownames %>% column_to_rownames(var="V2") -> heat_sig_es_HvN_RRR_gn

#get rid of isogroup designation
heat_sig_es_HvN_RRR_gn <- heat_sig_es_HvN_RRR_gn[,-1] 

#reorder columns based on oxygen category (see Carly's script)
#norm, hyp, an
heat_sig_es_HvN_RRR_gn_final<-heat_sig_es_HvN_RRR_gn[,c(11,13,14,15,8,9,1,2,3,4,5,6,7,10,12)]

#make heatmap
pheatmap(heat_sig_es_HvN_RRR_gn_final, color=col,cluster_cols=F,clustering_distance_rows="correlation")

#600X350
#for es AvN
df1 <-data.frame(res_es2_AvN)
table(df1$padj < 0.1)

df1 %>% filter(padj < 0.1) -> sig_es_AvN

length(rownames(sig_es_AvN))
# :)

heat_sig_es_AvN <- heat[rownames(heat) %in% rownames(sig_es_AvN),]

length(rownames(sig_es_AvN))
#1376 :)

#Get genes with desired KOG 
es_kog2gene %>%
  filter(V2 == 'Replication, recombination and repair') -> es_RRR

#get genes of this kog from desired list of significant genes  
heat_sig_es_AvN_RRR <- heat_sig_es_AvN[rownames(heat_sig_es_AvN) %in% es_RRR$V1,]

length(rownames(es_RRR))
#618

length(rownames(heat_sig_es_AvN_RRR))
#7 - should be very short

#get annotations
gn_es=read.table("E_sp2_sponge_ref_transcriptome_iso2geneName_final.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)

#make isogroups into rownames 
gn_es2 <- gn_es %>% remove_rownames %>% column_to_rownames(var="V1")
#merge by isogroup
heat_sig_es_AvN_RRR_gn <- merge(heat_sig_es_AvN_RRR, gn_es2, by=0)

#get rid of spaces in gene names
#get rid of spaces at the begining of gene names since they cause problems later 
heat_sig_es_AvN_RRR_gn$V2 <- trimws(heat_sig_es_AvN_RRR_gn$V2)
heat_sig_es_AvN_RRR_gn$V2 <- gsub(" ", "_", heat_sig_es_AvN_RRR_gn$V2)

#replace all special characters? e.g.
#rownames(counts_es) <- gsub("[()]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[:]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[/]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[.]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[,]", "_", rownames(counts_es))

#make rownames the annotation rather than the isogroup
#first rows have to have individual names (two genes with the same annotation need to be made unique)
#heat_sig_es_AvN_RRR_gn$V2 <- make.names(heat_sig_es_AvN_RRR_gn$V2, unique = TRUE)
#make gene names into rownames
heat_sig_es_AvN_RRR_gn %>% remove_rownames %>% column_to_rownames(var="V2") -> heat_sig_es_AvN_RRR_gn

#get rid of isogroup designation
heat_sig_es_AvN_RRR_gn <- heat_sig_es_AvN_RRR_gn[,-1] 

#reorder columns based on oxygen category (see Carly's script)
#norm, hyp, an
heat_sig_es_AvN_RRR_gn_final<-heat_sig_es_AvN_RRR_gn[,c(11,13,14,15,8,9,1,2,3,4,5,6,7,10,12)]

#make heatmap
pheatmap(heat_sig_es_AvN_RRR_gn_final, color=col,cluster_cols=F,clustering_distance_rows="correlation")

#600X200
#for es AvH
df1 <-data.frame(res_es2_HvA)
table(df1$padj < 0.1)

df1 %>% filter(padj < 0.1) -> sig_es_AvH

length(rownames(sig_es_AvH))
#1090 :)

heat_sig_es_AvH <- heat[rownames(heat) %in% rownames(sig_es_AvH),]

length(rownames(sig_es_AvH))
#1090 :)

#Get genes with desired KOG 
es_kog2gene %>%
  filter(V2 == 'Replication, recombination and repair') -> es_RRR

#get genes of this kog from desired list of significant genes  
heat_sig_es_AvH_RRR <- heat_sig_es_AvH[rownames(heat_sig_es_AvH) %in% es_RRR$V1,]

length(rownames(es_RRR))
#618

length(rownames(heat_sig_es_AvH_RRR))
#11 - should be very short

#get annotations
gn_es=read.table("E_sp2_sponge_ref_transcriptome_iso2geneName_final.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)

#make isogroups into rownames 
gn_es2 <- gn_es %>% remove_rownames %>% column_to_rownames(var="V1")
#merge by isogroup
heat_sig_es_AvH_RRR_gn <- merge(heat_sig_es_AvH_RRR, gn_es2, by=0)

#get rid of spaces in gene names
#get rid of spaces at the begining of gene names since they cause problems later 
heat_sig_es_AvH_RRR_gn$V2 <- trimws(heat_sig_es_AvH_RRR_gn$V2)
heat_sig_es_AvH_RRR_gn$V2 <- gsub(" ", "_", heat_sig_es_AvH_RRR_gn$V2)

#replace all special characters? e.g.
#rownames(counts_es) <- gsub("[()]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[:]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[/]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[.]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[,]", "_", rownames(counts_es))

#make rownames the annotation rather than the isogroup
#first rows have to have individual names (two genes with the same annotation need to be made unique)
heat_sig_es_AvH_RRR_gn$V2 <- make.names(heat_sig_es_AvH_RRR_gn$V2, unique = TRUE)
#make gene names into rownames
heat_sig_es_AvH_RRR_gn %>% remove_rownames %>% column_to_rownames(var="V2") -> heat_sig_es_AvH_RRR_gn

#get rid of isogroup designation
heat_sig_es_AvH_RRR_gn <- heat_sig_es_AvH_RRR_gn[,-1] 

#reorder columns based on oxygen category (see Carly's script)
#norm, hyp, an
heat_sig_es_AvH_RRR_gn_final<-heat_sig_es_AvH_RRR_gn[,c(11,13,14,15,8,9,1,2,3,4,5,6,7,10,12)]

#make heatmap
pheatmap(heat_sig_es_AvH_RRR_gn_final, color=col,cluster_cols=F,clustering_distance_rows="correlation")

#600X350

###save csvs for easy update to heatmaps 
#RRR 
#es
write.csv(heat_sig_es_AvH_RRR_gn_final, file='heat_es_HvA_RRR.csv', row.names = TRUE)
write.csv(heat_sig_es_AvN_RRR_gn_final, file='heat_es_AvN_RRR.csv', row.names = TRUE)
write.csv(heat_sig_es_HvN_RRR_gn_final, file='heat_es_HvN_RRR.csv', row.names = TRUE)

###################
#HS- genes in RRR KOG 
###################
#see https://support.bioconductor.org/p/133313/
vst <- vst(dds, blind=FALSE)

vst <- assay(vst)

vst <- as.data.frame(vst)

#pull out list of only significant genes

#vst_sig <- vst[rownames(vst) %in% significant_gene_names,]
#heat <- t(scale(t(vst_sig)))

#sets rowmeans to zero - heatmap = zscore above or below zero 
heat <- t(scale(t(vst)))

col=color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")),bias=0.75)(30)
#lower bias number gives more blues; higher bias gives more reds
#mess with bias to get white at 0

pheatmap(head(heat),color=col,cluster_cols=F,clustering_distance_rows="correlation")

#get only significantly different genes
#for es HvN
df1 <-data.frame(res_HvN)
table(df1$padj < 0.1)

df1 %>% filter(padj < 0.1) -> sig_hs_HvN

length(rownames(sig_hs_HvN))
#523 :)

heat_sig_hs_HvN <- heat[rownames(heat) %in% rownames(sig_hs_HvN),]

length(rownames(sig_hs_HvN))
#523 :)

#Get genes with desired KOG 
hs_kog2gene %>%
  filter(V2 == 'Replication, recombination and repair') -> hs_RRR

#get genes of this kog from desired list of significant genes  
heat_sig_hs_HvN_RRR <- heat_sig_hs_HvN[rownames(heat_sig_hs_HvN) %in% hs_RRR$V1,]

length(rownames(hs_RRR))
#730

length(rownames(heat_sig_hs_HvN_RRR))
#16 - should be very short

#get annotations
gn_hs=read.table("H_stellifera_sponge_ref_transcriptome_iso2geneName_final.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)

#make isogroups into rownames 
gn_hs2 <- gn_hs %>% remove_rownames %>% column_to_rownames(var="V1")
#merge by isogroup
heat_sig_hs_HvN_RRR_gn <- merge(heat_sig_hs_HvN_RRR, gn_hs2, by=0)

#get rid of spaces in gene names
#get rid of spaces at the begining of gene names since they cause problems later 
heat_sig_hs_HvN_RRR_gn$V2 <- trimws(heat_sig_hs_HvN_RRR_gn$V2)
heat_sig_hs_HvN_RRR_gn$V2 <- gsub(" ", "_", heat_sig_hs_HvN_RRR_gn$V2)

#replace all special characters? e.g.
#rownames(counts_es) <- gsub("[()]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[:]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[/]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[.]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[,]", "_", rownames(counts_es))

#make rownames the annotation rather than the isogroup
#first rows have to have individual names (two genes with the same annotation need to be made unique)
#heat_sig_hs_HvN_RRR_gn$V2 <- make.names(heat_sig_hs_HvN_RRR_gn$V2, unique = TRUE)
#make gene names into rownames
heat_sig_hs_HvN_RRR_gn %>% remove_rownames %>% column_to_rownames(var="V2") -> heat_sig_hs_HvN_RRR_gn

#get rid of isogroup designation
heat_sig_hs_HvN_RRR_gn <- heat_sig_hs_HvN_RRR_gn[,-1] 

#reorder columns based on oxygen category (see Carly's script)
#norm, hyp, an
heat_sig_hs_HvN_RRR_gn_final<-heat_sig_hs_HvN_RRR_gn[,c(5,1,2,3,4,6,7,8)]

#make heatmap
pheatmap(heat_sig_hs_HvN_RRR_gn_final, color=col,cluster_cols=F,clustering_distance_rows="correlation")

#600X150
#for hs AvN
df1 <-data.frame(res_AvN)
table(df1$padj < 0.1)

df1 %>% filter(padj < 0.1) -> sig_hs_AvN

length(rownames(sig_hs_AvN))
#284 :)

heat_sig_hs_AvN <- heat[rownames(heat) %in% rownames(sig_hs_AvN),]

length(rownames(sig_hs_AvN))
#284 :)

#Get genes with desired KOG 
hs_kog2gene %>%
  filter(V2 == 'Replication, recombination and repair') -> hs_RRR

#get genes of this kog from desired list of significant genes  
heat_sig_hs_AvN_RRR <- heat_sig_hs_AvN[rownames(heat_sig_hs_AvN) %in% hs_RRR$V1,]

length(rownames(hs_RRR))
#730

length(rownames(heat_sig_hs_AvN_RRR))
#0 genes!!! - should be very short
#makes sense with such only one sample from normoxia

#get annotations
gn_hs=read.table("H_stellifera_sponge_ref_transcriptome_iso2geneName_final.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)

#make isogroups into rownames 
gn_hs2 <- gn_hs %>% remove_rownames %>% column_to_rownames(var="V1")
#merge by isogroup
heat_sig_hs_AvN_RRR_gn <- merge(heat_sig_hs_AvN_RRR, gn_hs2, by=0)

#get rid of spaces in gene names
#get rid of spaces at the begining of gene names since they cause problems later 
heat_sig_hs_AvN_RRR_gn$V2 <- trimws(heat_sig_hs_AvN_RRR_gn$V2)
heat_sig_hs_AvN_RRR_gn$V2 <- gsub(" ", "_", heat_sig_hs_AvN_RRR_gn$V2)

#replace all special characters? e.g.
#rownames(counts_es) <- gsub("[()]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[:]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[/]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[.]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[,]", "_", rownames(counts_es))

#make rownames the annotation rather than the isogroup
#first rows have to have individual names (two genes with the same annotation need to be made unique)
#heat_sig_es_AvN_RRR_gn$V2 <- make.names(heat_sig_es_AvN_RRR_gn$V2, unique = TRUE)
#make gene names into rownames
heat_sig_hs_AvN_RRR_gn %>% remove_rownames %>% column_to_rownames(var="V2") -> heat_sig_hs_AvN_RRR_gn

#get rid of isogroup designation
heat_sig_hs_AvN_RRR_gn <- heat_sig_hs_AvN_RRR_gn[,-1] 

#reorder columns based on oxygen category (see Carly's script)
#norm, hyp, an
heat_sig_hs_AvN_RRR_gn_final<-heat_sig_hs_AvN_RRR_gn[,c(5,1,2,3,4,6,7,8)]

#make heatmap
pheatmap(heat_sig_hs_AvN_RRR_gn_final, color=col,cluster_cols=F,clustering_distance_rows="correlation")

#600X200
#for hs AvH
df1 <-data.frame(res_HvA)
table(df1$padj < 0.1)

df1 %>% filter(padj < 0.1) -> sig_hs_AvH

length(rownames(sig_hs_AvH))
#2154 :)

heat_sig_hs_AvH <- heat[rownames(heat) %in% rownames(sig_hs_AvH),]

length(rownames(sig_hs_AvH))
#2154 :)

#Get genes with desired KOG 
hs_kog2gene %>%
  filter(V2 == 'Replication, recombination and repair') -> hs_RRR

#get genes of this kog from desired list of significant genes  
heat_sig_hs_AvH_RRR <- heat_sig_hs_AvH[rownames(heat_sig_hs_AvH) %in% hs_RRR$V1,]

length(rownames(hs_RRR))
#730

length(rownames(heat_sig_hs_AvH_RRR))
#19 - should be very short

#get annotations
gn_hs=read.table("H_stellifera_sponge_ref_transcriptome_iso2geneName_final.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)

#make isogroups into rownames 
gn_hs2 <- gn_hs %>% remove_rownames %>% column_to_rownames(var="V1")
#merge by isogroup
heat_sig_hs_AvH_RRR_gn <- merge(heat_sig_hs_AvH_RRR, gn_hs2, by=0)

#get rid of spaces in gene names
#get rid of spaces at the begining of gene names since they cause problems later 
heat_sig_hs_AvH_RRR_gn$V2 <- trimws(heat_sig_hs_AvH_RRR_gn$V2)
heat_sig_hs_AvH_RRR_gn$V2 <- gsub(" ", "_", heat_sig_hs_AvH_RRR_gn$V2)

#replace all special characters? e.g.
#rownames(counts_es) <- gsub("[()]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[:]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[/]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[.]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[,]", "_", rownames(counts_es))

#make rownames the annotation rather than the isogroup
#first rows have to have individual names (two genes with the same annotation need to be made unique)
heat_sig_hs_AvH_RRR_gn$V2 <- make.names(heat_sig_hs_AvH_RRR_gn$V2, unique = TRUE)
#make gene names into rownames
heat_sig_hs_AvH_RRR_gn %>% remove_rownames %>% column_to_rownames(var="V2") -> heat_sig_hs_AvH_RRR_gn

#get rid of isogroup designation
heat_sig_hs_AvH_RRR_gn <- heat_sig_hs_AvH_RRR_gn[,-1] 

#reorder columns based on oxygen category (see Carly's script)
#norm, hyp, an
heat_sig_hs_AvH_RRR_gn_final<-heat_sig_hs_AvH_RRR_gn[,c(5,1,2,3,4,6,7,8)]

#make heatmap
pheatmap(heat_sig_hs_AvH_RRR_gn_final, border_color = "grey60", color=col,cluster_cols=F,clustering_distance_rows="correlation")
#wierdly differenot format for no reason
#600X350

###save csvs for easy update to heatmaps 
#RRR 
#hs
write.csv(heat_sig_hs_AvH_RRR_gn_final, file='heat_hs_HvA_RRR.csv', row.names = TRUE)
#write.csv(heat_sig_hs_AvN_RRR_gn_final, file='heat_hs_AvN_RRR.csv', row.names = TRUE)
write.csv(heat_sig_hs_HvN_RRR_gn_final, file='heat_hs_HvN_RRR.csv', row.names = TRUE)


#other two Kogs for es and the hs

########
#EPC KOG
########
vst <- vst(dds_es2, blind=FALSE)

vst <- assay(vst)

vst <- as.data.frame(vst)

#pull out list of only significant genes

#vst_sig <- vst[rownames(vst) %in% significant_gene_names,]
#heat <- t(scale(t(vst_sig)))

#sets rowmeans to zero - heatmap = zscore above or below zero 
heat <- t(scale(t(vst)))

col=color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")),bias=0.75)(30)
#lower bias number gives more blues; higher bias gives more reds
#mess with bias to get white at 0

pheatmap(head(heat),color=col,cluster_cols=F,clustering_distance_rows="correlation")

#get only significantly different genes
#for es HvN
df1 <-data.frame(res_es2_HvN)
table(df1$padj < 0.1)

df1 %>% filter(padj < 0.1) -> sig_es_HvN

length(rownames(sig_es_HvN))
#737 :)

heat_sig_es_HvN <- heat[rownames(heat) %in% rownames(sig_es_HvN),]

length(rownames(sig_es_HvN))
#737 :)

#Get genes with desired KOG 
es_kog2gene %>%
  filter(V2 == 'Energy production and conversion') -> es_EEE

#get genes of this kog from desired list of significant genes  
heat_sig_es_HvN_EEE <- heat_sig_es_HvN[rownames(heat_sig_es_HvN) %in% es_EEE$V1,]

length(rownames(es_EEE))
#712

length(rownames(heat_sig_es_HvN_EEE))
#16 - should be very short

#get annotations
gn_es=read.table("E_sp2_sponge_ref_transcriptome_iso2geneName_final.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)

#make isogroups into rownames 
gn_es2 <- gn_es %>% remove_rownames %>% column_to_rownames(var="V1")
#merge by isogroup
heat_sig_es_HvN_EEE_gn <- merge(heat_sig_es_HvN_EEE, gn_es2, by=0)

#get rid of spaces in gene names
#get rid of spaces at the begining of gene names since they cause problems later 
heat_sig_es_HvN_EEE_gn$V2 <- trimws(heat_sig_es_HvN_EEE_gn$V2)
heat_sig_es_HvN_EEE_gn$V2 <- gsub(" ", "_", heat_sig_es_HvN_EEE_gn$V2)

#make rownames the annotation rather than the isogroup
#first rows have to have individual names (two genes with the same annotation need to be made unique)
heat_sig_es_HvN_EEE_gn$V2 <- make.names(heat_sig_es_HvN_EEE_gn$V2, unique = TRUE)
#make gene names into rownames
heat_sig_es_HvN_EEE_gn %>% remove_rownames %>% column_to_rownames(var="V2") -> heat_sig_es_HvN_EEE_gn

#get rid of isogroup designation
heat_sig_es_HvN_EEE_gn <- heat_sig_es_HvN_EEE_gn[,-1] 

#reorder columns based on oxygen category (see Carly's script)
#norm, hyp, an
heat_sig_es_HvN_EEE_gn_final<-heat_sig_es_HvN_EEE_gn[,c(11,13,14,15,8,9,1,2,3,4,5,6,7,10,12)]

#shorten super long gene names 
rownames(heat_sig_es_HvN_EEE_gn_final) <- strtrim(rownames(heat_sig_es_HvN_EEE_gn_final), 75)

#make heatmap
pheatmap(heat_sig_es_HvN_EEE_gn_final, color=col,cluster_cols=F,clustering_distance_rows="correlation")

#600X350
#for es AvN
df1 <-data.frame(res_es2_AvN)
table(df1$padj < 0.1)

df1 %>% filter(padj < 0.1) -> sig_es_AvN

length(rownames(sig_es_AvN))
# :)

heat_sig_es_AvN <- heat[rownames(heat) %in% rownames(sig_es_AvN),]

length(rownames(sig_es_AvN))
#1376 :)

#Get genes with desired KOG 
es_kog2gene %>%
  filter(V2 == 'Energy production and conversion') -> es_EEE

#get genes of this kog from desired list of significant genes  
heat_sig_es_AvN_EEE <- heat_sig_es_AvN[rownames(heat_sig_es_AvN) %in% es_EEE$V1,]

length(rownames(es_EEE))
#712

length(rownames(heat_sig_es_AvN_EEE))
#19 - should be very short

#get annotations
gn_es=read.table("E_sp2_sponge_ref_transcriptome_iso2geneName_final.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)

#make isogroups into rownames 
gn_es2 <- gn_es %>% remove_rownames %>% column_to_rownames(var="V1")
#merge by isogroup
heat_sig_es_AvN_EEE_gn <- merge(heat_sig_es_AvN_EEE, gn_es2, by=0)

#get rid of spaces in gene names
#get rid of spaces at the begining of gene names since they cause problems later 
heat_sig_es_AvN_EEE_gn$V2 <- trimws(heat_sig_es_AvN_EEE_gn$V2)
heat_sig_es_AvN_EEE_gn$V2 <- gsub(" ", "_", heat_sig_es_AvN_EEE_gn$V2)

#replace all special characters? e.g.
#rownames(counts_es) <- gsub("[()]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[:]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[/]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[.]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[,]", "_", rownames(counts_es))

#make rownames the annotation rather than the isogroup
#first rows have to have individual names (two genes with the same annotation need to be made unique)
heat_sig_es_AvN_EEE_gn$V2 <- make.names(heat_sig_es_AvN_EEE_gn$V2, unique = TRUE)
#make gene names into rownames
heat_sig_es_AvN_EEE_gn %>% remove_rownames %>% column_to_rownames(var="V2") -> heat_sig_es_AvN_EEE_gn

#get rid of isogroup designation
heat_sig_es_AvN_EEE_gn <- heat_sig_es_AvN_EEE_gn[,-1] 

#reorder columns based on oxygen category (see Carly's script)
#norm, hyp, an
heat_sig_es_AvN_EEE_gn_final<-heat_sig_es_AvN_EEE_gn[,c(11,13,14,15,8,9,1,2,3,4,5,6,7,10,12)]

#make heatmap
pheatmap(heat_sig_es_AvN_EEE_gn_final, color=col,cluster_cols=F,clustering_distance_rows="correlation")

#600X200
#for es AvH
df1 <-data.frame(res_es2_HvA)
table(df1$padj < 0.1)

df1 %>% filter(padj < 0.1) -> sig_es_AvH

length(rownames(sig_es_AvH))
#1090 :)

heat_sig_es_AvH <- heat[rownames(heat) %in% rownames(sig_es_AvH),]

length(rownames(sig_es_AvH))
#1090 :)

#Get genes with desired KOG 
es_kog2gene %>%
  filter(V2 == 'Energy production and conversion') -> es_EEE

#get genes of this kog from desired list of significant genes  
heat_sig_es_AvH_EEE <- heat_sig_es_AvH[rownames(heat_sig_es_AvH) %in% es_EEE$V1,]

length(rownames(es_EEE))
#712

length(rownames(heat_sig_es_AvH_EEE))
#16 - should be very short

#get annotations
gn_es=read.table("E_sp2_sponge_ref_transcriptome_iso2geneName_final.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)

#make isogroups into rownames 
gn_es2 <- gn_es %>% remove_rownames %>% column_to_rownames(var="V1")
#merge by isogroup
heat_sig_es_AvH_EEE_gn <- merge(heat_sig_es_AvH_EEE, gn_es2, by=0)

#get rid of spaces in gene names
#get rid of spaces at the begining of gene names since they cause problems later 
heat_sig_es_AvH_EEE_gn$V2 <- trimws(heat_sig_es_AvH_EEE_gn$V2)
heat_sig_es_AvH_EEE_gn$V2 <- gsub(" ", "_", heat_sig_es_AvH_EEE_gn$V2)

#replace all special characters? e.g.
#rownames(counts_es) <- gsub("[()]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[:]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[/]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[.]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[,]", "_", rownames(counts_es))

#make rownames the annotation rather than the isogroup
#first rows have to have individual names (two genes with the same annotation need to be made unique)
heat_sig_es_AvH_EEE_gn$V2 <- make.names(heat_sig_es_AvH_EEE_gn$V2, unique = TRUE)
#make gene names into rownames
heat_sig_es_AvH_EEE_gn %>% remove_rownames %>% column_to_rownames(var="V2") -> heat_sig_es_AvH_EEE_gn

#get rid of isogroup designation
heat_sig_es_AvH_EEE_gn <- heat_sig_es_AvH_EEE_gn[,-1] 

#reorder columns based on oxygen category (see Carly's script)
#norm, hyp, an
heat_sig_es_AvH_EEE_gn_final<-heat_sig_es_AvH_EEE_gn[,c(11,13,14,15,8,9,1,2,3,4,5,6,7,10,12)]

#trim names
rownames(heat_sig_es_AvH_EEE_gn_final) <- strtrim(rownames(heat_sig_es_AvH_EEE_gn_final), 75)

#make heatmap
pheat_sig_es_AvH_EEE <- pheatmap(heat_sig_es_AvH_EEE_gn_final, color=col,cluster_cols=F,clustering_distance_rows="correlation")


#can I combine the heatmaps?
#no, pheatmap must not be ggplot based

#600X350

###save csvs for easy update to heatmaps 
#EEE 
#es
write.csv(heat_sig_es_AvH_EEE_gn_final, file='heat_es_HvA_EEE.csv', row.names = TRUE)
write.csv(heat_sig_es_AvN_EEE_gn_final, file='heat_es_AvN_EEE.csv', row.names = TRUE)
write.csv(heat_sig_es_HvN_EEE_gn_final, file='heat_es_HvN_EEE.csv', row.names = TRUE)

###################
#HS- genes in EEE KOG 
###################
#see https://support.bioconductor.org/p/133313/
vst <- vst(dds, blind=FALSE)

vst <- assay(vst)

vst <- as.data.frame(vst)

#pull out list of only significant genes

#vst_sig <- vst[rownames(vst) %in% significant_gene_names,]
#heat <- t(scale(t(vst_sig)))

#sets rowmeans to zero - heatmap = zscore above or below zero 
heat <- t(scale(t(vst)))

col=color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")),bias=0.75)(30)
#lower bias number gives more blues; higher bias gives more reds
#mess with bias to get white at 0

pheatmap(head(heat),color=col,cluster_cols=F,clustering_distance_rows="correlation")

#get only significantly different genes
#for es HvN
df1 <-data.frame(res_HvN)
table(df1$padj < 0.1)

df1 %>% filter(padj < 0.1) -> sig_hs_HvN

length(rownames(sig_hs_HvN))
#523 :)

heat_sig_hs_HvN <- heat[rownames(heat) %in% rownames(sig_hs_HvN),]

length(rownames(sig_hs_HvN))
#523 :)

#Get genes with desired KOG 
hs_kog2gene %>%
  filter(V2 == 'Energy production and conversion') -> hs_EEE

#get genes of this kog from desired list of significant genes  
heat_sig_hs_HvN_EEE <- heat_sig_hs_HvN[rownames(heat_sig_hs_HvN) %in% hs_EEE$V1,]

length(rownames(hs_EEE))
#803

length(rownames(heat_sig_hs_HvN_EEE))
#10 - should be very short

#get annotations
gn_hs=read.table("H_stellifera_sponge_ref_transcriptome_iso2geneName_final.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)

#make isogroups into rownames 
gn_hs2 <- gn_hs %>% remove_rownames %>% column_to_rownames(var="V1")
#merge by isogroup
heat_sig_hs_HvN_EEE_gn <- merge(heat_sig_hs_HvN_EEE, gn_hs2, by=0)

#get rid of spaces in gene names
#get rid of spaces at the begining of gene names since they cause problems later 
heat_sig_hs_HvN_EEE_gn$V2 <- trimws(heat_sig_hs_HvN_EEE_gn$V2)
heat_sig_hs_HvN_EEE_gn$V2 <- gsub(" ", "_", heat_sig_hs_HvN_EEE_gn$V2)

#replace all special characters? e.g.
#rownames(counts_es) <- gsub("[()]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[:]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[/]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[.]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[,]", "_", rownames(counts_es))

#make rownames the annotation rather than the isogroup
#first rows have to have individual names (two genes with the same annotation need to be made unique)
heat_sig_hs_HvN_EEE_gn$V2 <- make.names(heat_sig_hs_HvN_EEE_gn$V2, unique = TRUE)
#make gene names into rownames
heat_sig_hs_HvN_EEE_gn %>% remove_rownames %>% column_to_rownames(var="V2") -> heat_sig_hs_HvN_EEE_gn

#get rid of isogroup designation
heat_sig_hs_HvN_EEE_gn <- heat_sig_hs_HvN_EEE_gn[,-1] 

#reorder columns based on oxygen category (see Carly's script)
#norm, hyp, an
heat_sig_hs_HvN_EEE_gn_final<-heat_sig_hs_HvN_EEE_gn[,c(5,1,2,3,4,6,7,8)]

#shorten gene names
rownames(heat_sig_hs_HvN_EEE_gn_final) <- strtrim(rownames(heat_sig_hs_HvN_EEE_gn_final), 75)

#make heatmap
pheatmap(heat_sig_hs_HvN_EEE_gn_final, color=col,cluster_cols=F,clustering_distance_rows="correlation")

#700x250
#for hs AvN
df1 <-data.frame(res_AvN)
table(df1$padj < 0.1)

df1 %>% filter(padj < 0.1) -> sig_hs_AvN

length(rownames(sig_hs_AvN))
#284 :)

heat_sig_hs_AvN <- heat[rownames(heat) %in% rownames(sig_hs_AvN),]

length(rownames(sig_hs_AvN))
#284 :)

#Get genes with desired KOG 
hs_kog2gene %>%
  filter(V2 == 'Energy production and conversion') -> hs_EEE

#get genes of this kog from desired list of significant genes  
heat_sig_hs_AvN_EEE <- heat_sig_hs_AvN[rownames(heat_sig_hs_AvN) %in% hs_EEE$V1,]

length(rownames(hs_EEE))
#803

length(rownames(heat_sig_hs_AvN_EEE))
#4 genes!!! - should be very short
#makes sense with such only one sample from normoxia

#get annotations
gn_hs=read.table("H_stellifera_sponge_ref_transcriptome_iso2geneName_final.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)

#make isogroups into rownames 
gn_hs2 <- gn_hs %>% remove_rownames %>% column_to_rownames(var="V1")
#merge by isogroup
heat_sig_hs_AvN_EEE_gn <- merge(heat_sig_hs_AvN_EEE, gn_hs2, by=0)

#get rid of spaces in gene names
#get rid of spaces at the begining of gene names since they cause problems later 
heat_sig_hs_AvN_EEE_gn$V2 <- trimws(heat_sig_hs_AvN_EEE_gn$V2)
heat_sig_hs_AvN_EEE_gn$V2 <- gsub(" ", "_", heat_sig_hs_AvN_EEE_gn$V2)

#replace all special characters? e.g.
#rownames(counts_es) <- gsub("[()]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[:]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[/]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[.]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[,]", "_", rownames(counts_es))

#make rownames the annotation rather than the isogroup
#first rows have to have individual names (two genes with the same annotation need to be made unique)
#heat_sig_es_AvN_EEE_gn$V2 <- make.names(heat_sig_es_AvN_EEE_gn$V2, unique = TRUE)
#make gene names into rownames
heat_sig_hs_AvN_EEE_gn %>% remove_rownames %>% column_to_rownames(var="V2") -> heat_sig_hs_AvN_EEE_gn

#get rid of isogroup designation
heat_sig_hs_AvN_EEE_gn <- heat_sig_hs_AvN_EEE_gn[,-1] 

#reorder columns based on oxygen category (see Carly's script)
#norm, hyp, an
heat_sig_hs_AvN_EEE_gn_final<-heat_sig_hs_AvN_EEE_gn[,c(5,1,2,3,4,6,7,8)]
rownames(heat_sig_hs_AvN_EEE_gn_final) <- strtrim(rownames(heat_sig_hs_AvN_EEE_gn_final), 75)

#make heatmap
pheatmap(heat_sig_hs_AvN_EEE_gn_final, color=col,cluster_cols=F,clustering_distance_rows="correlation")

#600X200
#for hs AvH
df1 <-data.frame(res_HvA)
table(df1$padj < 0.1)

df1 %>% filter(padj < 0.1) -> sig_hs_AvH

length(rownames(sig_hs_AvH))
#2154 :)

heat_sig_hs_AvH <- heat[rownames(heat) %in% rownames(sig_hs_AvH),]

length(rownames(sig_hs_AvH))
#2154 :)

#Get genes with desired KOG 
hs_kog2gene %>%
  filter(V2 == 'Energy production and conversion') -> hs_EEE

#get genes of this kog from desired list of significant genes  
heat_sig_hs_AvH_EEE <- heat_sig_hs_AvH[rownames(heat_sig_hs_AvH) %in% hs_EEE$V1,]

length(rownames(hs_EEE))
#803

length(rownames(heat_sig_hs_AvH_EEE))
#13 - should be very short

#get annotations
gn_hs=read.table("H_stellifera_sponge_ref_transcriptome_iso2geneName_final.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)

#make isogroups into rownames 
gn_hs2 <- gn_hs %>% remove_rownames %>% column_to_rownames(var="V1")
#merge by isogroup
heat_sig_hs_AvH_EEE_gn <- merge(heat_sig_hs_AvH_EEE, gn_hs2, by=0)

#get rid of spaces in gene names
#get rid of spaces at the begining of gene names since they cause problems later 
heat_sig_hs_AvH_EEE_gn$V2 <- trimws(heat_sig_hs_AvH_EEE_gn$V2)
heat_sig_hs_AvH_EEE_gn$V2 <- gsub(" ", "_", heat_sig_hs_AvH_EEE_gn$V2)

#replace all special characters? e.g.
#rownames(counts_es) <- gsub("[()]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[:]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[/]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[.]", "_", rownames(counts_es))
#rownames(counts_es) <- gsub("[,]", "_", rownames(counts_es))

#make rownames the annotation rather than the isogroup
#first rows have to have individual names (two genes with the same annotation need to be made unique)
heat_sig_hs_AvH_EEE_gn$V2 <- make.names(heat_sig_hs_AvH_EEE_gn$V2, unique = TRUE)
#make gene names into rownames
heat_sig_hs_AvH_EEE_gn %>% remove_rownames %>% column_to_rownames(var="V2") -> heat_sig_hs_AvH_EEE_gn

#get rid of isogroup designation
heat_sig_hs_AvH_EEE_gn <- heat_sig_hs_AvH_EEE_gn[,-1] 

#reorder columns based on oxygen category (see Carly's script)
#norm, hyp, an
heat_sig_hs_AvH_EEE_gn_final<-heat_sig_hs_AvH_EEE_gn[,c(5,1,2,3,4,6,7,8)]

rownames(heat_sig_hs_AvH_EEE_gn_final) <- strtrim(rownames(heat_sig_hs_AvH_EEE_gn_final), 75)
#make heatmap
pheatmap(heat_sig_hs_AvH_EEE_gn_final, border_color = "grey60", color=col,cluster_cols=F,clustering_distance_rows="correlation")
#wierdly differenot format for no reason
#600X350

###save csvs for easy update to heatmaps 
#EEE 
#hs
write.csv(heat_sig_hs_AvH_EEE_gn_final, file='heat_hs_HvA_EEE.csv', row.names = TRUE)
#write.csv(heat_sig_hs_AvN_EEE_gn_final, file='heat_hs_AvN_EEE.csv', row.names = TRUE)
write.csv(heat_sig_hs_HvN_EEE_gn_final, file='heat_hs_HvN_EEE.csv', row.names = TRUE)

#########
#Prokaryotic data
#########

###########
#Es - thaumarchaeota
###########
counts_esProk=read.table("allcounts_es_bin3_nitroso_only.txt",header=TRUE,row.names=1) #Reading in the table of counts per isogroup by sample - you must open R in the same folder or change the working directory

counts_es=counts_esProk[,order(names(counts_esProk))]
names(counts_es)

##conditions (metadata)
conditions=read.table("Conditions_es.txt", header = T)

conds_es <- conditions[ order(row.names(conditions)), ]

####load data into deseq table format
#includes experimental design' since it is from the quick start method
DESeq2table_es_prok <- DESeqDataSetFromMatrix(countData = counts_es,
                                              colData = conds_es,
                                              design= ~ Oxy_cat)

#make normoxic the 'reference' level 
DESeq2table_es_prok$Oxy_cat <- relevel(DESeq2table_es_prok$Oxy_cat, "Normoxic")

#runs differential expression anaylsis quickly
dds_es_p_a <- DESeq(DESeq2table_es_prok)

resultsNames(dds_es_p_a) # lists the coefficients
res_es_AvN_p_a <- results(dds_es_p_a, name="Oxy_cat_Anoxic_vs_Normoxic")
# or to shrink log fold changes association with condition:

#shrunken results beased on LFC
res_es_AvN_apeglm_p_a <- lfcShrink(dds_es_p_a, coef="Oxy_cat_Anoxic_vs_Normoxic", type="apeglm")

#significantly expressed genes based on full results
table(res_es_AvN_p_a$padj < 0.1)

#FALSE  TRUE 
# 1540     3 
#3 differentially expressed genes between normoxic and anoxic samples :)

#significantly expressed genes based on LFC results
table(res_es_AvN_apeglm_p_a$padj < 0.1)

##3 differentially expressed genes between normoxic and anoxic samples = same as above :)

#standard method - hypoxic vs normoxic
res_es_HvN_p_a <- results(dds_es_p_a, name="Oxy_cat_Hypoxic_vs_Normoxic")

table(res_es_HvN_p_a$padj < 0.1)
#FALSE  TRUE 
#1523    20 
#20 differentially expressed genes between hypoxic and normoxic

# or to shrink log fold changes association with condition:
res_es_HvN_apeglm_p_a <- lfcShrink(dds_es_p_a, coef="Oxy_cat_Hypoxic_vs_Normoxic", type="apeglm")

table(res_es_HvN_apeglm_p_a$padj < 0.1)
#same as above

length(res_es_AvN_p_a[,1])

#1543 genes in total for es thaums

#comparison (contrast), between hypoxic and anoxic samples
res_es_HvA_p_a <- results(dds_es_p_a, contrast = list('Oxy_cat_Hypoxic_vs_Normoxic', 'Oxy_cat_Anoxic_vs_Normoxic'))

table(res_es_HvA_p_a$padj < 0.1)
#FALSE  TRUE 
# 1543

#0 genes differentially expressed between hypoxic and anoxic

#shrunken results 
#apeglm does not have contrast ability
#need to try 'normal shrinkage, or ashr

res_es_HvA_ashr_p_a <- lfcShrink(dds_es_p_a, contrast = list("Oxy_cat_Anoxic_vs_Normoxic", "Oxy_cat_Hypoxic_vs_Normoxic"), type="ashr")

table(res_es_HvA_ashr_p_a$padj < 0.1)
#same number as above :)

#prep for PCA/heatmaps
#transformation 
vsd <- vst(dds_es_p_a, blind=FALSE)

df <- as.data.frame(colData(dds_es_p_a)[,c("Oxy_cat")])

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

#not all anoxic or hypoxic samples group together. One anoxic sample and one hypoxic sample are very close to eachother 

#pca plot 
plotPCA(vsd, intgroup=c('Oxy_cat'))

#added ntop, flag, but one anoxic sample always grouped with the normoxic ones.

#lots of separation from one anoxic sample
#some separation with anoxia and normoxia..
#maybe would be better to use continuous oxygen rather than categorical

#could be effectively the same if top left outlier is removed from the PCA plot 

#MA plots

plotMA(res_es_AvN_p_a)

plotMA(res_es_AvN_apeglm_p_a)

plotMA(res_es_HvN_apeglm_p_a)

plotMA(res_es_HvA_ashr_p_a)

##After calling plotMA, one can use the function identify to interactively detect the row number of individual genes by clicking on the plot. One can then recover the gene identifiers by saving the resulting indices
#idx <- identify(res_es_AvN_apeglm$baseMean, res_es_AvN_apeglm$log2FoldChange)
##print gene names from points you clicked on 
#rownames(res_es)[idx]

####possible removing of outliers 
boxplot(log10(assays(dds_es_p_a)[["cooks"]]), range=0, las=2)

#24 and 56, are both outliers, but they are all the anoxic samples, so keep


#plotting dispersion estimates
plotDispEsts(dds_es_p_a)


#PCA to remove outlier
PCA <- plotPCA(vsd, intgroup=c('Oxy_cat'), returnData=TRUE)
PCA

#DC108 and 110 were outliers based on the PCA, but will keep

#########
#summary of all degs for prok using package ReportingTools
#########

#add annotations
#work around for getting longer gene names (not just isogroup)
gn=read.table("es_bin3_nitroso_contig2geneName.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)

#gene_id
#gene_name

names(gn)[names(gn) == "V1"] <- "gene_id"

names(gn)[names(gn) == "V2"] <- "gene_name"

head(gn)



#add annotation to dds_es object 
#check if they match in order 
#all(rownames(dds_es_p) == gn$gene_id)
#false, they do not match (expected)

#have to match them
#gn2 <- gn[match(rownames(dds_es2), gn$gene_id),]
#doesn'#t# work because many genes are not annotated

#try with left join? 
#need to get a count table with a column name for the row data to count 
counts_esProk2 <- cbind(rownames(counts_es), counts_es)
rownames(counts_esProk2) <- NULL
colnames(counts_esProk2) <- c("gene_id", 'DC108', 'DC110', 'DC112', 'DC114', 'DC118', 'DC123', 'DC131', 'DC135', 'DC136', 'DC24', 'DC33', 'DC39', 'DC56', 'DC85', 'DC93', 'DC97')

joined <- left_join(counts_esProk2, gn, by = "gene_id")
joined$full <-paste(joined$gene_id,joined$gene_name)

#reorder so that combined column is first (also drop extra columns not needed for count matrix)
#also got rid of names column
joined2 <-joined[, c(19, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17)]

#make first column back into rowname
joined3 <- joined2[,-1]
rownames(joined3) <- joined2[,1]

#repeat deseq anaylsis with annotated count table 
counts_es=joined3[,order(names(joined3))]
names(counts_es)

#replace all spaces in row names (gene names) with underscore
rownames(counts_es) <- gsub(" ", "_", rownames(counts_es))
#replace all special characters
rownames(counts_es) <- gsub("[()]", "_", rownames(counts_es))
rownames(counts_es) <- gsub("[:]", "_", rownames(counts_es))
rownames(counts_es) <- gsub("[/]", "_", rownames(counts_es))
rownames(counts_es) <- gsub("[.]", "_", rownames(counts_es))
rownames(counts_es) <- gsub("[,]", "_", rownames(counts_es))
rownames(counts_es) <- gsub("[[]", "_", rownames(counts_es))
rownames(counts_es) <- gsub("[]]", "_", rownames(counts_es))
#limit rownames to 50 characters
rownames(counts_es) <- strtrim(rownames(counts_es), 50)


DESeq2table_es_names_p <- DESeqDataSetFromMatrix(countData = counts_es,
                                                 colData = conds_es,
                                                 design= ~ Oxy_cat)


DESeq2table_es_names_p$Oxy_cat <- relevel(DESeq2table_es_names_p$Oxy_cat, "Normoxic")

dds_es_p2 <- DESeq(DESeq2table_es_names_p)

#with releveling, of course it doesn't fit
resultsNames(dds_es_p2)

#hypoxia vs. anoxia
res_es_HvA_p2 <- results(dds_es_p2, contrast = list('Oxy_cat_Anoxic_vs_Normoxic', 'Oxy_cat_Hypoxic_vs_Normoxic'))

table(res_es_HvA_p2$padj < 0.1)

#same as above once releveled with normoxia as baseline.


#res_es3_HvA_ashr <- lfcShrink(dds_es3, contrast = list("Oxy_cat_Anoxic_vs_Normoxic", "Oxy_cat_Hypoxic_vs_Normoxic"), type="ashr")

#table(res_es3_HvA_ashr$padj < 0.1)
#same as above :)

#Normoxia vs anoxia
res_es_AvN_p2 <- results(dds_es_p2, name="Oxy_cat_Anoxic_vs_Normoxic")

table(res_es_AvN_p2$padj < 0.1)
# 1540     3 
#3  differentially expressed genes between anoxic and normoxic samples :)

#shrunken results beased on LFC
#res_es3_AvN_apeglm <- lfcShrink(dds_es3, coef="Oxy_cat_Anoxic_vs_Normoxic", type="apeglm")

#significantly expressed genes based on full results
#table(res_es3_AvN_apeglm$padj < 0.1)


#significantly expressed genes based on LFC results
#table(res_es3_AvN_apeglm$padj < 0.1)
##need to run the other multiple comparisons
#standard method - hypoxic vs normoxic
res_es_HvN_p2 <- results(dds_es_p2, name="Oxy_cat_Hypoxic_vs_Normoxic")

table(res_es_HvN_p2$padj < 0.1)
#FALSE  TRUE 
#50 differentially expressed genes between hypoxic and normoxic

# or to shrink log fold changes association with condition:
#res_es3_HvN_apeglm <- lfcShrink(dds_es3, coef="Oxy_cat_Hypoxic_vs_Normoxic", type="apeglm")
#table(res_es3_HvN_apeglm$padj < 0.1)

##amo genes in es
###amo a,b and c genes?
#k141_1145589_26	Archaeal ammonia monooxygenase subunit A (AmoA)
#k141_1145589_27	-
#  k141_1145589_28	PFAM Ammonia monooxygenase methane monooxygenase, subunit C
#k141_1145589_29	Monooxygenase subunit B protein

#heatmap?
vst <- varianceStabilizingTransformation(dds_es_p2, blind=FALSE)

vst <- assay(vst)

vst <- as.data.frame(vst)

#pull out list of only significant genes

#vst_sig <- vst[rownames(vst) %in% significant_gene_names,]
#heat <- t(scale(t(vst_sig)))

#sets rowmeans to zero - heatmap = zscore above or below zero 
heat <- t(scale(t(vst)))

col=color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")),bias=0.75)(30)
#lower bias number gives more blues; higher bias gives more reds
#mess with bias to get white at 0

amo <- rbind(heat[grep('k141_1145589_26', rownames(heat)),], heat[grep('k141_1145589_29', rownames(heat)),], heat[grep('k141_1145589_28', rownames(heat)),])

rownames(amo)<- c('amoA', 'amoC', 'amoB') 



amo<-amo[,c(11,12,14,15,16,8,9,1,2,3,4,5,6,7,10,13)]


#reorder 
pheatmap(amo,color=col,cluster_cols=F,clustering_distance_rows="correlation")

#heatmap of sig genes
vst <- vst(dds_es_p2, blind=FALSE)

vst <- assay(vst)

vst <- as.data.frame(vst)

#pull out list of only significant genes

#vst_sig <- vst[rownames(vst) %in% significant_gene_names,]
#heat <- t(scale(t(vst_sig)))

#sets rowmeans to zero - heatmap = zscore above or below zero 
heat <- t(scale(t(vst)))

col=color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")),bias=0.75)(30)
#lower bias number gives more blues; higher bias gives more reds
#mess with bias to get white at 0

pheatmap(head(heat),color=col,cluster_cols=F,clustering_distance_rows="correlation")

#get only significantly different genes
#for es HvN
df1 <-data.frame(res_es_HvN_p2)
table(df1$padj < 0.1)

df1 %>% filter(padj < 0.1) -> sig_es_HvN

length(rownames(sig_es_HvN))
#20 :)

heat_sig_es_thaum_HvN <- heat[rownames(heat) %in% rownames(sig_es_HvN),]

length(rownames(heat_sig_es_thaum_HvN))
#20 :)


#reorder columns based on oxygen category (see Carly's script)
#norm, hyp, an
heat_sig_es_thaum_HvN_f<-heat_sig_es_thaum_HvN[,c(11,12,14,15,16,8,9,1,2,3,4,5,6,7,10,13)]

#shorten gene names
#rownames(heat_sig_hs_HvN_EEE_gn_final) <- strtrim(rownames(heat_sig_hs_HvN_EEE_gn_final), 75)

#make heatmap
pheatmap(heat_sig_es_thaum_HvN_f, color=col,cluster_cols=F,clustering_distance_rows="correlation")

#save csv to easily remake heatmap later
write.csv(heat_sig_es_thaum_HvN_f, 'es_thaum_heat_HvN.csv')

#AvN
df1 <-data.frame(res_es_AvN_p2)
table(df1$padj < 0.1)

df1 %>% filter(padj < 0.1) -> sig_es_AvN

length(rownames(sig_es_AvN))
#3 :)

heat_sig_es_thaum_AvN <- heat[rownames(heat) %in% rownames(sig_es_AvN),]

length(rownames(heat_sig_es_thaum_AvN))
#3 :)


#reorder columns based on oxygen category (see Carly's script)
#norm, hyp, an
heat_sig_es_thaum_AvN_f<-heat_sig_es_thaum_AvN[,c(11,12,14,15,16,8,9,1,2,3,4,5,6,7,10,13)]



#make heatmap
pheatmap(heat_sig_es_thaum_AvN_f, color=col,cluster_cols=F,clustering_distance_rows="correlation")

#save csv to easily remake heatmap later
write.csv(heat_sig_es_thaum_AvN_f, 'es_thaum_heat_AvN.csv')

######################
#make report with annotations
######################

#which comparisons is this running?
#how to specify with multiple factors 
#matches differnces between hypoxia and normoxia

des2Report <- HTMLReport(shortName = 'ES_prok_arch_RNAseq_analysis_with_DESeq2_long_names_HvN', title = 'ES prok arch Differential expression HvN',reportDirectory = "./reports")
publish(dds_es_p2,des2Report, pvalueCutoff=0.1, factor = colData(dds_es_p2)$Oxy_cat, reportDir="./reports")
finish(des2Report)

#try with specific contrasts 
#anoxia vs normoxia
#has a built in limit of 1000
des2Report <- HTMLReport(shortName = 'ES_prok_arch_RNAseq_analysis_with_DESeq2_long_names_AvN', title = 'ES prok arch Differential expression AvN',reportDirectory = "./reports")
publish(dds_es_p2,des2Report, pvalueCutoff=0.1, n=3000, factor = colData(dds_es_p2)$Oxy_cat, contrast =c("Oxy_cat", "Normoxic", "Anoxic"), reportDir="./reports")
finish(des2Report)


#anoxia vs hypoxia
des2Report <- HTMLReport(shortName = 'ES_prok_arch_RNAseq_analysis_with_DESeq2_longer_names_AvH', title = 'ES prok arch filtered Differential expression AvH',reportDirectory = "./reports")
publish(dds_es_p2,des2Report, pvalueCutoff=0.1, n=3000, factor = colData(dds_es_p2)$Oxy_cat, contrast =c("Oxy_cat", "Anoxic", "Hypoxic"), reportDir="./reports")
finish(des2Report)

#get KOGMWU input files for Es, arch
#just contig names for KOG anaylsis


#KOG directional file
df1 <-data.frame(res_es_HvA_p2)
df2 <- df1 

table(df2$padj < 0.1)
#avh LFC >0 -> upregulated in anoxia
df2$lfc <- df2$log2FoldChange

#make csv with just gene names and lfc
kog <- data.frame(rownames(df2))
kog$lfc <-df2$lfc 
names(kog)[1] <- 'gene_id'

#write as csv
write.csv(kog, file='es_kog_arch_HvA_lfc_final.csv', row.names = FALSE)

#kog HvN
df1 <-data.frame(res_es_HvN_p2)
df2 <- df1 

table(df2$padj < 0.1)
#LFC >0 -> upregulated in hypoxia
df2$lfc <- df2$log2FoldChange

#make csv with just gene names and lfc
kog <- data.frame(rownames(df2))
kog$lfc <-df2$lfc 
names(kog)[1] <- 'gene_id'

#write as csv
write.csv(kog, file='es_kog_arch_HvN_lfc_final.csv', row.names = FALSE)

#AvN
df1 <-data.frame(res_es_AvN_p2)
df2 <- df1 

table(df2$padj < 0.1)
#LFC >0 -> upregulated in normoxia (downreg in anoxia)
df2$lfc <- df2$log2FoldChange*-1

#make csv with just gene names and lfc
kog <- data.frame(rownames(df2))
kog$lfc <-df2$lfc 
names(kog)[1] <- 'gene_id'

#write as csv
write.csv(kog, file='es_kog_arch_AvN_lfc_final.csv', row.names = FALSE)

#add contig names + trim out unannotated, non-thaum (for sure) genes from gene2kog file
es_a_kog2gene=read.table("es_bin3_nitroso_contig2kogClass_final.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)

names(es_a_kog2gene)[names(es_a_kog2gene) == "V1"] <- "gene_id"

gn=read.table("es_bin3_nitroso_contig2geneName.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)

#gene_id
#gene_name

names(gn)[names(gn) == "V1"] <- "gene_id"

names(gn)[names(gn) == "V2"] <- "gene_name"

joined <- left_join(es_a_kog2gene, gn, by = "gene_id")
joined$full <-paste(joined$gene_id,joined$gene_name)

#reorder so that combined column is first (also dropp extra columns not needed for count matrix)
#also got rid of names column
joined2 <-joined[, c(4, 2)]

#shorten beforehand
#joined2$full 
joined2$full <- gsub(" ", "_", joined2$full)
#replace all special characters
joined2$full <- gsub("[()]", "_", joined2$full)
joined2$full <- gsub("[:]", "_", joined2$full)
joined2$full <- gsub("[/]", "_", joined2$full)
joined2$full <- gsub("[.]", "_", joined2$full)
joined2$full <- gsub("[,]", "_", joined2$full)
joined2$full <- gsub("[[]", "_", joined2$full)
joined2$full <- gsub("[]]", "_", joined2$full)
#limit rownames to 50 characters
joined2$full <- strtrim(joined2$full, 50)

es_arch_ann_kog <- joined2


#KOGMWU 
es_AvN_arch.lfc <- read.csv('es_kog_arch_AvN_lfc_final.csv', header =TRUE)

es_AvN_arch.lth=kog.mwu(es_AvN_arch.lfc,es_arch_ann_kog) 
es_AvN_arch.lth 

es_HvN_arch.lfc <- read.csv('es_kog_arch_HvN_lfc_final.csv', header =TRUE)

es_HvN_arch.lth=kog.mwu(es_HvN_arch.lfc,es_arch_ann_kog) 
es_HvN_arch.lth 

es_AvH_arch.lfc <- read.csv('es_kog_arch_HvA_lfc_final.csv', header =TRUE)

es_AvH_arch.lth=kog.mwu(es_AvH_arch.lfc,es_arch_ann_kog) 
es_AvH_arch.lth 
#no significantly enriched KOGs


#summarize annotated genes 


#seems like some of the SDEGs are driven by the outliers in HvN.
#I'll try getting rid of them just in case
##############

#######################
#Es gammas
######################
#Es - gammas
counts_esProk=read.table("allcounts_es_bin2_gamma_only.txt",header=TRUE,row.names=1) #Reading in the table of counts per isogroup by sample - you must open R in the same folder or change the working directory

counts_es=counts_esProk[,order(names(counts_esProk))]
names(counts_es)

##conditions (metadata)
conditions=read.table("Conditions_es.txt", header = T)

conds_es <- conditions[ order(row.names(conditions)), ]

####load data into deseq table format
#includes experimental design' since it is from the quick start method
DESeq2table_es_prok <- DESeqDataSetFromMatrix(countData = counts_es,
                                              colData = conds_es,
                                              design= ~ Oxy_cat)

#make normoxic the 'reference' level 
DESeq2table_es_prok$Oxy_cat <- relevel(DESeq2table_es_prok$Oxy_cat, "Normoxic")

#DESeq2Table_no_out2 <- DESeq2table_es_prok[, !(colnames(DESeq2table_es_prok) %in% c('DC108', 'DC110'))]

#change input below to keep outliers
dds_es_p_g <- DESeq(DESeq2table_es_prok)

resultsNames(dds_es_p_g) # lists the coefficients
res_es_AvN_p_g <- results(dds_es_p_g, name="Oxy_cat_Anoxic_vs_Normoxic")
# or to shrink log fold changes association with condition:

#shrunken results beased on LFC
res_es_AvN_apeglm_p_g <- lfcShrink(dds_es_p_g, coef="Oxy_cat_Anoxic_vs_Normoxic", type="apeglm")

#significantly expressed genes based on full results
table(res_es_AvN_p_g$padj < 0.1)

#FALSE  TRUE 
# 2062    38
#38 differentially expressed genes between normoxic and anoxic samples :)
#no outs

table(res_es_AvN_apeglm_p_g$padj < 0.1)

##57 differentially expressed genes between normoxic and anoxic samples = same as above :)

#standard method - hypoxic vs normoxic
res_es_HvN_p_g <- results(dds_es_p_g, name="Oxy_cat_Hypoxic_vs_Normoxic")

table(res_es_HvN_p_g$padj < 0.1)
#FALSE  TRUE 
#8608     
#0 differentially expressed genes between hypoxic and normoxic


# or to shrink log fold changes association with condition:
res_es_HvN_apeglm_p_g <- lfcShrink(dds_es_p_g, coef="Oxy_cat_Hypoxic_vs_Normoxic", type="apeglm")

table(res_es_HvN_apeglm_p_g$padj < 0.1)
#same as above

length(res_es_AvN_p_g[,1])

#8609 genes in total for es gamma
#too many for a bac, will take out all of the non gamma unannotated genes in next step

#comparison (contrast), between hypoxic and anoxic samples
res_es_HvA_p_g <- results(dds_es_p_g, contrast = list('Oxy_cat_Hypoxic_vs_Normoxic', 'Oxy_cat_Anoxic_vs_Normoxic'))

table(res_es_HvA_p_g$padj < 0.1)
#FALSE  TRUE 
#   2232    34 


#shrunken results 
#apeglm does not have contrast ability
#need to try 'normal shrinkage, or ashr

res_es_HvA_ashr_p_g <- lfcShrink(dds_es_p_g, contrast = list("Oxy_cat_Anoxic_vs_Normoxic", "Oxy_cat_Hypoxic_vs_Normoxic"), type="ashr")

table(res_es_HvA_ashr_p_g$padj < 0.1)
#same number as above :)

#more about contrasts here: https://github.com/tavareshugo/tutorial_DESeq2_contrasts/blob/main/DESeq2_contrasts.md
#also includes instructions for combining levels using a matrix structure


#intercept does not work for normoxic vs. hypoxic... how do we get these data?!
#checked https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
#install.packages('pheatmap')

#prep for PCA/heatmaps
#transformation 
vsd <- vst(dds_es_p_g, blind=FALSE)

df <- as.data.frame(colData(dds_es_p_g)[,c("Oxy_cat")])

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

#now anoxic samples group well together
#pca plot 
plotPCA(vsd, intgroup=c('Oxy_cat'))

#added ntop, flag, but one anoxic (again DC24) sample always grouped with the normoxic ones.

#lots of separation from one anoxic sample
#some separation with anoxia and normoxia..
#maybe would be better to use continuous oxygen rather than categorical

#could be effectively the same if top left outlier is removed from the PCA plot 

#MA plots

plotMA(res_es_AvN_p_g)

plotMA(res_es_AvN_apeglm_p_g)

plotMA(res_es_HvN_apeglm_p_g)

plotMA(res_es_HvA_ashr_p_g)

#dispersions also look better without outliers

##After calling plotMA, one can use the function identify to interactively detect the row number of individual genes by clicking on the plot. One can then recover the gene identifiers by saving the resulting indices
#idx <- identify(res_es_AvN_apeglm$baseMean, res_es_AvN_apeglm$log2FoldChange)
##print gene names from points you clicked on 
#rownames(res_es)[idx]

####possible removing of outliers 
boxplot(log10(assays(dds_es_p_g)[["cooks"]]), range=0, las=2)

#anoxic samples both look like outliers, but biologically are not (same as for euks)

#plotting dispersion estimates
plotDispEsts(dds_es_p)


#########
#summary of all degs for prok using package ReportingTools
#########

#add annotations
#work around for getting longer gene names (not just isogroup)
gn=read.table("es_bin2_gamma_contig2geneName.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)

#gene_id
#gene_name

names(gn)[names(gn) == "V1"] <- "gene_id"

names(gn)[names(gn) == "V2"] <- "gene_name"

head(gn)

#add annotation to dds_es object 
#check if they match in order 
#all(rownames(dds_es_p) == gn$gene_id)
#false, they do not match (expected)

#have to match them
#gn2 <- gn[match(rownames(dds_es2), gn$gene_id),]
#doesn'#t# work because many genes are not annotated

#try with left join? 
#need to get a count table with a column name for the row data to count 
counts_esProk2 <- cbind(rownames(counts_es), counts_es)
rownames(counts_esProk2) <- NULL
colnames(counts_esProk2) <- c("gene_id", 'DC108', 'DC110', 'DC112', 'DC114', 'DC118', 'DC123', 'DC131', 'DC135', 'DC136', 'DC24', 'DC33', 'DC39', 'DC56', 'DC85', 'DC93', 'DC97')

joined <- left_join(counts_esProk2, gn, by = "gene_id")
joined$full <-paste(joined$gene_id,joined$gene_name)

#reorder so that combined column is first (also dropp extra columns not needed for count matrix)
#also got rid of names column
joined2 <-joined[, c(19, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17)]

#make first column back into rowname
joined3 <- joined2[,-1]
rownames(joined3) <- joined2[,1]

#repeat deseq anaylsis with annotated count table 
counts_es=joined3[,order(names(joined3))]
names(counts_es)

#replace all spaces in row names (gene names) with underscore
rownames(counts_es) <- gsub(" ", "_", rownames(counts_es))
#replace all special characters
rownames(counts_es) <- gsub("[()]", "_", rownames(counts_es))
rownames(counts_es) <- gsub("[:]", "_", rownames(counts_es))
rownames(counts_es) <- gsub("[/]", "_", rownames(counts_es))
rownames(counts_es) <- gsub("[.]", "_", rownames(counts_es))
rownames(counts_es) <- gsub("[,]", "_", rownames(counts_es))
rownames(counts_es) <- gsub("[[]", "_", rownames(counts_es))
rownames(counts_es) <- gsub("[]]", "_", rownames(counts_es))
#limit rownames to 50 characters
rownames(counts_es) <- strtrim(rownames(counts_es), 50)

#get rid of no match, unannotated genes
counts_es_f <- counts_es[!grepl("NA", rownames(counts_es)),]


#includes experimental design' since it is from the quick start method
DESeq2table_es_names_p <- DESeqDataSetFromMatrix(countData = counts_es_f,
                                                 colData = conds_es,
                                                 design= ~ Oxy_cat)


DESeq2table_es_names_p$Oxy_cat <- relevel(DESeq2table_es_names_p$Oxy_cat, "Normoxic")

dds_es_g2 <- DESeq(DESeq2table_es_names_p)

#with releveling, of course it doesn't fit
resultsNames(dds_es_g2)

#hypoxia vs. anoxia
res_es_HvA_p2 <- results(dds_es_g2, contrast = list('Oxy_cat_Anoxic_vs_Normoxic', 'Oxy_cat_Hypoxic_vs_Normoxic'))

table(res_es_HvA_p2$padj < 0.1)

#same as above once releveled with normoxia as baseline.
#FALSE  TRUE 
#2494

#res_es3_HvA_ashr <- lfcShrink(dds_es3, contrast = list("Oxy_cat_Anoxic_vs_Normoxic", "Oxy_cat_Hypoxic_vs_Normoxic"), type="ashr")

#table(res_es3_HvA_ashr$padj < 0.1)
#same as above :)

#Normoxia vs anoxia
res_es_AvN_p2 <- results(dds_es_g2, name="Oxy_cat_Anoxic_vs_Normoxic")

table(res_es_AvN_p2$padj < 0.1)
# 2494
#0  differentially expressed genes between anoxic and normoxic samples :)

#shrunken results beased on LFC
#res_es3_AvN_apeglm <- lfcShrink(dds_es3, coef="Oxy_cat_Anoxic_vs_Normoxic", type="apeglm")

#significantly expressed genes based on full results
#table(res_es3_AvN_apeglm$padj < 0.1)


#significantly expressed genes based on LFC results
#table(res_es3_AvN_apeglm$padj < 0.1)
##need to run the other multiple comparisons
#standard method - hypoxic vs normoxic
res_es_HvN_p2 <- results(dds_es_g2, name="Oxy_cat_Hypoxic_vs_Normoxic")

table(res_es_HvN_p2$padj < 0.1)
#FALSE  TRUE 
#818    32 
#32 differentially expressed genes between hypoxic and normoxic

# or to shrink log fold changes association with condition:
#res_es3_HvN_apeglm <- lfcShrink(dds_es3, coef="Oxy_cat_Hypoxic_vs_Normoxic", type="apeglm")
#table(res_es3_HvN_apeglm$padj < 0.1)

#heatmap of sig genes
vst <- varianceStabilizingTransformation(dds_es_g2, blind=FALSE)

vst <- assay(vst)

vst <- as.data.frame(vst)

#pull out list of only significant genes

#vst_sig <- vst[rownames(vst) %in% significant_gene_names,]
#heat <- t(scale(t(vst_sig)))

#sets rowmeans to zero - heatmap = zscore above or below zero 
heat <- t(scale(t(vst)))

col=color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")),bias=0.75)(30)
#lower bias number gives more blues; higher bias gives more reds
#mess with bias to get white at 0

pheatmap(head(heat),color=col,cluster_cols=F,clustering_distance_rows="correlation")

#get only significantly different genes
#for es HvN
df1 <-data.frame(res_es_HvN_p2)
table(df1$padj < 0.1)

df1 %>% filter(padj < 0.1) -> sig_es_HvN

length(rownames(sig_es_HvN))
#32 :)

heat_sig_es_gamma_HvN <- heat[rownames(heat) %in% rownames(sig_es_HvN),]

length(rownames(heat_sig_es_gamma_HvN))
#32 :)


#reorder columns based on oxygen category (see Carly's script)
#norm, hyp, an
heat_sig_es_gamma_HvN_f<-heat_sig_es_gamma_HvN[,c(11,12,14,15,16,8,9,1,2,3,4,5,6,7,10,13)]

#shorten gene names
#rownames(heat_sig_hs_HvN_EEE_gn_final) <- strtrim(rownames(heat_sig_hs_HvN_EEE_gn_final), 75)

#make heatmap
pheatmap(heat_sig_es_gamma_HvN_f, color=col,cluster_cols=F,clustering_distance_rows="correlation")

#save csv to easily remake heatmap later
write.csv(heat_sig_es_gamma_HvN_f, 'es_gamma_heat_HvN.csv')

#no sig difs in the other two comparisons 

######################
#make report with annotations
######################

#which comparisons is this running?
#how to specify with multiple factors 
#only shows 608 entries
#matches differnces between hypoxia and normoxia

des2Report <- HTMLReport(shortName = 'ES_prok_gamma_RNAseq_analysis_with_DESeq2_long_names_HvN_no_outs', title = 'ES prok gamma Differential expression HvN_no_outs',reportDirectory = "./reports")
publish(dds_es_g2,des2Report, pvalueCutoff=0.1, factor = colData(dds_es_g2)$Oxy_cat, reportDir="./reports")
finish(des2Report)

#try with specific contrasts 
#anoxia vs normoxia
#get visual for lfc signs by accounting for NAs as well
des2Report <- HTMLReport(shortName = 'ES_prok_gamma_RNAseq_analysis_with_DESeq2_long_names_AvN_no_outs', title = 'ES prok gamma Differential expression AvN_no_outs',reportDirectory = "./reports")
publish(dds_es_g2,des2Report, pvalueCutoff='NA', n=3000, factor = colData(dds_es_g2)$Oxy_cat, contrast =c("Oxy_cat", "Normoxic", "Anoxic"), reportDir="./reports")
finish(des2Report)


#anoxia vs hypoxia
#has a built in limit of 1000
des2Report <- HTMLReport(shortName = 'ES_prok_gamma_RNAseq_analysis_with_DESeq2_longer_names_AvH_no_outs', title = 'ES prok gamma filtered Differential expression AvH_no_outs',reportDirectory = "./reports")
publish(dds_es_g2,des2Report, pvalueCutoff=0.9999999999999, n=3000, factor = colData(dds_es_g2)$Oxy_cat, contrast =c("Oxy_cat", "Anoxic", "Hypoxic"), reportDir="./reports")
finish(des2Report)

#tomorrow: start by removing hypoxic outliers? -makes proks uniform
#also check annotations

####KOG prep ES gamma
#KOG directional file
df1 <-data.frame(res_es_HvA_p2)
df2 <- df1 

table(df2$padj < 0.1)
#avh LFC >0 -> upregulated in anoxia
df2$lfc <- df2$log2FoldChange

#make csv with just gene names and lfc
kog <- data.frame(rownames(df2))
kog$lfc <-df2$lfc 
names(kog)[1] <- 'gene_id'

#write as csv
write.csv(kog, file='es_kog_gamma_HvA_lfc_final.csv', row.names = FALSE)

#kog HvN
df1 <-data.frame(res_es_HvN_p2)
df2 <- df1 

table(df2$padj < 0.1)
#LFC >0 -> upregulated in hypoxia
df2$lfc <- df2$log2FoldChange

#make csv with just gene names and lfc
kog <- data.frame(rownames(df2))
kog$lfc <-df2$lfc 
names(kog)[1] <- 'gene_id'

#write as csv
write.csv(kog, file='es_kog_gamma_HvN_lfc_final.csv', row.names = FALSE)

#AvN
df1 <-data.frame(res_es_AvN_p2)
df2 <- df1 

table(df2$padj < 0.1)
#LFC >0 -> upregulated in normoxia (downreg in anoxia)
df2$lfc <- df2$log2FoldChange*-1

#make csv with just gene names and lfc
kog <- data.frame(rownames(df2))
kog$lfc <-df2$lfc 
names(kog)[1] <- 'gene_id'

#write as csv
write.csv(kog, file='es_kog_gamma_AvN_lfc_final.csv', row.names = FALSE) 


#KOG_MWU es gamma
es_g_kog2gene=read.table("es_bin2_gamma_contig2kogClass_final.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)

names(es_g_kog2gene)[names(es_g_kog2gene) == "V1"] <- "gene_id"

gn=read.table("es_bin2_gamma_contig2geneName.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)

#gene_id
#gene_name

names(gn)[names(gn) == "V1"] <- "gene_id"

names(gn)[names(gn) == "V2"] <- "gene_name"

joined <- left_join(es_g_kog2gene, gn, by = "gene_id")
joined$full <-paste(joined$gene_id,joined$gene_name)

#reorder so that combined column is first (also dropp extra columns not needed for count matrix)
#also got rid of names column
joined2 <-joined[, c(4, 2)]

#shorten beforehand
#joined2$full 
joined2$full <- gsub(" ", "_", joined2$full)
#replace all special characters
joined2$full <- gsub("[()]", "_", joined2$full)
joined2$full <- gsub("[:]", "_", joined2$full)
joined2$full <- gsub("[/]", "_", joined2$full)
joined2$full <- gsub("[.]", "_", joined2$full)
joined2$full <- gsub("[,]", "_", joined2$full)
joined2$full <- gsub("[[]", "_", joined2$full)
joined2$full <- gsub("[]]", "_", joined2$full)
#limit rownames to 50 characters
joined2$full <- strtrim(joined2$full, 50)

es_gamma_ann_kog <- joined2


#KOGMWU 
es_AvN_gamma.lfc <- read.csv('es_kog_gamma_AvN_lfc_final.csv', header =TRUE)

es_AvN_gamma.lth=kog.mwu(es_AvN_gamma.lfc,es_gamma_ann_kog) 
es_AvN_gamma.lth 

es_HvN_gamma.lfc <- read.csv('es_kog_gamma_HvN_lfc_final.csv', header =TRUE)

es_HvN_gamma.lth=kog.mwu(es_HvN_gamma.lfc,es_gamma_ann_kog) 
es_HvN_gamma.lth 

es_AvH_gamma.lfc <- read.csv('es_kog_gamma_HvA_lfc_final.csv', header =TRUE)

es_AvH_gamma.lth=kog.mwu(es_AvH_gamma.lfc,es_gamma_ann_kog) 
es_AvH_gamma.lth 


###################
###hs archaea
###################
#setwd("C:/Users/strehlow/OneDrive - Syddansk Universitet/SDU/Transcriptomics/DGE_lough_hyne")

countsHost=read.table("allcounts_hs_prok_bin1_nitroso_only.txt",header=TRUE,row.names=1) #Reading in the table of counts per isogroup by sample - you must open R in the same folder or change the working directory

counts=countsHost[,order(names(countsHost))]
names(counts)

##conditions (metadata)
conditions=read.table("Conditions_hs.txt", header = T)

conds <- conditions[ order(row.names(conditions)), ]

####load data into deseq table format
#includes experimental design' since it is from the quick start method
DESeq2Table <- DESeqDataSetFromMatrix(countData = countsHost,
                                      colData = conds,
                                      design= ~ Oxy_cat)

#make normoxic the 'reference' level 
DESeq2Table$Oxy_cat <- relevel(DESeq2Table$Oxy_cat, "Normoxic")

#runs differential expression anaylsis quickly
dds <- DESeq(DESeq2Table)

resultsNames(dds) # lists the coefficients
res_AvN <- results(dds, name="Oxy_cat_Anoxic_vs_Normoxic")
# or to shrink log fold changes association with condition:

#shrunken results beased on LFC
res_AvN_apeglm <- lfcShrink(dds, coef="Oxy_cat_Anoxic_vs_Normoxic", type="apeglm")

#significantly expressed genes based on full results
table(res_AvN$padj < 0.1)

#FALSE  TRUE 
#1478   0 
#0 differentially expressed genes between normoxic and anoxic samples :)

#significantly expressed genes based on LFC results
table(res_AvN_apeglm$padj < 0.1)

#standard method - hypoxic vs normoxic
res_HvN <- results(dds, name="Oxy_cat_Hypoxic_vs_Normoxic")

table(res_HvN$padj < 0.1)
#FALSE  TRUE 
#  399    17
#17 differentially expressed genes between hypoxic and normoxic

# or to shrink log fold changes association with condition:
res_HvN_apeglm <- lfcShrink(dds, coef="Oxy_cat_Hypoxic_vs_Normoxic", type="apeglm")

table(res_HvN_apeglm$padj < 0.1)
#same as above

length(res_HvN[,1])

#1481 genes in total for HS

#comparison (contrast), between hypoxic and anoxic samples
res_HvA <- results(dds, contrast = list('Oxy_cat_Hypoxic_vs_Normoxic', 'Oxy_cat_Anoxic_vs_Normoxic'))

table(res_HvA$padj < 0.1)
#FALSE  TRUE 
#  445 28

#28 genes differentially expressed between hypoxic and anoxic

#shrunken results 
#apeglm does not have contrast ability
#need to try 'normal shrinkage, or ashr

res_HvA_ashr <- lfcShrink(dds, contrast = list("Oxy_cat_Anoxic_vs_Normoxic", "Oxy_cat_Hypoxic_vs_Normoxic"), type="ashr")

table(res_HvA_ashr$padj < 0.1)

#same number as above :)

#more about contrasts here: https://github.com/tavareshugo/tutorial_DESeq2_contrasts/blob/main/DESeq2_contrasts.md
#also includes instructions for combining levels using a matrix structure


#intercept does not work for normoxic vs. hypoxic... how do we get these data?!
#checked https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
#install.packages('pheatmap')

#prep for PCA/heatmaps
#transformation 
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

df <- as.data.frame(colData(dds)[,c("Oxy_cat")])

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
plotPCA(vsd, intgroup=c('Oxy_cat'))

#good separation of the anoxic sample from the hypoxic and normoxic ones
#unfortunately only one normoxic sample :(
#however, given that the sponges were behaving normally under hypoxia (pumping), hypoxic could be pretty close to normal
#so hypoxic dot on the left could be an outlier...

#could be effectively the same if top left outlier is removed from the PCA plot 

#MA plots

plotMA(res_AvN)

plotMA(res_AvN_apeglm)

plotMA(res_HvN_apeglm)

plotMA(res_HvA_ashr)

##After calling plotMA, one can use the function identify to interactively detect the row number of individual genes by clicking on the plot. One can then recover the gene identifiers by saving the resulting indices
#idx <- identify(res_AvN_apeglm$baseMean, res_AvN_apeglm$log2FoldChange)
##print gene names from points you clicked on 
#rownames(res)[idx]

####possible removing of outliers 
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)

#plotting dispersion estimates
plotDispEsts(dds)


#PCA to remove outlier
PCA <- plotPCA(vsd, intgroup=c('Oxy_cat'), returnData=TRUE)
PCA

#get a report for direction of log change reference
des2Report <- HTMLReport(shortName = 'HS_prok_arch_RNAseq_analysis_with_DESeq2_long_names_HvN_alldata', title = 'HS arch Differential expression HvN',reportDirectory = "./reports")
publish(dds,des2Report, pvalueCutoff=0.1, factor = colData(dds)$Oxy_cat, contrast =c("Oxy_cat", "Normoxic", "Hypoxic"), reportDir="./reports_scratch")
finish(des2Report)

#try with specific contrasts 
#anoxia vs normoxia
#zero degs! 
des2Report <- HTMLReport(shortName = 'HS_prok_arch_RNAseq_analysis_with_DESeq2_long_names_AvN', title = 'HS prok arch Differential expression AvN',reportDirectory = "./reports")
publish(dds,des2Report, pvalueCutoff=0.9999999999, n=100, factor = colData(dds)$Oxy_cat, reportDir="./reports_scratch")
finish(des2Report)
#made report just to see the direction of change


#anoxia vs hypoxia
#has a built in limit of 1000
des2Report <- HTMLReport(shortName = 'HS_prok_arch_RNAseq_analysis_with_DESeq2_longer_names_AvH', title = 'HS prok arch filtered Differential expression AvH',reportDirectory = "./reports")
publish(dds,des2Report, pvalueCutoff=0.1, n=3000, factor = colData(dds)$Oxy_cat, contrast =c("Oxy_cat", "Anoxic", "Hypoxic"), reportDir="./reports_scratch")
finish(des2Report)




#########
#summary of all degs using package ReportingTools
#########

#add annotations
#work around for getting longer gene names (not just isogroup)
gn=read.table("hs_bin1_nitroso_contig2geneName.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)

#gene_id
#gene_name

names(gn)[names(gn) == "V1"] <- "gene_id"

names(gn)[names(gn) == "V2"] <- "gene_name"


#try with left join? 
#need to get a count table with a column name for the row data to count 
countsHost2 <- cbind(rownames(countsHost), countsHost)
rownames(countsHost2) <- NULL
colnames(countsHost2) <- c("gene_id", "DC105", "DC107", "DC116", "DC121", "DC129", "DC133", "DC54",  "DC55")

joined <- left_join(countsHost2, gn, by = "gene_id")
joined$full <-paste(joined$gene_id,joined$gene_name)

#reorder so that combined column is first (also dropp extra columns not needed for count matrix)
joined2 <-joined[, c(11, 2, 3, 4, 5, 6, 7, 8, 9)]

#make first column back into rowname
joined3 <- joined2[,-1]
rownames(joined3) <- joined2[,1]

#repeat deseq anaylsis with annotated count table 
counts=joined3[,order(names(joined3))]
names(counts)

#replace all spaces in row names (gene names) with underscore
rownames(counts) <- gsub(" ", "_", rownames(counts))
#replace all special characters
rownames(counts) <- gsub("[()]", "_", rownames(counts))
rownames(counts) <- gsub("[:]", "_", rownames(counts))
rownames(counts) <- gsub("[/]", "_", rownames(counts))
rownames(counts) <- gsub("[.]", "_", rownames(counts))
rownames(counts) <- gsub("[,]", "_", rownames(counts))

#limit rownames to 50 characters
rownames(counts) <- strtrim(rownames(counts), 50)


DESeq2Table_names <- DESeqDataSetFromMatrix(countData = counts,
                                            colData = conds,
                                            design= ~ Oxy_cat)
#remove outliers again
#drop normoxic as level 
#DESeq2Table_names$Oxy_cat <- droplevels.factor(DESeq2Table_names$Oxy_cat, "Normoxic")

#relevel with hypoxia as reference 
DESeq2Table_names$Oxy_cat <- relevel(DESeq2Table_names$Oxy_cat, "Normoxic")

#DESeq2Table_no_out2 <- DESeq2Table_names[, !(colnames(DESeq2Table_names) %in% c('DC129'))]
dds3 <- DESeq(DESeq2Table_names)

resultsNames(dds3) # lists the coefficients
res_AvN <- results(dds3, name="Oxy_cat_Anoxic_vs_Normoxic")
# or to shrink log fold changes association with condition:

#significantly expressed genes based on full results
table(res_AvN$padj < 0.1)

#FALSE  TRUE 
#1478   0 
#0 differentially expressed genes between normoxic and anoxic samples :)

#standard method - hypoxic vs normoxic
res_HvN <- results(dds3, name = 'Oxy_cat_Hypoxic_vs_Normoxic')

table(res_HvN$padj < 0.1)
#FALSE  TRUE 
#  399    17
#17 differentially expressed genes between hypoxic and normoxic

length(res_HvN[,1])

#comparison (contrast), between hypoxic and anoxic samples
res_HvA <- results(dds3, contrast = list("Oxy_cat_Anoxic_vs_Normoxic", "Oxy_cat_Hypoxic_vs_Normoxic"))

table(res_HvA$padj < 0.1)
#FALSE  TRUE 
#  445  28



#make report with annotations

des2Report <- HTMLReport(shortName = 'HS_prok_arch_RNAseq_analysis_with_DESeq2_long_names_HvN_alldata', title = 'HS arch Differential expression HvN',reportDirectory = "./reports")
publish(dds3,des2Report, pvalueCutoff=0.1, factor = colData(dds3)$Oxy_cat, contrast =c("Oxy_cat", "Normoxic", "Hypoxic"), reportDir="./reports")
finish(des2Report)

#try with specific contrasts 
#anoxia vs normoxia
#zero degs! 
des2Report <- HTMLReport(shortName = 'HS_prok_arch_RNAseq_analysis_with_DESeq2_long_names_AvN', title = 'HS prok arch Differential expression AvN',reportDirectory = "./reports")
publish(dds3,des2Report, pvalueCutoff=0.1, factor = colData(dds3)$Oxy_cat, reportDir="./reports")
finish(des2Report)
#made report just to see the direction of change


#anoxia vs hypoxia
#has a built in limit of 1000
des2Report <- HTMLReport(shortName = 'HS_prok_arch_RNAseq_analysis_with_DESeq2_longer_names_AvH', title = 'HS prok arch filtered Differential expression AvH',reportDirectory = "./reports")
publish(dds3,des2Report, pvalueCutoff=0.1, n=3000, factor = colData(dds3)$Oxy_cat, contrast =c("Oxy_cat", "Anoxic", "Hypoxic"), reportDir="./reports")
finish(des2Report)

###amo a,b and c genes?
#k141_490446_10	Archaeal ammonia monooxygenase subunit A (AmoA)
#k141_490446_12	PFAM Ammonia monooxygenase methane monooxygenase, subunit C
#k141_490446_13	Monooxygenase subunit B protein

#heatmap?
vst <- varianceStabilizingTransformation(dds3, blind=FALSE)

vst <- assay(vst)

vst <- as.data.frame(vst)

#pull out list of only significant genes

#vst_sig <- vst[rownames(vst) %in% significant_gene_names,]
#heat <- t(scale(t(vst_sig)))

#sets rowmeans to zero - heatmap = zscore above or below zero 
heat <- t(scale(t(vst)))

col=color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")),bias=0.75)(30)
#lower bias number gives more blues; higher bias gives more reds
#mess with bias to get white at 0

amo <- rbind(heat[grep('k141_490446_10', rownames(heat)),], heat[grep('k141_490446_12', rownames(heat)),], heat[grep('k141_490446_13', rownames(heat)),])

rownames(amo)<- c('amoA', 'amoC', 'amoB') 

amo<-amo[,c(5,1,2,3,4,6,7,8)]

#reorder 
pheatmap(amo,color=col,cluster_cols=F,clustering_distance_rows="correlation")


#heatmap of sig genes
vst <- varianceStabilizingTransformation(dds3, blind=FALSE)

vst <- assay(vst)

vst <- as.data.frame(vst)

#pull out list of only significant genes

#vst_sig <- vst[rownames(vst) %in% significant_gene_names,]
#heat <- t(scale(t(vst_sig)))

#sets rowmeans to zero - heatmap = zscore above or below zero 
heat <- t(scale(t(vst)))

col=color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")),bias=0.75)(30)
#lower bias number gives more blues; higher bias gives more reds
#mess with bias to get white at 0

pheatmap(head(heat),color=col,cluster_cols=F,clustering_distance_rows="correlation")

#get only significantly different genes
#for es HvN
df1 <-data.frame(res_HvN)
table(df1$padj < 0.1)

df1 %>% filter(padj < 0.1) -> sig_hs_HvN

length(rownames(sig_hs_HvN))
#17 :)

heat_sig_hs_thaum_HvN <- heat[rownames(heat) %in% rownames(sig_hs_HvN),]

length(rownames(heat_sig_hs_thaum_HvN))
#17 :)


#reorder columns based on oxygen category (see Carly's script)
#norm, hyp, an
heat_sig_hs_thaum_HvN_f<-heat_sig_hs_thaum_HvN[,c(5,1,2,3,4,6,7,8)]

#shorten gene names
#rownames(heat_sig_hs_HvN_EEE_gn_final) <- strtrim(rownames(heat_sig_hs_HvN_EEE_gn_final), 75)

#make heatmap
pheatmap(heat_sig_hs_thaum_HvN_f, color=col,cluster_cols=F,clustering_distance_rows="correlation")

#save csv to easily remake heatmap later
write.csv(heat_sig_hs_thaum_HvN_f, 'hs_thaum_heat_HvN.csv')

#AvN
#no sig differences 
df1 <-data.frame(res_HvA)
table(df1$padj < 0.1)

df1 %>% filter(padj < 0.1) -> sig_hs_AvH

length(rownames(sig_hs_AvH))
#28 :)

heat_sig_hs_thaum_AvH <- heat[rownames(heat) %in% rownames(sig_hs_AvH),]

length(rownames(heat_sig_hs_thaum_AvH))
#28 :)


#reorder columns based on oxygen category (see Carly's script)
#norm, hyp, an
heat_sig_hs_thaum_AvH_f<-heat_sig_hs_thaum_AvH[,c(5,1,2,3,4,6,7,8)]


#make heatmap
pheatmap(heat_sig_hs_thaum_AvH_f, color=col,cluster_cols=F,clustering_distance_rows="correlation")

#save csv to easily remake heatmap later
write.csv(heat_sig_hs_thaum_AvH_f, 'hs_thaum_heat_AvH.csv')

##kogs
#kog_mwu with all HS arch data based on log fold change, without no match genes 
#KOG directional file
df1 <-data.frame(res_HvA)
df2 <- df1 

table(df2$padj < 0.1)
#avh LFC >0 -> upregulated in anoxia
df2$lfc <- df2$log2FoldChange

#make csv with just gene names and lfc
kog <- data.frame(rownames(df2))
kog$lfc <-df2$lfc 
names(kog)[1] <- 'gene_id'

#write as csv
write.csv(kog, file='hs_kog_arch_HvA_lfc_final.csv', row.names = FALSE)

#kog HvN
df1 <-data.frame(res_HvN)
df2 <- df1 

table(df2$padj < 0.1)
#LFC <0 -> upregulated in hypoxia
df2$lfc <- df2$log2FoldChange*(-1)

#make csv with just gene names and lfc
kog <- data.frame(rownames(df2))
kog$lfc <-df2$lfc 
names(kog)[1] <- 'gene_id'

#write as csv
write.csv(kog, file='hs_kog_arch_HvN_lfc_final.csv', row.names = FALSE)

#AvN
df1 <-data.frame(res_AvN)
df2 <- df1 

table(df2$padj < 0.1)
#LFC >0 -> upregulated in normoxia (downreg in anoxia)
df2$lfc <- df2$log2FoldChange*(-1)

#make csv with just gene names and lfc
kog <- data.frame(rownames(df2))
kog$lfc <-df2$lfc 
names(kog)[1] <- 'gene_id'

#write as csv
write.csv(kog, file='hs_kog_arch_AvN_lfc_final.csv', row.names = FALSE)

# Import KOG table Hymeraphia arch

hs_arch_kog2gene=read.table("hs_bin1_nitroso_contig2kogClass_final.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)

###add gene names to kog class file to match

names(hs_arch_kog2gene)[names(hs_arch_kog2gene) == "V1"] <- "gene_id"

gn=read.table("hs_bin1_nitroso_contig2geneName.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)

#gene_id
#gene_name

names(gn)[names(gn) == "V1"] <- "gene_id"

names(gn)[names(gn) == "V2"] <- "gene_name"

joined <- left_join(hs_arch_kog2gene, gn, by = "gene_id")
joined$full <-paste(joined$gene_id,joined$gene_name)

#reorder so that combined column is first (also dropp extra columns not needed for count matrix)
#also got rid of names column
joined2 <-joined[, c(4, 2)]

#shorten beforehand
#joined2$full 
joined2$full <- gsub(" ", "_", joined2$full)
#replace all special characters
joined2$full <- gsub("[()]", "_", joined2$full)
joined2$full <- gsub("[:]", "_", joined2$full)
joined2$full <- gsub("[/]", "_", joined2$full)
joined2$full <- gsub("[.]", "_", joined2$full)
joined2$full <- gsub("[,]", "_", joined2$full)
joined2$full <- gsub("[[]", "_", joined2$full)
joined2$full <- gsub("[]]", "_", joined2$full)
#limit rownames to 50 characters
joined2$full <- strtrim(joined2$full, 50)

hs_arch_ann_kog <- joined2


#lfc for comparison
hs_AvN.lfc <- read.csv('hs_kog_arch_AvN_lfc_final.csv', header =TRUE)

hs_AvN.lth=kog.mwu(hs_AvN.lfc,hs_arch_ann_kog) 
hs_AvN.lth 

#lfc (can switch them out to see the diffrences)
hs_HvN.lfc <- read.csv('hs_kog_arch_HvN_lfc_final.csv', header =TRUE)
hs_HvN.lth=kog.mwu(hs_HvN.lfc,hs_arch_ann_kog)
hs_HvN.lth

# hs HvA
hs_HvA.lfc <- read.csv('hs_kog_arch_HvA_lfc_final.csv', header =TRUE)
hs_HvA.lth=kog.mwu(hs_HvA.lfc,hs_arch_ann_kog)
hs_HvA.lth



# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("Es_Anox"=es_AvN_arch.lth, "Es_Hypox"=es_HvN_arch.lth, 'Es_Hypox_vs_Anox'=es_AvH_arch.lth,   "Hs_Anox"=hs_AvN.lth,"Hs_Hypox"=hs_HvN.lth,"Hs_Hypox_vs_Anox"=hs_HvA.lth))

#plot KOG differences
pheatmap(as.matrix(ktable), clustering_distance_cols="correlation") 
pheatmap(as.matrix(ktable), cluster_cols=F,clustering_distance_rows="correlation", clustering_distance_cols="correlation")


pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor)
# p-values of these correlations in the upper panel:
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor.pval)

#only 18 KOGs (COGs) b/c proks don't have things like chromatin, etc

#save KOG csv files
write.csv(hs_HvA.lth, 'hs_arch_KOG_HvA.csv')
write.csv(hs_HvN.lth, 'hs_arch_KOG_HvN.csv')
write.csv(hs_AvN.lth, 'hs_arch_KOG_AvN.csv')

#es (move up)
write.csv(es_AvH_arch.lth, 'es_arch_KOG_HvA.csv')
write.csv(es_HvN_arch.lth, 'es_arch_KOG_HvN.csv')
write.csv(es_AvN_arch.lth, 'es_arch_KOG_AvN.csv')



###########
#hymeraphia gammaproteobacteria
###########

countsHost=read.table("allcounts_hs_prok_bin2_gamma_only.txt",header=TRUE,row.names=1) #Reading in the table of counts per isogroup by sample - you must open R in the same folder or change the working directory

counts=countsHost[,order(names(countsHost))]
names(counts)

##conditions (metadata)
conditions=read.table("Conditions_hs.txt", header = T)

conds <- conditions[ order(row.names(conditions)), ]

####load data into deseq table format
#includes experimental design' since it is from the quick start method
DESeq2Table <- DESeqDataSetFromMatrix(countData = counts,
                                      colData = conds,
                                      design= ~ Oxy_cat)

#make normoxic the 'reference' level 
DESeq2Table$Oxy_cat <- relevel(DESeq2Table$Oxy_cat, "Normoxic")

#runs differential expression anaylsis quickly
dds <- DESeq(DESeq2Table)

resultsNames(dds) # lists the coefficients
res_AvN <- results(dds, name="Oxy_cat_Anoxic_vs_Normoxic")
# or to shrink log fold changes association with condition:

#shrunken results beased on LFC
res_AvN_apeglm <- lfcShrink(dds, coef="Oxy_cat_Anoxic_vs_Normoxic", type="apeglm")

#significantly expressed genes based on full results
table(res_AvN$padj < 0.1)

#FALSE  TRUE 
# 8409
#0 differentially expressed genes between normoxic and anoxic samples :)

#significantly expressed genes based on LFC results
table(res_AvN_apeglm$padj < 0.1)

#standard method - hypoxic vs normoxic
res_HvN <- results(dds, name="Oxy_cat_Hypoxic_vs_Normoxic")

table(res_HvN$padj < 0.1)
#FALSE  TRUE 
#4023 304

# or to shrink log fold changes association with condition:
res_HvN_apeglm <- lfcShrink(dds, coef="Oxy_cat_Hypoxic_vs_Normoxic", type="apeglm")

table(res_HvN_apeglm$padj < 0.1)
#same as above

length(res_HvN[,1])

#1828 genes in total for HS

#comparison (contrast), between hypoxic and anoxic samples
res_HvA <- results(dds, contrast = list('Oxy_cat_Hypoxic_vs_Normoxic', 'Oxy_cat_Anoxic_vs_Normoxic'))

table(res_HvA$padj < 0.1)
#FALSE  TRUE 
# 4694  1429 

res_HvA_ashr <- lfcShrink(dds, contrast = list("Oxy_cat_Anoxic_vs_Normoxic", "Oxy_cat_Hypoxic_vs_Normoxic"), type="ashr")

table(res_HvA_ashr$padj < 0.1)

#same number as above :)

#more about contrasts here: https://github.com/tavareshugo/tutorial_DESeq2_contrasts/blob/main/DESeq2_contrasts.md
#also includes instructions for combining levels using a matrix structure


#intercept does not work for normoxic vs. hypoxic... how do we get these data?!
#checked https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
#install.packages('pheatmap')

#prep for PCA/heatmaps
#transformation 
vsd <- vst(dds, blind=FALSE)

df <- as.data.frame(colData(dds)[,c("Oxy_cat")])

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
plotPCA(vsd, intgroup=c('Oxy_cat'))


#could be effectively the same if top left outlier is removed from the PCA plot 

#MA plots

plotMA(res_AvN)

plotMA(res_AvN_apeglm)

plotMA(res_HvN_apeglm)

plotMA(res_HvA_ashr)

##After calling plotMA, one can use the function identify to interactively detect the row number of individual genes by clicking on the plot. One can then recover the gene identifiers by saving the resulting indices
#idx <- identify(res_AvN_apeglm$baseMean, res_AvN_apeglm$log2FoldChange)
##print gene names from points you clicked on 
#rownames(res)[idx]

####possible removing of outliers 
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)

#sample DC129 has consistantly higher outliers 
#this is the normoxic sample though... so it is probably best to include it.

#plotting dispersion estimates
plotDispEsts(dds)


#PCA to remove outlier
PCA <- plotPCA(vsd, intgroup=c('Oxy_cat'), returnData=TRUE)
PCA

#add annotations
#work around for getting longer gene names (not just isogroup)
gn=read.table("hs_bin2_gamma_final_contig2geneName.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)

#gene_id
#gene_name

names(gn)[names(gn) == "V1"] <- "gene_id"

names(gn)[names(gn) == "V2"] <- "gene_name"

head(gn)


#try with left join? 
#need to get a count table with a column name for the row data to count 
countsHost2 <- cbind(rownames(countsHost), countsHost)
rownames(countsHost2) <- NULL
colnames(countsHost2) <- c("gene_id", "DC105", "DC107", "DC116", "DC121", "DC129", "DC133", "DC54",  "DC55")

joined <- left_join(countsHost2, gn, by = "gene_id")
joined$full <-paste(joined$gene_id,joined$gene_name)

#reorder so that combined column is first (also dropp extra columns not needed for count matrix)
joined2 <-joined[, c(11, 2, 3, 4, 5, 6, 7, 8, 9)]

#make first column back into rowname
joined3 <- joined2[,-1]
rownames(joined3) <- joined2[,1]

#repeat deseq anaylsis with annotated count table 
counts=joined3[,order(names(joined3))]
names(counts)

#replace all spaces in row names (gene names) with underscore
rownames(counts) <- gsub(" ", "_", rownames(counts))
#replace all special characters
rownames(counts) <- gsub("[()]", "_", rownames(counts))
rownames(counts) <- gsub("[:]", "_", rownames(counts))
rownames(counts) <- gsub("[/]", "_", rownames(counts))
rownames(counts) <- gsub("[.]", "_", rownames(counts))
rownames(counts) <- gsub("[,]", "_", rownames(counts))

#limit rownames to 50 characters
rownames(counts) <- strtrim(rownames(counts), 50)

###remove no matches to gamma with no annotation either
counts_hs_f <- counts[!grepl("NA", rownames(counts)),]

#adding longer names
####load data into deseq table format
#includes experimental design' since it is from the quick start method
DESeq2Table_names <- DESeqDataSetFromMatrix(countData = counts_hs_f,
                                            colData = conds,
                                            design= ~ Oxy_cat)
#remove outliers again
#drop normoxic as level 
#DESeq2Table_names$Oxy_cat <- droplevels.factor(DESeq2Table_names$Oxy_cat, "Normoxic")

#relevel with hypoxia as reference 
DESeq2Table_names$Oxy_cat <- relevel(DESeq2Table_names$Oxy_cat, "Normoxic")

#DESeq2Table_no_out2 <- DESeq2Table_names[, !(colnames(DESeq2Table_names) %in% c('DC129'))]
dds3 <- DESeq(DESeq2Table_names)

boxplot(log10(assays(dds3)[["cooks"]]), range=0, las=2)
#only 129 (normoxia) is potential outlier, as expected

resultsNames(dds3) # lists the coefficients
res3_AvH <- results(dds3, contrast = list('Oxy_cat_Hypoxic_vs_Normoxic', 'Oxy_cat_Anoxic_vs_Normoxic'))

#significantly expressed genes based on full results
table(res3_AvH$padj < 0.1)

#FALSE  TRUE 
# 974   436
#436 differentially expressed genes between hypoxic and anoxic samples :)

res3_HvN <- results(dds3, name="Oxy_cat_Hypoxic_vs_Normoxic")
table(res3_HvN$padj < 0.1)
#0

res3_AvN <- results(dds3, name="Oxy_cat_Anoxic_vs_Normoxic")
table(res3_AvN$padj < 0.1)
#0 

#heatmap of sig genes
vst <- varianceStabilizingTransformation(dds3, blind=FALSE)

vst <- assay(vst)

vst <- as.data.frame(vst)

#pull out list of only significant genes

#vst_sig <- vst[rownames(vst) %in% significant_gene_names,]
#heat <- t(scale(t(vst_sig)))

#sets rowmeans to zero - heatmap = zscore above or below zero 
#changed sign because reports and heatmap did not aggree in their direction of expression
heat <- -t(scale(t(vst)))

col=color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")),bias=0.75)(30)
#lower bias number gives more blues; higher bias gives more reds
#mess with bias to get white at 0

pheatmap(head(heat),color=col,cluster_cols=F,clustering_distance_rows="correlation")

#only sig diffs in AvH 

#AvN
#no sig differences 
df1 <-data.frame(res3_AvH)
table(df1$padj < 0.1)

df1 %>% filter(padj < 0.1) -> sig_hs_AvH

length(rownames(sig_hs_AvH))
#436 :)

heat_sig_hs_gamma_AvH <- heat[rownames(heat) %in% rownames(sig_hs_AvH),]

length(rownames(heat_sig_hs_gamma_AvH))



#reorder columns based on oxygen category (see Carly's script)
#norm, hyp, an
heat_sig_hs_gamma_AvH_f<-heat_sig_hs_gamma_AvH[,c(5,1,2,3,4,6,7,8)]


#make heatmap
pheatmap(heat_sig_hs_gamma_AvH_f, color=col, cluster_cols=F,clustering_distance_rows="correlation")

#just first 20 (upregulated)?
pheatmap(head(heat_sig_hs_gamma_AvH_f), color=col, cluster_cols=F,clustering_distance_rows="correlation")

#only the top genes based on p values

top <- sig_hs_AvH %>%
  arrange(padj) %>%
  slice_head(n = 10)

#only upregulated
up <- sig_hs_AvH %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 20)



top2 <- heat[rownames(heat) %in% rownames(top),]

length(rownames(top2))
#10 :)

up2 <- heat[rownames(heat) %in% rownames(up),]

length(rownames(up))

#reorder columns based on oxygen category (see Carly's script)
#norm, hyp, an
top_f<-top2[,c(5,1,2,3,4,6,7,8)]

up_f<-up2[,c(5,1,2,3,4,6,7,8)]

#make heatmap
pheatmap(top_f, color=col, cluster_cols=F,clustering_distance_rows="correlation")

pheatmap(up_f, color=col, cluster_cols=F,clustering_distance_rows="correlation")



#save csv to easily remake heatmap later
write.csv(heat_sig_hs_gamma_AvH_f, 'hs_gamma_heat_AvH.csv')



#make report with annotations
#HvA
des2Report <- HTMLReport(shortName = 'HS prok gamma RNAseq_analysis_with_DESeq2_longer_namesHvA', title = 'HS prok gamma Differential expression HvA_names',reportDirectory = "./reports")
publish(dds3,des2Report, pvalueCutoff=0.1, factor = colData(dds3)$Oxy_cat, contrast =c("Oxy_cat", "Anoxic", "Hypoxic"), reportDir="./reports")
finish(des2Report)

des2Report <- HTMLReport(shortName = 'HS_prok_gamma_RNAseq_analysis_with_DESeq2_long_names_HvN_alldata', title = 'HS gamma Differential expression HvN',reportDirectory = "./reports")
publish(dds3,des2Report, pvalueCutoff=2, factor = colData(dds3)$Oxy_cat, contrast =c("Oxy_cat", "Normoxic", "Hypoxic"), reportDir="./reports")
finish(des2Report)

#try with specific contrasts 
#anoxia vs normoxia
#zero degs! 
des2Report <- HTMLReport(shortName = 'HS_prok_gamma_RNAseq_analysis_with_DESeq2_long_names_AvN', title = 'HS prok gamma Differential expression AvN',reportDirectory = "./reports")
publish(dds3,des2Report, pvalueCutoff=2, n=100, factor = colData(dds3)$Oxy_cat, reportDir="./reports")
finish(des2Report)
#made report just to see the direction of change

##kogs
#kog_mwu with all HS arch data based on log fold change, without no match genes 
#KOG directional file
df1 <-data.frame(res3_AvH)
df2 <- df1 

table(df2$padj < 0.1)
#avh LFC >0 -> upregulated in anoxia
df2$lfc <- df2$log2FoldChange

#make csv with just gene names and lfc
kog <- data.frame(rownames(df2))
kog$lfc <-df2$lfc 
names(kog)[1] <- 'gene_id'

#write as csv
write.csv(kog, file='hs_kog_gamma_HvA_lfc_final.csv', row.names = FALSE)

#kog HvN
df1 <-data.frame(res3_HvN)
df2 <- df1 

table(df2$padj < 0.1)
#LFC <0 -> upregulated in hypoxia
df2$lfc <- df2$log2FoldChange*(-1)

#make csv with just gene names and lfc
kog <- data.frame(rownames(df2))
kog$lfc <-df2$lfc 
names(kog)[1] <- 'gene_id'

#write as csv
write.csv(kog, file='hs_kog_gamma_HvN_lfc_final.csv', row.names = FALSE)

#AvN
df1 <-data.frame(res3_AvN)
df2 <- df1 

table(df2$padj < 0.1)
#LFC >0 -> upregulated in anoxia
df2$lfc <- df2$log2FoldChange

#make csv with just gene names and lfc
kog <- data.frame(rownames(df2))
kog$lfc <-df2$lfc 
names(kog)[1] <- 'gene_id'

#write as csv
write.csv(kog, file='hs_kog_gamma_AvN_lfc_final.csv', row.names = FALSE)

# Import KOG table Hymeraphia gamma

hs_gamma_kog2gene=read.table("hs_bin2_gammas_contig2kogClass_final.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)

###add gene names to kog class file to match

names(hs_gamma_kog2gene)[names(hs_gamma_kog2gene) == "V1"] <- "gene_id"

gn=read.table("hs_bin2_gamma_final_contig2geneName.tab",header=F,sep="	",quote="",row.names=NULL,stringsAsFactors=FALSE)

#gene_id
#gene_name

names(gn)[names(gn) == "V1"] <- "gene_id"

names(gn)[names(gn) == "V2"] <- "gene_name"

joined <- left_join(hs_gamma_kog2gene, gn, by = "gene_id")
joined$full <-paste(joined$gene_id,joined$gene_name)

#reorder so that combined column is first (also dropp extra columns not needed for count matrix)
#also got rid of names column
joined2 <-joined[, c(4, 2)]

#shorten beforehand
#joined2$full 
joined2$full <- gsub(" ", "_", joined2$full)
#replace all special characters
joined2$full <- gsub("[()]", "_", joined2$full)
joined2$full <- gsub("[:]", "_", joined2$full)
joined2$full <- gsub("[/]", "_", joined2$full)
joined2$full <- gsub("[.]", "_", joined2$full)
joined2$full <- gsub("[,]", "_", joined2$full)
joined2$full <- gsub("[[]", "_", joined2$full)
joined2$full <- gsub("[]]", "_", joined2$full)
#limit rownames to 50 characters
joined2$full <- strtrim(joined2$full, 50)

hs_gamma_ann_kog <- joined2


#lfc for comparison
hs_AvN.lfc <- read.csv('hs_kog_gamma_AvN_lfc_final.csv', header =TRUE)

hs_AvN_gamma.lth=kog.mwu(hs_AvN.lfc,hs_gamma_ann_kog) 
hs_AvN_gamma.lth 

#lfc (can switch them out to see the diffrences)
hs_HvN.lfc <- read.csv('hs_kog_gamma_HvN_lfc_final.csv', header =TRUE)
hs_HvN_gamma.lth=kog.mwu(hs_HvN.lfc,hs_gamma_ann_kog)
hs_HvN_gamma.lth

# hs HvA
hs_HvA.lfc <- read.csv('hs_kog_gamma_HvA_lfc_final.csv', header =TRUE)
hs_HvA_gamma.lth=kog.mwu(hs_HvA.lfc,hs_gamma_ann_kog)
hs_HvA_gamma.lth

#es_kogs

#need es kog names 

es_AvN_gamma.lfc <- read.csv('es_kog_gamma_AvN_lfc_final.csv', header =TRUE)

es_AvN_gamma.lth=kog.mwu(es_AvN_gamma.lfc,es_gamma_ann_kog) 
es_AvN_gamma.lth 

es_HvN_gamma.lfc <- read.csv('es_kog_gamma_HvN_lfc_final.csv', header =TRUE)

es_HvN_gamma.lth=kog.mwu(es_HvN_gamma.lfc,es_gamma_ann_kog) 
es_HvN_gamma.lth 

es_AvH_gamma.lfc <- read.csv('es_kog_gamma_HvA_lfc_final.csv', header =TRUE)

es_AvH_gamma.lth=kog.mwu(es_AvH_gamma.lfc,es_gamma_ann_kog) 
es_AvH_gamma.lth 


# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list('Es_Hypox_vs_Anox'=es_AvH_gamma.lth, "Es_Hypox"=es_HvN_gamma.lth, "Es_Anox"=es_AvN_gamma.lth, "Hs_Anox"=hs_AvN_gamma.lth,"Hs_Hypox"=hs_HvN_gamma.lth,"Hs_Hypox_vs_Anox"=hs_HvA_gamma.lth))

#hs missing RNA processing and modification COG...
#drop this from matrix 
matrix <- as.matrix(ktable)
matrix[-7,]
#plot KOG differences
pheatmap(matrix[-7,], clustering_distance_cols="correlation") 

#drop from each individual KOG table based on their row
ktable=makeDeltaRanksTable(list('Es_Hypox_vs_Anox'=es_AvH_gamma.lth[-7,], "Es_Hypox"=es_HvN_gamma.lth[-6,], "Es_Anox"=es_AvN_gamma.lth[-4,], "Hs_Anox"=hs_AvN_gamma.lth,"Hs_Hypox"=hs_HvN_gamma.lth,"Hs_Hypox_vs_Anox"=hs_HvA_gamma.lth))


pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor)
# p-values of these correlations in the upper panel:
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor.pval)
#only 18 KOGs (COGs) b/c proks don't have things like chromatin, etc

#save KOG csv files
write.csv(hs_HvA_gamma.lth, 'hs_gamma_KOG_HvA.csv')
write.csv(hs_HvN_gamma.lth, 'hs_gamma_KOG_HvN.csv')
write.csv(hs_AvN_gamma.lth, 'hs_gamma_KOG_AvN.csv')

#es (move up)
write.csv(es_AvH_gamma.lth, 'es_gamma_KOG_HvA.csv')
write.csv(es_HvN_gamma.lth, 'es_gamma_KOG_HvN.csv')
write.csv(es_AvN_gamma.lth, 'es_gamma_KOG_AvN.csv')


##bring in tables from files
hs_AvH_gamma.lth <- read.csv('hs_gamma_KOG_HvA.csv', header =TRUE)
hs_AvH_gamma.lth2 <- hs_AvH_gamma.lth %>% remove_rownames %>% column_to_rownames(var="X")
hs_AvN_gamma.lth <- read.csv('hs_gamma_KOG_AvN.csv', header =TRUE)
hs_AvN_gamma.lth2 <- hs_AvN_gamma.lth %>% remove_rownames %>% column_to_rownames(var="X")
hs_HvN_gamma.lth <- read.csv('hs_gamma_KOG_HvN.csv', header =TRUE)
hs_HvN_gamma.lth2 <- hs_HvN_gamma.lth %>% remove_rownames %>% column_to_rownames(var="X")

es_AvH_gamma.lth <- read.csv('es_gamma_KOG_HvA.csv', header =TRUE)
es_AvH_gamma.lth2 <- es_AvH_gamma.lth %>% remove_rownames %>% column_to_rownames(var="X")
es_AvN_gamma.lth <- read.csv('es_gamma_KOG_AvN.csv', header =TRUE)
es_AvN_gamma.lth2 <- es_AvN_gamma.lth %>% remove_rownames %>% column_to_rownames(var="X")
es_HvN_gamma.lth <- read.csv('es_gamma_KOG_HvN.csv', header =TRUE)
es_HvN_gamma.lth2 <- es_HvN_gamma.lth %>% remove_rownames %>% column_to_rownames(var="X")

ktable=makeDeltaRanksTable(list("Es_Anox"=es_AvN_gamma.lth2[-4,],
                                "Es_Hypox"=es_HvN_gamma.lth2[-6,],
                                'Es_Hypox_vs_Anox'=es_AvH_gamma.lth2[-7,],
                                "Hs_Anox"=hs_AvN_gamma.lth2,
                                "Hs_Hypox"=hs_HvN_gamma.lth2,
                                "Hs_Hypox_vs_Anox"=hs_AvH_gamma.lth2))


pheatmap(ktable, cluster_cols = F) 


##mitochondrial data
##########
#es mitos
##########
counts_esHost=read.table("allcounts_es_prok_mito_only.txt",header=TRUE,row.names=1) #Reading in the table of counts per isogroup by sample - you must open R in the same folder or change the working directory

counts_es=counts_esHost[,order(names(counts_esHost))]
names(counts_es)

##conditions (metadata)
conditions=read.table("Conditions_es.txt", header = T)

conds_es <- conditions[ order(row.names(conditions)), ]

####load data into deseq table format
#includes experimental design' since it is from the quick start method
DESeq2table_es <- DESeqDataSetFromMatrix(countData = counts_es,
                                         colData = conds_es,
                                         design= ~ Oxy_cat)

#make normoxic the 'reference' level 
DESeq2table_es$Oxy_cat <- relevel(DESeq2table_es$Oxy_cat, "Normoxic")

#runs differential expression anaylsis quickly
dds_es <- DESeq(DESeq2table_es)

resultsNames(dds_es) # lists the coefficients
res_es_AvN <- results(dds_es, name="Oxy_cat_Anoxic_vs_Normoxic")
# or to shrink log fold changes association with condition:

#shrunken results beased on LFC
res_es_AvN_apeglm <- lfcShrink(dds_es, coef="Oxy_cat_Anoxic_vs_Normoxic", type="apeglm")

#significantly expressed genes based on full results
table(res_es_AvN$padj < 0.1)

#FALSE  TRUE 
#14
#0 differentially expressed genes between normoxic and anoxic samples in mitos
#same as hs

#standard method - hypoxic vs normoxic
res_es_HvN <- results(dds_es, name="Oxy_cat_Hypoxic_vs_Normoxic")

table(res_es_HvN$padj < 0.1)

length(res_es_AvN[,1])

#14 genes in total for es

#comparison (contrast), between hypoxic and anoxic samples
res_es_HvA <- results(dds_es, contrast = list('Oxy_cat_Hypoxic_vs_Normoxic', 'Oxy_cat_Anoxic_vs_Normoxic'))

table(res_es_HvA$padj < 0.1)
#FALSE  TRUE 
#14

#0 genes differentially expressed between hypoxic and anoxic

#prep for PCA/heatmaps
#transformation 
vsd <- vst(dds_es, blind=FALSE, nsub=14)

df <- as.data.frame(colData(dds_es)[,c("Oxy_cat")])

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
plotPCA(vsd, intgroup=c('Oxy_cat'))

#lots of separation from one anoxic sample
#some separation with anoxia and normoxia..
#maybe would be better to use continuous oxygen rather than categorical


##After calling plotMA, one can use the function identify to interactively detect the row number of individual genes by clicking on the plot. One can then recover the gene identifiers by saving the resulting indices
#idx <- identify(res_es_AvN_apeglm$baseMean, res_es_AvN_apeglm$log2FoldChange)
##print gene names from points you clicked on 
#rownames(res_es)[idx]

####possible removing of outliers 
boxplot(log10(assays(dds_es)[["cooks"]]), range=0, las=2)

#anoxic samples both look like outliers, but biologically are not
#maybe DC110 is an outlier for mitos

#plotting dispersion estimates
plotDispEsts(dds_es)


#PCA to remove outlier
PCA <- plotPCA(vsd, intgroup=c('Oxy_cat'), returnData=TRUE)
PCA

#DC110 is the outlier based in the PCA, distorting the whole view, so we can try pulling it out
#not DC33 as in other sample sets
###########################
##removing 110
#############################
#because it removes a whole factor, have to change the design matrix or reload the data
DESeq2table_es2 <- DESeq2table_es
#drop normoxic as level 
#DESeq2table_es2$Oxy_cat <- droplevels.factor(DESeq2table_es2$Oxy_cat, "Normoxic")

#relevel with hypoxia as reference 
#DESeq2table_es2$Oxy_cat <- relevel(DESeq2table_es2$Oxy_cat, "Hypoxic")

DESeq2table_es_no_out2 <- DESeq2table_es2[, !(colnames(DESeq2table_es2) %in% c('DC110'))]
dds_es2 <- DESeq(DESeq2table_es_no_out2)

resultsNames(dds_es2) # lists the coefficients
res_es2_AvN <- results(dds_es2, name="Oxy_cat_Anoxic_vs_Normoxic")
# or to shrink log fold changes association with condition:

#shrunken results beased on LFC
res_es2_AvN_apeglm <- lfcShrink(dds_es2, coef="Oxy_cat_Anoxic_vs_Normoxic", type="apeglm")

#significantly expressed genes based on full results
table(res_es2_AvN$padj < 0.1)

#FALSE  TRUE 
#14
#0 differentially expressed genes between anoxic and normoxic samples, even with outlier removed

#significantly expressed genes based on LFC results
table(res_es2_AvN_apeglm$padj < 0.1)
##need to run the other multiple comparisons
#standard method - hypoxic vs normoxic
res_es2_HvN <- results(dds_es2, name="Oxy_cat_Hypoxic_vs_Normoxic")

table(res_es2_HvN$padj < 0.1)
#FALSE  TRUE 
#14
#0 differentially expressed genes between hypoxic and normoxic

# or to shrink log fold changes association with condition:
res_es2_HvN_apeglm <- lfcShrink(dds_es2, coef="Oxy_cat_Hypoxic_vs_Normoxic", type="apeglm")

table(res_es2_HvN_apeglm$padj < 0.1)
#same as above

length(res_es2_AvN[,1])

#14 genes in total for es

#comparison (contrast), between hypoxic and anoxic samples
res_es2_HvA <- results(dds_es2, contrast = list('Oxy_cat_Hypoxic_vs_Normoxic', 'Oxy_cat_Anoxic_vs_Normoxic'))

table(res_es2_HvA$padj < 0.1)
#FALSE  TRUE 
#13   1
#one gene significantly differentially expressed in hypoxia vs anoxia
#= nad2

res_es2_HvA_ashr <- lfcShrink(dds_es2, contrast = list("Oxy_cat_Anoxic_vs_Normoxic", "Oxy_cat_Hypoxic_vs_Normoxic"), type="ashr")

table(res_es2_HvA_ashr$padj < 0.1)

#prep for PCA/heatmaps
#transformation 
vsd <- vst(dds_es2, blind=FALSE, nsub=14)

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
plotPCA(vsd, intgroup=c('Oxy_cat'))

#no clear separation for mitos based on oxygen condition

#MA plots

plotMA(res_es2_AvN)

plotMA(res_es2_AvN_apeglm)

plotMA(res_es2_HvA_ashr)

plotMA(res_es2_HvN_apeglm)
##After calling plotMA, one can use the function identify to interactively detect the row number of individual genes by clicking on the plot. One can then recover the gene identifiers by saving the resulting indices
#idx <- identify(res_es_AvN_apeglm$baseMean, res_es_AvN_apeglm$log2FoldChange)
##print gene names from points you clicked on 
#rownames(res_es)[idx]

####possible removing of outliers 
boxplot(log10(assays(dds_es2)[["cooks"]]), range=0, las=2)

#anoxic samples still look like outliers, but we need to keep them in

#plotting dispersion estimates
plotDispEsts(dds_es2)

#plot based on minimum p value
plotCounts(dds_es2, gene=which.min(res_es2_AvN$padj), intgroup=c('Oxy_cat'))

#diagnostic histogram 
hist(res_es2_AvN$pvalue[res_es2_AvN$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")
########bokplot of all mito genes for es 
tcounts <- t(log2((counts(dds_es2, normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(dds_es2), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-14+1):ncol(.))

theme_set(theme_bw(base_size=14) + theme(strip.background = element_blank()))

es_mitos_gg <- ggplot(tcounts, aes(Oxy_cat, expression)) + 
  geom_boxplot() + 
  facet_wrap(~factor(gene, levels=c('atp6','atp8','atp9','cob', 'cox1', 'cox2', 'cox3', 'nad1', 'nad2', 'nad3', 'nad4', 'nad4l', 'nad5', 'nad6'))) +
  labs(x="Oxygen condition", 
       y="Expression (log normalized counts)") 
#       fill="(Some made \nup variable)", 
#       title="Top Results")
####
#hs
####
countsHost=read.table("allcounts_hs_prok_mitos_only.txt",header=TRUE,row.names=1) #Reading in the table of counts per isogroup by sample - you must open R in the same folder or change the working directory

counts=countsHost[,order(names(countsHost))]
names(counts)

##conditions (metadata)
conditions=read.table("Conditions_hs.txt", header = T)

conds <- conditions[ order(row.names(conditions)), ]

####load data into deseq table format
#includes experimental design' since it is from the quick start method
DESeq2Table <- DESeqDataSetFromMatrix(countData = countsHost,
                                      colData = conds,
                                      design= ~ Oxy_cat)

#make normoxic the 'reference' level 
DESeq2Table$Oxy_cat <- relevel(DESeq2Table$Oxy_cat, "Normoxic")

#runs differential expression anaylsis quickly
dds <- DESeq(DESeq2Table)

resultsNames(dds) # lists the coefficients
res_AvN <- results(dds, name="Oxy_cat_Anoxic_vs_Normoxic")
# or to shrink log fold changes association with condition:

#shrunken results beased on LFC
res_AvN_apeglm <- lfcShrink(dds, coef="Oxy_cat_Anoxic_vs_Normoxic", type="apeglm")

#significantly expressed genes based on full results
table(res_AvN$padj < 0.1)

#FALSE  TRUE 
#14
#no significantly differentially expressed genes

#significantly expressed genes based on LFC results
table(res_AvN_apeglm$padj < 0.1)

#standard method - hypoxic vs normoxic
res_HvN <- results(dds, name="Oxy_cat_Hypoxic_vs_Normoxic")

table(res_HvN$padj < 0.1)
#FALSE  TRUE 
#14

# or to shrink log fold changes association with condition:
res_HvN_apeglm <- lfcShrink(dds, coef="Oxy_cat_Hypoxic_vs_Normoxic", type="apeglm")

table(res_HvN_apeglm$padj < 0.1)
#same as above

length(res_HvN[,1])

#14 genes in total for HS (one must have an N.A. significance rather than a numerical value= only 13)

#comparison (contrast), between hypoxic and anoxic samples
res_HvA <- results(dds, contrast = list('Oxy_cat_Hypoxic_vs_Normoxic', 'Oxy_cat_Anoxic_vs_Normoxic'))

table(res_HvA$padj < 0.1)
#FALSE  TRUE 
#14

#shrunken results 
#apeglm does not have contrast ability
#need to try 'normal shrinkage, or ashr

res_HvA_ashr <- lfcShrink(dds, contrast = list("Oxy_cat_Anoxic_vs_Normoxic", "Oxy_cat_Hypoxic_vs_Normoxic"), type="ashr")

table(res_HvA_ashr$padj < 0.1)

#same number as above :)

#more about contrasts here: https://github.com/tavareshugo/tutorial_DESeq2_contrasts/blob/main/DESeq2_contrasts.md
#also includes instructions for combining levels using a matrix structure


#intercept does not work for normoxic vs. hypoxic... how do we get these data?!
#checked https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
#install.packages('pheatmap')

#prep for PCA/heatmaps
#transformation 
vsd <- vst(dds, blind=FALSE, nsub=14)

df <- as.data.frame(colData(dds)[,c("Oxy_cat")])

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
plotPCA(vsd, intgroup=c('Oxy_cat'))

#good separation of the anoxic sample from the hypoxic and normoxic ones
#unfortunately only one normoxic sample :(
#however, given that the sponges were behaving normally under hypoxia (pumping), hypoxic could be pretty close to normal
#so hypoxic dot on the left could be an outlier...

#could be effectively the same if top left outlier is removed from the PCA plot 

#MA plots

plotMA(res_AvN)

plotMA(res_AvN_apeglm)

plotMA(res_HvN_apeglm)

plotMA(res_HvA_ashr)

##After calling plotMA, one can use the function identify to interactively detect the row number of individual genes by clicking on the plot. One can then recover the gene identifiers by saving the resulting indices
#idx <- identify(res_AvN_apeglm$baseMean, res_AvN_apeglm$log2FoldChange)
##print gene names from points you clicked on 
#rownames(res)[idx]

####possible removing of outliers 
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)

#sample DC129 still consistantly higher outliers 
#this is the normoxic sample though... so it is probably best to include it.

#plotting dispersion estimates
plotDispEsts(dds)

###########################
##removing DC129
#############################
#because it removes a whole factor, have to change the design matrix or reload the data
DESeq2Table2 <- DESeq2Table
#drop normoxic as level 
DESeq2Table2$Oxy_cat <- droplevels.factor(DESeq2Table2$Oxy_cat, "Normoxic")

#relevel with hypoxia as reference 
DESeq2Table2$Oxy_cat <- relevel(DESeq2Table2$Oxy_cat, "Hypoxic")

DESeq2Table_no_out2 <- DESeq2Table2[, !(colnames(DESeq2Table2) %in% c('DC129'))]
dds2 <- DESeq(DESeq2Table_no_out2)

resultsNames(dds2) # lists the coefficients
res2_AvH <- results(dds2, name="Oxy_cat_Anoxic_vs_Hypoxic")
# or to shrink log fold changes association with condition:

#shrunken results beased on LFC
res2_AvH_apeglm <- lfcShrink(dds2, coef="Oxy_cat_Anoxic_vs_Hypoxic", type="apeglm")

#significantly expressed genes based on full results
table(res2_AvH$padj < 0.1)

#FALSE  TRUE 
#13
#0 differentially expressed genes between hypoxic and anoxic samples :)

#significantly expressed genes based on LFC results
table(res2_AvH_apeglm$padj < 0.1)

##same as above :)

#prep for PCA/heatmaps
#transformation 
vsd <- vst(dds2, blind=FALSE, nsub=14)

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
PCA <- plotPCA(vsd, intgroup=c('Oxy_cat'))
PCA

#could be effectively the same if top left outlier is removed from the PCA plot 

#MA plots

plotMA(res2_AvH)

plotMA(res2_AvH_apeglm)

##After calling plotMA, one can use the function identify to interactively detect the row number of individual genes by clicking on the plot. One can then recover the gene identifiers by saving the resulting indices
#idx <- identify(res_AvN_apeglm$baseMean, res_AvN_apeglm$log2FoldChange)
##print gene names from points you clicked on 
#rownames(res)[idx]

####possible removing of outliers 
boxplot(log10(assays(dds2)[["cooks"]]), range=0, las=2)

#this is the normoxic sample though... so it is probably best to include it.

#plotting dispersion estimates
plotDispEsts(dds2)

#plot based on minimum p value
plotCounts(dds2, gene=which.min(res2_AvH$padj), intgroup=c('Oxy_cat'))

#diagnostic histogram 
hist(res2_AvH$pvalue[res2_AvH$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")

#has same freguency with high on the left :) 
#see http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html

#make some nice plots with all the mito genes.
#see https://rpubs.com/turnersd/plot-deseq-results-multipage-pdf
tcounts <- t(log2((counts(dds, normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(dds), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-14+1):ncol(.))

theme_set(theme_bw(base_size=14) + theme(strip.background = element_blank()))

hs_mitos_gg <- ggplot(tcounts, aes(Oxy_cat, expression)) + 
  geom_boxplot() + 
  facet_wrap(~factor(gene, levels=c('atp6','atp8','atp9','cytb', 'cox1', 'cox2', 'cox3', 'nad1', 'nad2', 'nad3', 'nad4', 'nad4L', 'nad5', 'nad6'))) +
  labs(x="Oxygen condition", 
       y="Expression (log normalized counts)") 
#       fill="(Some made \nup variable)", 
#       title="Top Results")

#figure mitos
Figure_mitos <- ggarrange(es_mitos_gg, hs_mitos_gg, 
                            labels = c('A', 'B'),
                            common.legend = TRUE,
                            ncol = 2, nrow = 1, 
                            legend = 'bottom')
#sec





sessionInfo()


