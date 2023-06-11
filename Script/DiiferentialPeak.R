
rm(list = ls())
getwd()
setwd("/md01/Chenzh275/ori/Data/scATAC/cardiacmyocyte_regeneration/IDR/CellSubType/cardio/FindDAR/subtype-onlyD3-qval001/countreads-shuf/count")

library("DESeq")
library("RColorBrewer")
library("IDPmisc")
library(ComplexHeatmap)
library(pheatmap)
library(ggplot2)
library(preprocessCore)
library(circlize)

countsTable<-read.table("cardiamyocyte_subtype_per_sample_count.txt",header=T)
head(countsTable)
rownames(countsTable)<-countsTable$peak_name
countsTable <- countsTable[,5:9];
head(countsTable)
dim(countsTable)
colnames(countsTable) <- c("CM1","CM2", "CM3", "CM4", "CM5")
head(countsTable)

cds<-newCountDataSet(countData = countsTable,condition=as.factor(c("CM1","CM2","CM3","CM4","CM5")))
#dim(countsTable)
cds <- estimateSizeFactors(cds)
sizeFactors(cds)
cds <- estimateDispersions(cds, method="blind", sharingMode="fit-only",fitType="local" )

res.12 <- nbinomTest(cds,"CM1","CM2")
res.13 <- nbinomTest(cds,"CM1","CM3")
res.14 <- nbinomTest(cds,"CM1","CM4")
res.15 <- nbinomTest(cds,"CM1","CM5")
res.23 <- nbinomTest(cds,"CM2","CM3")
res.24 <- nbinomTest(cds,"CM2","CM4")
res.25 <- nbinomTest(cds,"CM2","CM5")
res.34 <- nbinomTest(cds,"CM3","CM4")
res.35 <- nbinomTest(cds,"CM3","CM5")
res.45 <- nbinomTest(cds,"CM4","CM5")

CM12<-as.data.frame(res.12)
CM13<-as.data.frame(res.13)
CM14<-as.data.frame(res.14)
CM15<-as.data.frame(res.15)
CM23<-as.data.frame(res.23)
CM24<-as.data.frame(res.24)
CM25<-as.data.frame(res.25)
CM34<-as.data.frame(res.34)
CM35<-as.data.frame(res.35)
CM45<-as.data.frame(res.45)

#仅通过pvalue和log2FC来挑选
CM12<-subset(CM12,abs(CM12$log2FoldChange)>1.5| abs(CM12$log2FoldChange) =="Inf")
CM13<-subset(CM13,abs(CM13$log2FoldChange)>1.5|  abs(CM13$log2FoldChange) =="Inf")
CM14<-subset(CM14,abs(CM14$log2FoldChange)>2|  abs(CM14$log2FoldChange) =="Inf")
CM15<-subset(CM15,abs(CM15$log2FoldChange)>1|  abs(CM15$log2FoldChange) =="Inf")
CM23<-subset(CM23,abs(CM23$log2FoldChange)>1|  abs(CM23$log2FoldChange) =="Inf")
CM24<-subset(CM24,abs(CM24$log2FoldChange)>1.2|  abs(CM24$log2FoldChange) =="Inf")
CM25<-subset(CM25,abs(CM25$log2FoldChange)>1|  abs(CM25$log2FoldChange) =="Inf")
CM34<-subset(CM34,abs(CM34$log2FoldChange)>1.85|  abs(CM34$log2FoldChange) =="Inf")
CM35<-subset(CM35,abs(CM35$log2FoldChange)>1|  abs(CM35$log2FoldChange) =="Inf")
CM45<-subset(CM45,abs(CM45$log2FoldChange)>2| abs(CM45$log2FoldChange) =="Inf")

CM12<-subset(CM12,log2(CM12$baseMean)>=5.3)
CM13<-subset(CM13,log2(CM13$baseMean)>=5.3)
CM14<-subset(CM14,log2(CM14$baseMean)>=4.8)
CM15<-subset(CM15,log2(CM15$baseMean)>=5.3)
CM23<-subset(CM23,log2(CM23$baseMean)>=5.5)
CM24<-subset(CM24,log2(CM24$baseMean)>=4.45)
CM25<-subset(CM25,log2(CM25$baseMean)>=5.8)
CM34<-subset(CM34,log2(CM34$baseMean)>=5.2)
CM35<-subset(CM35,log2(CM35$baseMean)>=4.95)
CM45<-subset(CM45,log2(CM45$baseMean)>=3.8)

CM14<-subset(CM14,log2(CM14$baseMean)>=4.8)
CM24<-subset(CM24,log2(CM24$baseMean)>=4.45)
CM34<-subset(CM34,log2(CM34$baseMean)>=5.2)
CM45<-subset(CM45,log2(CM45$baseMean)>=3.8)

CM12<-subset(CM12,CM12$pval <= 0.1)
CM13<-subset(CM13,CM13$pval <= 0.1)
CM14<-subset(CM14,CM14$pval <= 0.1)
CM15<-subset(CM15,CM15$pval <= 0.1)
CM23<-subset(CM23,CM23$pval <= 0.1)
CM24<-subset(CM24,CM24$pval <= 0.08)
CM25<-subset(CM25,CM25$pval <= 0.1)
CM34<-subset(CM34,CM34$pval <= 0.09)
CM35<-subset(CM35,CM35$pval <= 0.1)
CM45<-subset(CM45,CM45$pval <= 0.07)

CM12<-CM12$id
CM13<-CM13$id
CM14<-CM14$id
CM15<-CM15$id
CM23<-CM23$id
CM24<-CM24$id
CM25<-CM25$id
CM34<-CM34$id
CM35<-CM35$id
CM45<-CM45$id

filter<-union(CM12,CM13)
filter<-union(filter,CM14)
filter<-union(filter,CM15)
filter<-union(filter,CM23)
filter<-union(filter,CM24)
filter<-union(filter,CM25)
filter<-union(filter,CM34)
filter<-union(filter,CM35)
filter<-union(filter,CM45)
length(filter) #17035
#15391

filtercounts<-countsTable[which(rownames(countsTable)%in%filter),]
head(filtercounts)
dim(filtercounts)
filtermat=t(scale(t(filtercounts),scale = T,center = T))
#boxplot(log2(filtermat),main="Nor")
filtermat <- normalize.quantiles(filtermat)
rownames(filtermat) <- rownames(filtercounts)
colnames(filtermat) <- c("CM1","CM2","CM3","CM4","CM5")
#boxplot(log2(filtermat),main="Nor")
dim(filtermat)
#head(filtermat)

set.seed(123)
list<-draw(Heatmap(filtermat,km = 10
                   , show_row_names = F, cluster_columns = FALSE, show_heatmap_legend = T 
                   , col = c("blue","white","red"), # show_parent_dend_line = FALSE
                   , column_title="Cardiomyocyte Subtype Chromatin Accessibility"))

cluster1 = filtermat[row_order(list)[[6]],] 
head(cluster1)
cluster2 = filtermat[row_order(list)[[7]],] 
dim(cluster2)
cluster3 = filtermat[row_order(list)[[5]],] 
cluster4 = filtermat[row_order(list)[[4]],] 
cluster5 = filtermat[row_order(list)[[2]],]  
cluster6 = filtermat[row_order(list)[[3]],]
cluster7 = filtermat[row_order(list)[[1]],]

cluster4 = filtermat[row_order(list)[[1]],]
cluster5 = filtermat[row_order(list)[[4]],]
cluster3 = filtermat[row_order(list)[[4]],]
cluster6 = filtermat[row_order(list)[[3]],]
cluster7 = filtermat[row_order(list)[[3]],]

cluster6 = filtermat[row_order(list)[[6]],]
cluster8 = filtermat[row_order(list)[[8]],]
cluster9 = filtermat[row_order(list)[[10]],]
cluster10 = filtermat[row_order(list)[[10]],]

coordinate<-read.table("cardiamyocyte_subtype_per_sample_count.txt",header=T)
coordinate<-coordinate[,1:4]
dim(coordinate)
head(coordinate)
recluster1<-coordinate[which(coordinate$peak_name%in%rownames(cluster1)),]
head(recluster1)
recluster2<-coordinate[which(coordinate$peak_name%in%rownames(cluster2)),]
recluster3<-coordinate[which(coordinate$peak_name%in%rownames(cluster3)),]
recluster4<-coordinate[which(coordinate$peak_name%in%rownames(cluster4)),]
recluster5<-coordinate[which(coordinate$peak_name%in%rownames(cluster5)),]
dim(recluster5)
recluster6<-coordinate[which(coordinate$peak_name%in%rownames(cluster6)),]
recluster7<-coordinate[which(coordinate$peak_name%in%rownames(cluster7)),]
dim(recluster7)

recluster8<-coordinate[which(coordinate$peak_name%in%rownames(cluster8)),]
recluster9<-coordinate[which(coordinate$peak_name%in%rownames(cluster9)),]
recluster10<-coordinate[which(coordinate$peak_name%in%rownames(cluster10)),]

total_peak <- coordinate[which(coordinate$peak_name%in%rownames(filtermat)),]
dim(total_peak)

library(clusterProfiler)
library("org.Mm.eg.db")
library("ChIPseeker")
library(AnnotationDbi)
library(GenomicFeatures)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

library(GenomicRanges)
recluster3_gene <- makeGRangesFromDataFrame(recluster3,keep.extra.columns=TRUE)
recluster4_gene <- makeGRangesFromDataFrame(recluster4,keep.extra.columns=TRUE)
recluster5_gene <- makeGRangesFromDataFrame(recluster5,keep.extra.columns=TRUE)
recluster6_gene <- makeGRangesFromDataFrame(recluster6,keep.extra.columns=TRUE)
recluster7_gene <- makeGRangesFromDataFrame(recluster7,keep.extra.columns=TRUE)
recluster8_gene <- makeGRangesFromDataFrame(recluster8,keep.extra.columns=TRUE)
recluster9_gene <- makeGRangesFromDataFrame(recluster9,keep.extra.columns=TRUE)
recluster10_gene <- makeGRangesFromDataFrame(recluster10,keep.extra.columns=TRUE)
total_peak_gene <- makeGRangesFromDataFrame(total_peak,keep.extra.columns=TRUE)

markerpeak_gene <- annotatePeak(recluster6_gene, TxDb=txdb,tssRegion=c(-3000, 3000), 
                                addFlankGeneInfo=TRUE, flankDistance=5000, annoDb="org.Mm.eg.db")
markerpeak_gene
markerpeak_gene <- data.frame(markerpeak_gene)

GO_gene <- markerpeak_gene$geneId

GO_geneID <- bitr(GO_gene , fromType = "ENTREZID", toType = c("ENSEMBL", "SYMBOL"), 
  OrgDb = org.Mm.eg.db, drop = T)

select_gene <- c("Erbb4","Esrrg","Gata4","Mef2a","Mef2c","Tead3","E2f3","Zfpm2")

slect <- GO_geneID[which(GO_geneID$SYMBOL%in%select_gene),]
slect

peak_selct <- markerpeak_gene[which(markerpeak_gene$SYMBOL %in%select_gene ),]
table(peak_selct$SYMBOL)
peak_selct_uniq <- markerpeak_gene[!duplicated(markerpeak_gene$SYMBOL),]
dim(peak_selct_uniq)

 E2f3 Erbb4 Esrrg Mef2a Zfpm2 10
Gata4 tead3 Mef2c 
peak_selct_uniq10 <- peak_selct_uniq$SYMBOL
peak_selct_uniq6 <- peak_selct_uniq$SYMBOL
peak_selct_uniq <- c(peak_selct_uniq10,peak_selct_uniq6)

peak_selct_uniq2 <- peak_selct_uniq[which(peak_selct_uniq %in%select_gene )]
peak_selct_uniq3 <- peak_selct_uniq[!duplicated(peak_selct_uniq)]
length(peak_selct_uniq3)

write.csv(peak_selct_uniq[,c(1:3,17)],"CM4_cluster4_gene.csv",quote=F)

peak_selct <- markerpeak_gene[which(markerpeak_gene$SYMBOL =="Tead3" ),]
peak <- peak_selct$peak_name

filtermat_select <- filtermat[which(rownames(filtermat)%in% peak),]
filtermat_select

