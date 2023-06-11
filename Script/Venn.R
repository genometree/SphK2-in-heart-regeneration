rm(list = ls())
getwd()
setwd("C:/Users/NoName/OneDrive - ÖÐÉ½´óÑ§/Project/Cardiac regeneration/CM/Analysis-Day3-only/Countreads/")

#GO term 
library(clusterProfiler)
library("org.Mm.eg.db")
library("ChIPseeker")
library(AnnotationDbi)
library(GenomicFeatures)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

CM1_peak <- readPeakFile("./shuf/CM1-homer.bed")
CM2_peak <- readPeakFile("./shuf/CM2-homer.bed")
CM3_peak <- readPeakFile("./shuf/CM3-homer.bed")
CM4_peak <- readPeakFile("./shuf/CM4-homer.bed")
CM4_2_peak <- readPeakFile("./shuf/CM4.2-homer.bed")
CM5_peak <- readPeakFile("./shuf/CM5-homer.bed")

##########CM1##########
CM1_peakAnnoList<- annotatePeak(CM1_peak, TxDb=txdb,tssRegion=c(-1000, 1000), 
                                addFlankGeneInfo=TRUE, flankDistance=2000, annoDb="org.Mm.eg.db")
CM1_peakAnnoList
CM1_gene = as.data.frame(CM1_peakAnnoList)$geneId
CM1_geneID <- bitr(CM1_gene, fromType = "ENTREZID", toType = c("ENSEMBL", "SYMBOL"), OrgDb = org.Mm.eg.db,drop = T)
CM1_ego <- enrichGO(CM1_geneID$ENTREZID,
                    keyType = 'ENTREZID',
                    OrgDb = org.Mm.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.01,
                    qvalueCutoff = 0.01,
                    readable = TRUE)
barplot(CM1_ego, showCategory = 15,font.size = 15 )
dotplot(CM1_ego, showCategory = 15,orderBy = "pvalue",font.size = 15)
CM1_ego <- as.data.frame(CM1_ego)
CM1_ego$Description
write.table(CM1_ego,"./shuf/CM1_ego.txt",sep="\t",quote = FALSE, row.names=F)

##########CM2##########
CM2_peakAnnoList<- annotatePeak(CM2_peak, TxDb=txdb,tssRegion=c(-1000, 1000), 
                                addFlankGeneInfo=TRUE, flankDistance=2000, annoDb="org.Mm.eg.db")
CM2_peakAnnoList
plotAnnoPie(CM2_peakAnnoList)
CM2_gene = as.data.frame(CM2_peakAnnoList)$geneId
CM2_geneID <- bitr(CM2_gene, fromType = "ENTREZID", toType = c("ENSEMBL", "SYMBOL"), OrgDb = org.Mm.eg.db,drop = T)
CM2_ego <- enrichGO(CM2_geneID$ENTREZID,
                    keyType = 'ENTREZID',
                    OrgDb = org.Mm.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.01,
                    qvalueCutoff = 0.01,
                    readable = TRUE)
barplot(CM2_ego, showCategory = 15 ,font.size = 15 )
dotplot(CM2_ego, showCategory = 15,orderBy = "pvalue",font.size = 15 )
CM2_ego <- as.data.frame(CM2_ego)
CM2_ego$Description
write.table(CM2_ego,"./shuf/CM2_ego.txt",sep="\t",quote = FALSE, row.names=F)

##########CM3##########
CM3_peakAnnoList<- annotatePeak(CM3_peak, TxDb=txdb,tssRegion=c(-1000, 1000), 
                                addFlankGeneInfo=TRUE, flankDistance=2000, annoDb="org.Mm.eg.db")
CM3_peakAnnoList
plotAnnoPie(CM3_peakAnnoList)
CM3_gene = as.data.frame(CM3_peakAnnoList)$geneId
CM3_geneID <- bitr(CM3_gene, fromType = "ENTREZID", toType = c("ENSEMBL", "SYMBOL"), OrgDb = org.Mm.eg.db,drop = T)
CM3_ego <- enrichGO(CM3_geneID$ENTREZID,
                    keyType = 'ENTREZID',
                    OrgDb = org.Mm.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.01,
                    qvalueCutoff = 0.01,
                    readable = TRUE)
barplot(CM3_ego, showCategory = 15 ,font.size = 15  )
dotplot(CM3_ego, showCategory = 15,orderBy = "pvalue",font.size = 15 )
CM3_ego <- as.data.frame(CM3_ego)
CM3_ego$Description
write.table(CM3_ego,"./shuf/CM3_ego.txt",sep="\t",quote = FALSE, row.names=F)

##########CM4##########
CM4_peakAnnoList<- annotatePeak(CM4_peak, TxDb=txdb,tssRegion=c(-1000, 1000), 
                                addFlankGeneInfo=TRUE, flankDistance=2000, annoDb="org.Mm.eg.db")
CM4_peakAnnoList
plotAnnoPie(CM4_peakAnnoList)
CM4_gene = as.data.frame(CM4_peakAnnoList)$geneId
CM4_geneID <- bitr(CM4_gene, fromType = "ENTREZID", toType = c("ENSEMBL", "SYMBOL"), OrgDb = org.Mm.eg.db,drop = T)
CM4_ego <- enrichGO(CM4_geneID$ENTREZID,
                    keyType = 'ENTREZID',
                    OrgDb = org.Mm.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.01,
                    qvalueCutoff = 0.01,
                    readable = TRUE)
barplot(CM4_ego, showCategory = 15 ,font.size = 15 )
dotplot(CM4_ego, showCategory = 15,orderBy = "pvalue",font.size = 15 )
CM4_ego <- as.data.frame(CM4_ego)
CM4_ego$Description
write.table(CM4_ego,"./shuf/CM4_ego.txt",sep="\t",quote = FALSE, row.names=F)

##########CM4-2##########
CM4_2_peakAnnoList<- annotatePeak(CM4_2_peak, TxDb=txdb,tssRegion=c(-1000, 1000), 
                                addFlankGeneInfo=TRUE, flankDistance=2000, annoDb="org.Mm.eg.db")
CM4_2_peakAnnoList
plotAnnoPie(CM4_2_peakAnnoList)
CM4_2_gene = as.data.frame(CM4_2_peakAnnoList)$geneId
CM4_2_geneID <- bitr(CM4_2_gene, fromType = "ENTREZID", toType = c("ENSEMBL", "SYMBOL"), OrgDb = org.Mm.eg.db,drop = T)
CM4_2_ego <- enrichGO(CM4_2_geneID$ENTREZID,
                    keyType = 'ENTREZID',
                    OrgDb = org.Mm.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.01,
                    qvalueCutoff = 0.01,
                    readable = TRUE)
barplot(CM4_2_ego, showCategory = 15 ,font.size = 15 )
dotplot(CM4_2_ego, showCategory = 15,orderBy = "pvalue",font.size = 15 )
CM4_2_ego <- as.data.frame(CM4_2_ego)
CM4_2_ego$Description
write.table(CM4_2_ego,"./shuf/CM4_2_ego.txt",sep="\t",quote = FALSE, row.names=F)

#######merge CM4 and CM4¡ª¡ª2##########
CM4_peak_2 <- read.table("./shuf/CM4-homer.bed")
CM4_2_peak_2 <- read.table("./shuf/CM4.2-homer.bed")
All_CM4_peak <-rbind(CM4_peak_2,CM4_2_peak_2)
head(All_CM4_peak)
colnames(All_CM4_peak) =c("Chr","Start","End","Peak","Strand")
All_CM4_peak <- GRanges(All_CM4_peak)

All_CM4_peakAnnoList<- annotatePeak(All_CM4_peak, TxDb=txdb,tssRegion=c(-1000, 1000), 
                                  addFlankGeneInfo=TRUE, flankDistance=2000, annoDb="org.Mm.eg.db")
All_CM4_peakAnnoList
plotAnnoPie(All_CM4_peakAnnoList)
All_CM4_gene = as.data.frame(All_CM4_peakAnnoList)$geneId
All_CM4_geneID <- bitr(All_CM4_gene, fromType = "ENTREZID", toType = c("ENSEMBL", "SYMBOL"), OrgDb = org.Mm.eg.db,drop = T)
All_CM4_ego <- enrichGO(All_CM4_geneID$ENTREZID,
                      keyType = 'ENTREZID',
                      OrgDb = org.Mm.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.01,
                      readable = TRUE)
barplot(All_CM4_ego, showCategory = 15 ,font.size = 15 )
dotplot(All_CM4_ego, showCategory = 15,orderBy = "pvalue",font.size = 15 )
All_CM4_ego <- as.data.frame(All_CM4_ego)
All_CM4_ego$Description
write.table(All_CM4_ego,"./shuf/All_CM4_ego.txt",sep="\t",quote = FALSE, row.names=F)

#######CM4 artifitial select genes##########
genes <- read.table("./shuf/CM4-artifitial-selectGenens.txt")
tail(genes)
CM4_arti_gene = as.data.frame(genes)$V1
CM4_arti_geneID <- bitr(CM4_arti_gene, fromType = "SYMBOL", toType = c("ENSEMBL","ENTREZID"), OrgDb = org.Mm.eg.db,drop = T)
CM4_arti_ego <- enrichGO(CM4_arti_geneID$ENTREZID,
                    keyType = 'ENTREZID',
                    OrgDb = org.Mm.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.01,
                    qvalueCutoff = 0.01,
                    readable = TRUE)
barplot(CM4_arti_ego, showCategory = 15 ,font.size = 15  )
dotplot(CM4_arti_ego, showCategory = 15,orderBy = "pvalue",font.size = 15 )
CM4_arti_ego <- as.data.frame(CM4_arti_ego)
CM4_arti_ego$Description
write.table(CM4_arti_ego,"./shuf/CM4_arti_ego.txt",sep="\t",quote = FALSE, row.names=F)

data =read.csv("./shuf/CM4-artifitial-GOterm.csv",header = T)
head(data)
data<- edit(data)
my_levels <-data$Description
data$Description <- factor(data$Description, levels= c("muscle tissue development"
                                                       ,"striated muscle tissue development"
                                                       ,"muscle cell differentiation"
                                                       ,"striated muscle cell differentiation"
                                                       ,"heart process"
                                                       ,"cardiac muscle contraction"
                                                       , "cardiac cell development"
                                                       ,"cardiac muscle cell differentiation"
                                                       ,"axonogenesis"
                                                       ,"mitotic cell cycle phase transition"
                                                       ,"glycogen metabolic process"
                                                       ,"cell cycle G1/S phase transition"
                                                       ,"cell cycle phase transition"
                                                       ,"regulation of mitotic cell cycle phase transition"
                                                       ,"positive regulation of cell cycle" ))

p<-ggplot(data, aes(x = GeneRatio, y = Description , size = Count, color= p.adjust)) + geom_point()+
  scale_size_continuous(range=c(5,10))+
  theme(axis.text = element_text(face="bold", color="black", size=15),axis.text.x = element_text(angle = 60, hjust = 1),
        panel.border = element_rect(colour="black",fill=NA)) +
  scale_color_gradient(low="#f42d23",high = "#441383") 
p

##########CM5##########
CM5_peakAnnoList<- annotatePeak(CM5_peak, TxDb=txdb,tssRegion=c(-1000, 1000), 
                                addFlankGeneInfo=TRUE, flankDistance=2000, annoDb="org.Mm.eg.db")
CM5_peakAnnoList
plotAnnoPie(CM5_peakAnnoList)
CM5_gene = as.data.frame(CM5_peakAnnoList)$geneId
CM5_geneID <- bitr(CM5_gene, fromType = "ENTREZID", toType = c("ENSEMBL", "SYMBOL"), OrgDb = org.Mm.eg.db,drop = T)
CM5_ego <- enrichGO(CM5_geneID$ENTREZID,
                    keyType = 'ENTREZID',
                    OrgDb = org.Mm.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.01,
                    qvalueCutoff = 0.01,
                    readable = TRUE)
barplot(CM5_ego, showCategory = 15 ,font.size = 15 )
dotplot(CM5_ego, showCategory = 15,orderBy = "pvalue",font.size = 15 )
CM5_ego <- as.data.frame(CM5_ego)
CM5_ego$Description
write.table(CM5_ego,"./shuf/CM5_ego.txt",sep="\t",quote = FALSE, row.names=F)
