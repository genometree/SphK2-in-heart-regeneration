rm(list = ls())
getwd()
setwd("/md01/Chenzh275/ori/Data/scATAC/cardiacmyocyte_regeneration/IDR/CellSubType/cardio")

library(SnapATAC);
library(ggplot2);
library(GenomicRanges);
library(viridisLite);

rm(list = ls())
getwd()
setwd("/public/home/Chenzh275/Data/scATAC/cardiacmyocyte_regeneration/IDR/CellSubType/cardio")

library(SnapATAC);
library(ggplot2);
library(GenomicRanges);
library(viridisLite);

x.sp = readRDS("/md01/nieyg/scATAC/all/all.rds");
table(x.sp@sample)
x.sp=x.sp[which(x.sp@sample== "AR3" | x.sp@sample== "AR7" |x.sp@sample== "NAR3" |x.sp@sample== "NAR7" |
                  x.sp@sample== "P4" |x.sp@sample== "P8" ),]

x.sp = runViz(
  obj=x.sp, 
  tmp.folder=tempdir(),
  dims=2,
  eigs.dims=1:10, 
  method="umap",
  seed.use=18
);

###celltype annotation by umap
genes = read.table("/md01/Chenzh275/ori/Data/scATAC/cardiacmyocyte_regeneration/IDR/CellSubType/cardio/gene_site/gencode.vM16.gene.bed");
genes.gr = GRanges(genes[,1], IRanges(genes[,2], genes[,3]), name=genes[,4]  );
marker.genes = c(
  "Ckap4",
  "Fbln2",
  "Col3a1",
  "Mmp2",
  "Col1a2",
  "Fstl1",
  "Gsn",
  "Sparc",
  "Vim",
  "Thy1",
  "Pdgfra"
 );

marker.genes <- c(
  "Tnnt2","Myl2","Actc1","Nppa","Myh6","Atp2a2","Acta1", ####cardio
  "Fcgr1","Adgre1","Cd14","Csf1r","Cd163","Cd68",
  "Itgam","Lgals3","Lyzl1","Mrc1","Fabp5","Mertk", ###Macrophages
  "Cd3e","Cd3d","Cd8a","Cd8b1","Nkg7","Igfbp4","Lat", ###Tcell
  "Cd79a","Cd79b","Mzb1","Ly6d", ###Bcell
  "Cd74","Cd83","Cd86","Flt3","Cd209a", ####DC
  "Cdh5","Pecam1","Ednrb","Aqp7","Emcn","Madcam1","Epas1","Flt1",
  "Tie1","Fabp4","Esam","Kdr","Tek","Eng",###endothelial
  "Col3a1","Mmp2","Col1a2","Fstl1","Gsn","Thy1","Pdgfra","Lama2",
  "Dclk1","Fkbp5","Akt3",##fibroblast
  "Rgs5","Ano1","Acta2"####smoothmuscle
)

marker.genes <- c("Rgs5","Abcc9")  ####smoothmuscle

genes.sel.gr <- genes.gr[which(genes.gr$name %in% marker.genes)];
genes.sel.gr
#write.table(genes.sel.gr,"gegeggggggg.txt")
x.sp = createGmatFromMat(
    obj=x.sp, 
    input.mat="bmat",
    genes=genes.sel.gr,
    do.par=TRUE,
    num.cores=5
  ); 

x.sp = scaleCountMatrix(
    obj=x.sp, 
    cov=x.sp@metaData$passed_filters + 1,
    mat="gmat",
    method = "RPM"
  );
x.sp = runMagic(
    obj=x.sp,
    input.mat="gmat",
    step.size=3
  );
x.sp
#pdf("Fibroblast-marker-gene-annotation-in-umap.pdf");
pdf("allcetype-marker-gene-annotation-in-umap-ver2.pdf")
pdf("SMC-marker-gene-annotation-in-umap.pdf")
#par(mfrow = c(2,2));
for(i in 1:2){
    plotFeatureSingle(
        obj=x.sp,
        feature.value=x.sp@gmat[, marker.genes[i]],
        method="umap", 
        main=marker.genes[i],
        point.size=0.4, 
        point.shape=19, 
        down.sample=10000,
        quantiles=c(0, 1)
  )};
dev.off()

#total = readRDS("cardio.rds")
x.sp = readRDS("/md01/Chenzh275/ori/Data/scATAC/cardiacmyocyte_regeneration/IDR/CellSubType/cardio/cardiomyocyte_addpmat_allinfo2.rds");
x.sp=x.sp[which(x.sp@sample== "AR3" |x.sp@sample== "NAR3"| x.sp@sample== "P4" ),]

#########annotate cluster into subtype (without seperate cluster9)####
meta<-matrix(data=NA,nrow=1925,ncol=2);
rownames(meta) = x.sp@barcode
colnames(meta)=c("cluster","subtype");
meta[,1] = x.sp@cluster
meta <- as.data.frame(meta)
cl <- meta$cluster
cl <- as.numeric(cl)
for ( i in 1:length(cl))
{
  if (cl[i] == 3 |cl[i] == 5 | cl[i] == 6 |cl[i] ==  7 | cl[i] == 12  )
    meta[i,2] = "CM1"
  else if (cl[i] == 1)
    meta[i,2]  = "CM2"
  else if (cl[i] == 2)
    meta[i,2]  = "CM3"
  else if (cl[i] == 8 | cl[i] == 9 | cl[i] == 10 | cl[i] == 11 | cl[i] == 13 )
    meta[i,2]  = "CM4"
  else
    meta[i,2]  = "CM5"
}
head(meta)
x.sp@metaData$subtype = meta$subtype

#####annotate cluster into subtype (seperate cluster9)####
meta<-matrix(data=NA,nrow=1925,ncol=3);
rownames(meta) = x.sp@barcode
colnames(meta)=c("cluster","sample","subtype");
meta[,1] = x.sp@cluster
meta[,2] = x.sp@sample
for ( i in 1:nrow(meta))
{
  if (meta[i,1] == 1)
    meta[i,3]  = "CM2"
  else if (meta[i,1] == 2)
    meta[i,3]  = "CM3"
  else if (meta[i,1] == 8  | meta[i,1] == 10 |meta[i,1] == 11 | meta[i,1] == 13 ) 
    meta[i,3]  = "CM4"
  else if (meta[i,1] == 4)
    meta[i,3] = "CM5"
  else if (meta[i,1] == 9 & meta[i,2]== "AR3")
    meta[i,3] = "CM4"
  else
    meta[i,3]  = "CM1"
}
meta <- as.data.frame(meta)
head(meta)
table(meta$subtype)
table(meta$sample)
table(x.sp@metaData$subtype2)
x.sp@metaData$subtype2 = meta$subtype

#saveRDS(x.sp,"cardiomyocyte_addpmat_allinfo2.rds")

####umap
df = data.frame(x.sp@umap, x.sp@metaData$timepoint)
df = data.frame(x.sp@umap,x.sp@metaData$subtype2)
df = data.frame(x.sp@umap,x.sp@metaData$Sample2)
Timepoint = x.sp@metaData$timepoint
Subtype = x.sp@metaData$subtype2
Sample = x.sp@metaData$Sample2
###timepoint
ggplot(df, aes(x = umap.1, y = umap.2, color = Timepoint)) +
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = c("#FFA500","#0073C2FF", "#EFC000FF", "#868686FF")) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("Timepoint")

##subtype umap
ggplot(df, aes(x = umap.1, y = umap.2, color = Subtype)) +
  geom_point(size=0.8) + 
  theme_bw() + 
  scale_color_manual(values = c("#BC3C29FF","#0072B5FF", "#E18727FF","#20854EFF","#7876B1FF")) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("Subtype")

###sample
library(ggsci)
ggplot(df, aes(x = umap.1, y = umap.2, color = Sample)) +
  geom_point(alpha=0.7) + 
  theme_bw() + 
  scale_color_manual(values = c("#E64B35B2","#00468BB2","#42B540B2" )) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("Sample")


table(x.sp@sample)
table(x.sp@cluster)

x.sp = runViz(
  obj=x.sp,
  tmp.folder=tempdir(),
  dims=2,
  eigs.dims=1:10,
  method="umap",
  seed.use=18
)

x.sp=x.sp[which(x.sp@sample== "AR3" |x.sp@sample== "NAR3" | x.sp@sample== "P4" ),]

####umap
df = data.frame(x.sp@umap,x.sp@metaData$subtype)
Subtype = x.sp@metaData$subtype

##subtype umap
pdf("Subtype-umap-only-D3.pdf")
ggplot(df, aes(x = umap.1, y = umap.2, color = Subtype)) +
  geom_point(size=0.8) + 
  theme_bw() + 
  scale_color_manual(values = c("#BC3C29FF","#0072B5FF", "#E18727FF","#20854EFF","#7876B1FF")) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("Subtype")
dev.off()

####umap
df = data.frame(x.sp@umap,x.sp@metaData$sample)
Sample = x.sp@metaData$sample

library(ggsci)
pdf("Sample-umap-only-D3.pdf")
ggplot(df, aes(x = umap.1, y = umap.2, color = Sample)) +
  geom_point(alpha=0.7) + 
  theme_bw() + 
  scale_color_manual(values = c("#E64B35B2","#00468BB2","#42B540B2" )) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("Sample")
dev.off()

#total cell
total = readRDS("/md01/nieyg/scATAC/all/all.rds");
table(total@sample)
total=total[which(total@sample== "AR3" |total@sample== "NAR3" |  total@sample== "P4" ),]
table(total@sample)

total = runViz(
  obj=total, 
  tmp.folder=tempdir(),
  dims=2,
  eigs.dims=1:10, 
  method="umap",
  seed.use=18
)

markergene <- c(
  "Tnnt2","Myl2","Actc1","Nppa","Myh6","Atp2a2","Acta1", ####cardio
  "Fcgr1","Adgre1","Cd14","Csf1r","Cd163","Cd68",
  "Itgam","Lgals3","Lyzl1","Mrc1","Fabp5","Mertk", ###Macrophages
  "Cd3e","Cd3d","Cd8a","Cd8b1","Nkg7","Igfbp4","Lat", ###Tcell
  "Cd79a","Cd79b","Mzb1","Ly6d", ###Bcell
  "Cd74","Cd83","Cd86","Flt3","Cd209a", ####DC
  "Cdh5","Pecam1","Ednrb","Aqp7","Emcn","Madcam1","Epas1","Flt1",
  "Tie1","Fabp4","Esam","Kdr","Tek","Eng",###endothelial
  "Col3a1","Mmp2","Col1a2","Fstl1","Gsn","Thy1","Pdgfra","Lama2",
  "Dclk1","Fkbp5","Akt3",##fibroblast
  "Rgs5","Ano1","Acta2"####smoothmuscle
)

markergene <- c("Cd74","Cd22","Cd79b","Cxcr4", #B cell
	"Tnnt2","Myl2","Actc1","Nppa","Myh6","Atp2a2","Acta1", ####cardio
  "Cdh5","Pecam1","Ednrb","Aqp7","Emcn","Madcam1","Epas1","Flt1",
  "Tie1","Fabp4","Esam","Kdr","Tek","Eng","Adgrf5","Cldn5",###EC
	"Col3a1","Col6a3","Col6a1","Ddr2","Col5a1", "Pdgfra", #FB
	"Samsn1","Clec4d","S100a9", "S100a8",#Granulocy
	"Cd86","Mpeg1","Cd68","Ccl3","Cd44",# MP
	"Rgs5","Pdgfrb","Accsl","Abcc9","Rgs4", #SMC
	"Cd3e","Cd3d", "Cd8a", "Cd8b1","nkg7","Igfbp4","Lat","Cd53","Itk" #T cell
	)

genes = read.table("/md01/Chenzh275/ori/Data/scATAC/cardiacmyocyte_regeneration/IDR/CellSubType/cardio/gene_site/gencode.vM16.gene.bed");
genes.gr = GRanges(genes[,1], IRanges(genes[,2], genes[,3]), name=genes[,4]  );

genes.sel.gr <- genes.gr[which(genes.gr$name %in% markergene)];
genes.sel.gr

markergene <- markergene[which(markergene %in% genes.sel.gr$name)]
total = createGmatFromMat(
    obj=total, 
    input.mat="bmat",
    genes=genes.sel.gr,
    do.par=TRUE,
    num.cores=5
  )

total = scaleCountMatrix(
    obj=total, 
    cov=total@metaData$passed_filters + 1,
    mat="gmat",
    method = "RPM"
  )

total = runMagic(
    obj=total,
    input.mat="gmat",
    step.size=3
  )

total

pdf("allcetype-marker-gene-annotation-in-umap-splitIC.pdf")
for(i in 1:55){
    plotFeatureSingle(
        obj=total,
        feature.value=total@gmat[,markergene[i]],
        method="umap", 
        main=markergene[i],
        point.size=0.4, 
        point.shape=19, 
        down.sample=10000,
#        quantiles=c(0, 1)
  )
}
dev.off()

pdf("allcetype-cluster-in-umap-splitIC.pdf")
plotViz(
  obj=total,
  method="umap", 
  main=" Cluster",
  point.color=total@cluster, 
  point.size=0.5, 
  point.shape=19, 
  point.alpha=0.6, 
  text.add=TRUE,
  text.size=0.5,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  down.sample=10000,
  legend.add= TRUE
)
dev.off()


total = total[which(total@cluster != "15"),]
table(total@cluster)

meta<-matrix(data=NA,nrow=12578,ncol=2);
rownames(meta) = total@barcode
colnames(meta)=c("cluster","celltype");
meta[,1] = total@cluster
meta <- as.data.frame(meta)
cl <- meta$cluster
cl <- as.numeric(cl)
for ( i in 1:length(cl))
{
  if (cl[i] == 9 | cl[i] == 19  )
    meta[i,2] = "Immune_Cell"
  else if (cl[i] == 12)
    meta[i,2]  = "Cardiomyocyte"
  else if (cl[i] == 11 | cl[i] == 18)
    meta[i,2]  = "Smooth_Muscle_Cell"
  else if (cl[i] == 1 | cl[i] == 3 | cl[i] == 5 | cl[i] == 6 | cl[i] == 17 )
    meta[i,2]  = "Fibroblast"
  else
    meta[i,2]  = "Endothelial_Cell"
}
head(meta)
total@metaData$CellType = meta$celltype
table(total@metaData$CellType)

df = data.frame(x.sp@umap,x.sp@metaData$CellType)
df = data.frame(x.sp@umap,x.sp@metaData$cluster)
Celltype = x.sp@metaData$CellType
Cluster = x.sp@metaData$cluster

##cell type umap
ggplot(df, aes(x = umap.1, y = umap.2, color = Celltype)) +
  geom_point(size=0.8) + 
  theme_bw() + 
  scale_color_manual(values = c("#BC3C29FF","#0072B5FF", "#E18727FF","#20854EFF","#7876B1FF")) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("Celltype")

###Cluster umap
library(ggsci)
ggplot(df, aes(x = umap.1, y = umap.2, color = Cluster)) +
  geom_point(alpha=0.7) + 
  theme_bw() + 
  scale_color_manual(values = c("#E64B35B2","#00468BB2","#42B540B2" )) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("Cluster")

#
pdf("allcetype-celltype-in-umap-splitIC.pdf")
plotViz(
  obj=total,
  method="umap", 
  main="Cell type",
  point.color=total@metaData$CellType, 
  point.size=0.5, 
  point.shape=19, 
  point.alpha=0.6, 
  text.add=TRUE,
  text.size=0.8,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  down.sample=10000,
  legend.add= TRUE
)
dev.off()

#extract barcode
CM1_BC <- x.sp[which(x.sp@metaData$subtype == "CM1")]@barcode
CM2_BC <- x.sp[which(x.sp@metaData$subtype == "CM2")]@barcode
CM3_BC <- x.sp[which(x.sp@metaData$subtype == "CM3")]@barcode
CM4_BC <- x.sp[which(x.sp@metaData$subtype == "CM4")]@barcode
CM5_BC <- x.sp[which(x.sp@metaData$subtype == "CM5")]@barcode

setwd("/md01/Chenzh275/ori/Data/scATAC/cardiacmyocyte_regeneration/IDR/CellSubType/cardio/FindDAR/subtype-onlyD3-qval001/subtype-BC")
write.table(CM1_BC,"CM1-onlyD3.txt",row.names=F,quote =F)
write.table(CM2_BC,"CM2-onlyD3.txt",row.names=F,quote =F)
write.table(CM3_BC,"CM3-onlyD3.txt",row.names=F,quote =F)
write.table(CM4_BC,"CM4-onlyD3.txt",row.names=F,quote =F)
write.table(CM5_BC,"CM5-onlyD3.txt",row.names=F,quote =F)

#bioc process heatmap
setwd("/public/home/Chenzh275/Data/scATAC/cardiacmyocyte_regeneration/IDR/CellSubType/cardio/FindDAR/subtype-onlyD3-qval001/countreads-shuf/count/")
####事实证明用Gmat做热图不太理想，无论是promoter区还是genebody都比较乱，最后用将基因的TSS区提取出来进行count reads来做热图
#genes = read.table("/md01/Chenzh275/ori/Data/scATAC/cardiacmyocyte_regeneration/IDR/CellSubType/cardio/gene_site/gencode.vM16.gene.bed");
genes = read.table("~/ori/Data/scATAC/cardiacmyocyte_regeneration/IDR/CellSubType/cardio/gene_site/mm10.refseq.tss500bp.bed");
genes.gr = GRanges(genes[,1], IRanges(genes[,2], genes[,3]), name=genes[,4]  );

marker=read.table("../../../subtype_qval001/seperat_cluster9_subtype/peak/BiocProcess_Gene/cell_cycle_gene.txt",header=TRUE);
marker=read.table("../../../subtype_qval001/seperat_cluster9_subtype/peak/BiocProcess_Gene/glycolysis_gene.txt",header=TRUE);
marker=read.table("../../../subtype_qval001/seperat_cluster9_subtype/peak/BiocProcess_Gene/Antioxidant.txt",header=TRUE);
marker=read.table("../../../subtype_qval001/seperat_cluster9_subtype/peak/BiocProcess_Gene/ATPsysthesis.txt",header=TRUE);
marker=read.table("../../../subtype_qval001/seperat_cluster9_subtype/peak/BiocProcess_Gene/DNAdamage.txt",header=TRUE);
marker=read.table("../../../subtype_qval001/seperat_cluster9_subtype/peak/BiocProcess_Gene/Mitochondrial_electron_transport.txt",header=TRUE);
marker=read.table("../../../subtype_qval001/seperat_cluster9_subtype/peak/BiocProcess_Gene/PI3K-Akt_Singnal_pathway.txt",header=TRUE);
marker=read.table("../../../subtype_qval001/seperat_cluster9_subtype/peak/BiocProcess_Gene/Extracellular_Matrix.txt",header=TRUE);

marker.genes<-as.character(marker$Name)

marker.genes
genes.sel.gr <- genes.gr[which(genes.gr$name %in% marker.genes)];
genes.sel.gr

gene<-as.character(genes.sel.gr$name);
gene
x.sp = createGmatFromMat(
  obj=x.sp,
  input.mat="bmat",
  genes=genes.sel.gr,
  do.par=TRUE,
  num.cores=9
);

x.sp = scaleCountMatrix(
  obj=x.sp,
  cov=x.sp@metaData$passed_filters + 1,
  mat="gmat",
  method = "RPM"
);

x.sp = runMagic(
  obj=x.sp,
  input.mat="gmat",
  step.size=3
);
x.sp

amat<-matrix(data=NA,nrow=1590,ncol=80);
rownames(amat)=x.sp@metaData$subtype;
colnames(amat)=gene;
dim(amat)
for(j in 1:80){
  sam=x.sp@gmat[,gene[j]];
  amat[,j]=sam;           	
};
head(amat);

#promoter
write.csv(amat,"./BiocAccessibility/cardiomyocyte-CellCycle-promoter500bp[onlyD3].csv")
write.csv(amat,"./BiocAccessibility/cardiomyocyte-glycolysis-promoter500bp[onlyD3].csv")
write.csv(amat,"./BiocAccessibility/cardiomyocyte-Antioxidant-promoter500bp[onlyD3].csv")
write.csv(amat,"./BiocAccessibility/cardiomyocyte-ATPsysthesis-promoter500bp[onlyD3].csv")
write.csv(amat,"./BiocAccessibility/cardiomyocyte-DNAdamage-promoter500bp[onlyD3].csv")
write.csv(amat,"./BiocAccessibility/cardiomyocyte-Mitochondrial_electron_transport-promoter500bp[onlyD3].csv")
write.csv(amat,"./BiocAccessibility/cardiomyocyte-PI3K-Akt_Singnal_pathway-promoter500bp[onlyD3].csv")
write.csv(amat,"./BiocAccessibility/cardiomyocyte-Extracellular_Matrix-promoter500bp[onlyD3].csv")

###combine process and count mean
gly <-read.csv("./BiocAccessibility/cardiomyocyte-glycolysis-promoter500bp[onlyD3].csv")
anti <- read.csv("./BiocAccessibility/cardiomyocyte-Antioxidant-promoter500bp[onlyD3].csv")
ATP <- read.csv("./BiocAccessibility/cardiomyocyte-ATPsysthesis-promoter500bp[onlyD3].csv")
ddr <- read.csv("./BiocAccessibility/cardiomyocyte-DNAdamage-promoter500bp[onlyD3].csv")
mit <- read.csv("./BiocAccessibility/cardiomyocyte-Mitochondrial_electron_transport-promoter500bp[onlyD3].csv")
ECM <- read.csv("./BiocAccessibility/cardiomyocyte-Extracellular_Matrix-promoter500bp[onlyD3].csv")
PI3K <- read.csv("./BiocAccessibility/cardiomyocyte-PI3K-Akt_Singnal_pathway-promoter500bp[onlyD3].csv")
CC <- read.csv("./BiocAccessibility/cardiomyocyte-CellCycle-promoter500bp[onlyD3].csv")

head(CC)
dim(CC)
head(anti)

###将CM1-5按顺序排列
gly.cm1<-subset(gly,gly$X =="CM1")
gly.cm2<-subset(gly,gly$X =="CM2")
gly.cm3<-subset(gly,gly$X =="CM3")
gly.cm4<-subset(gly,gly$X =="CM4")
gly.cm5<-subset(gly,gly$X =="CM5")
gly.mean <-rbind(gly.cm1,gly.cm2)
gly.mean <-rbind(gly.mean,gly.cm3)
gly.mean <-rbind(gly.mean,gly.cm4)
gly.mean <-rbind(gly.mean,gly.cm5)
head(gly.mean)
gly.mean <- t(gly.mean)
write.csv(gly.mean,"./BiocAccessibility/tcardiomyocyte-glycolysis-promoter500bp[onlyD3].csv")

anti.cm1<-subset(anti,anti$X =="CM1")
anti.cm2<-subset(anti,anti$X =="CM2")
anti.cm3<-subset(anti,anti$X =="CM3")
anti.cm4<-subset(anti,anti$X =="CM4")
anti.cm5<-subset(anti,anti$X =="CM5")
anti.mean <-rbind(anti.cm1,anti.cm2)
anti.mean <-rbind(anti.mean,anti.cm3)
anti.mean <-rbind(anti.mean,anti.cm4)
anti.mean <-rbind(anti.mean,anti.cm5)
head(anti.mean)
anti.mean <- t(anti.mean)
write.csv(anti.mean,"./BiocAccessibility/tcardiomyocyte-Antioxidant-promoter500bp[onlyD3].csv")

ATP.cm1<-subset(ATP,ATP$X =="CM1")
ATP.cm2<-subset(ATP,ATP$X =="CM2")
ATP.cm3<-subset(ATP,ATP$X =="CM3")
ATP.cm4<-subset(ATP,ATP$X =="CM4")
ATP.cm5<-subset(ATP,ATP$X =="CM5")
ATP.mean <-rbind(ATP.cm1,ATP.cm2)
ATP.mean <-rbind(ATP.mean,ATP.cm3)
ATP.mean <-rbind(ATP.mean,ATP.cm4)
ATP.mean <-rbind(ATP.mean,ATP.cm5)
head(ATP.mean)
ATP.mean <- t(ATP.mean)
write.csv(ATP.mean,"./BiocAccessibility/tcardiomyocyte-ATPsysthesis-promoter500bp[onlyD3].csv")

ddr.cm1<-subset(ddr,ddr$X =="CM1")
ddr.cm2<-subset(ddr,ddr$X =="CM2")
ddr.cm3<-subset(ddr,ddr$X =="CM3")
ddr.cm4<-subset(ddr,ddr$X =="CM4")
ddr.cm5<-subset(ddr,ddr$X =="CM5")
ddr.mean <-rbind(ddr.cm1,ddr.cm2)
ddr.mean <-rbind(ddr.mean,ddr.cm3)
ddr.mean <-rbind(ddr.mean,ddr.cm4)
ddr.mean <-rbind(ddr.mean,ddr.cm5)
head(ddr.mean)
ddr.mean <- t(ddr.mean)
write.csv(ddr.mean,"./BiocAccessibility/tcardiomyocyte-DNAdamage-promoter500bp[onlyD3].csv")

mit.cm1<-subset(mit,mit$X =="CM1")
mit.cm2<-subset(mit,mit$X =="CM2")
mit.cm3<-subset(mit,mit$X =="CM3")
mit.cm4<-subset(mit,mit$X =="CM4")
mit.cm5<-subset(mit,mit$X =="CM5")
mit.mean <-rbind(mit.cm1,mit.cm2)
mit.mean <-rbind(mit.mean,mit.cm3)
mit.mean <-rbind(mit.mean,mit.cm4)
mit.mean <-rbind(mit.mean,mit.cm5)
head(mit.mean)
mit.mean <- t(mit.mean)
write.csv(mit.mean,"./BiocAccessibility/tcardiomyocyte-Mitochondrial_electron_transport-promoter500bp[onlyD3].csv")

ECM.cm1<-subset(ECM,ECM$X =="CM1")
ECM.cm2<-subset(ECM,ECM$X =="CM2")
ECM.cm3<-subset(ECM,ECM$X =="CM3")
ECM.cm4<-subset(ECM,ECM$X =="CM4")
ECM.cm5<-subset(ECM,ECM$X =="CM5")
ECM.mean <-rbind(ECM.cm1,ECM.cm2)
ECM.mean <-rbind(ECM.mean,ECM.cm3)
ECM.mean <-rbind(ECM.mean,ECM.cm4)
ECM.mean <-rbind(ECM.mean,ECM.cm5)
head(ECM.mean)
ECM.mean <- t(ECM.mean)
write.csv(ECM.mean,"./BiocAccessibility/tcardiomyocyte-ExtraCelluarMatrix-promoter500bp[onlyD3].csv")

PI3K.cm1<-subset(PI3K,PI3K$X =="CM1")
PI3K.cm2<-subset(PI3K,PI3K$X =="CM2")
PI3K.cm3<-subset(PI3K,PI3K$X =="CM3")
PI3K.cm4<-subset(PI3K,PI3K$X =="CM4")
PI3K.cm5<-subset(PI3K,PI3K$X =="CM5")
PI3K.mean <-rbind(PI3K.cm1,PI3K.cm2)
PI3K.mean <-rbind(PI3K.mean,PI3K.cm3)
PI3K.mean <-rbind(PI3K.mean,PI3K.cm4)
PI3K.mean <-rbind(PI3K.mean,PI3K.cm5)
head(PI3K.mean)
PI3K.mean <- t(PI3K.mean)
write.csv(PI3K.mean,"./BiocAccessibility/tcardiomyocyte-PI3K-Akt-Singnal-pathway-promoter500bp[onlyD3].csv")

CC.cm1<-subset(CC,CC$X =="CM1")
CC.cm2<-subset(CC,CC$X =="CM2")
CC.cm3<-subset(CC,CC$X =="CM3")
CC.cm4<-subset(CC,CC$X =="CM4")
CC.cm5<-subset(CC,CC$X =="CM5")
CC.mean <-rbind(CC.cm1,CC.cm2)
CC.mean <-rbind(CC.mean,CC.cm3)
CC.mean <-rbind(CC.mean,CC.cm4)
CC.mean <-rbind(CC.mean,CC.cm5)
head(CC.mean)
CC.mean <- t(CC.mean)
write.csv(CC.mean,"./BiocAccessibility/tcardiomyocyte-CellCycle-promoter500bp[onlyD3].csv")

#####combine all subtype and genes
###excel 处理过这些csv文件，把第一行去掉,并且添加一列term信息 在最后
gly <-read.csv("./BiocAccessibility/tcardiomyocyte-glycolysis-promoter500bp[onlyD3].csv")
anti <- read.csv("./BiocAccessibility/tcardiomyocyte-Antioxidant-promoter500bp[onlyD3].csv")
ATP <- read.csv("./BiocAccessibility/tcardiomyocyte-ATPsysthesis-promoter500bp[onlyD3].csv")
ddr <- read.csv("./BiocAccessibility/tcardiomyocyte-DNAdamage-promoter500bp[onlyD3].csv")
mit <- read.csv("./BiocAccessibility/tcardiomyocyte-Mitochondrial_electron_transport-promoter500bp[onlyD3].csv")
ECM <- read.csv("./BiocAccessibility/tcardiomyocyte-ExtraCelluarMatrix-promoter500bp[onlyD3].csv")
PI3K <- read.csv("./BiocAccessibility/tcardiomyocyte-PI3K-Akt-Singnal-pathway-promoter500bp[onlyD3].csv")
CC <- read.csv("./BiocAccessibility/tcardiomyocyte-CellCycle-promoter500bp[onlyD3].csv")

head(CC[,1:4])

all_sub_gene <- rbind(anti,gly)
all_sub_gene <- rbind(all_sub_gene ,ddr)
all_sub_gene <- rbind(all_sub_gene ,ECM)
all_sub_gene <- rbind(all_sub_gene ,PI3K)
all_sub_gene <- rbind(all_sub_gene ,CC)
all_sub_gene <- rbind(all_sub_gene ,mit)
all_sub_gene <- rbind(all_sub_gene ,ATP)

dim(all_sub_gene)
head(all_sub_gene[,1:10])

all_sub_gene =all_sub_gene[!duplicated(all_sub_gene$X),]
rownames(all_sub_gene) <- all_sub_gene$X
all_sub_gene <- all_sub_gene[,-1]
head(all_sub_gene[,1:10])
dim(all_sub_gene)

head(all_sub_gene[,1589:1590])

CM1<-all_sub_gene[,1:661]
CM2<-all_sub_gene[,662:838]
CM3<-all_sub_gene[,839:1023]
CM4<-all_sub_gene[,1024:1480]
CM5<-all_sub_gene[,1481:1590]

######col means matrix###########
CM1 <- data.frame(CM1= rowMeans(CM1))
CM2 <- data.frame(CM2= rowMeans(CM2))
CM3 <- data.frame(CM3= rowMeans(CM3))
CM4 <- data.frame(CM4= rowMeans(CM4))
CM5 <- data.frame(CM5= rowMeans(CM5))
dim(CM2)
all_sub_gene_mean <- cbind(CM1,CM2)
all_sub_gene_mean <- cbind(all_sub_gene_mean ,CM3)
all_sub_gene_mean <- cbind(all_sub_gene_mean ,CM4)
all_sub_gene_mean <- cbind(all_sub_gene_mean ,CM5)
all_sub_gene_mean$term <- all_sub_gene$term
head(all_sub_gene_mean)
write.csv(all_sub_gene_mean,"./BiocAccessibility/all_sub_gene_mean.newterm.csv")

##############colmeans########
library(preprocessCore)

all_sub_gene <-read.csv("./BiocAccessibility/all_sub_gene_mean.newterm.csv")
rownames(all_sub_gene) <- all_sub_gene$X
all_sub_gene <- all_sub_gene[,-1]
head(all_sub_gene)
all_sub_gene_sample <- all_sub_gene

all_sub_gene[,1:5] =t(scale(t(all_sub_gene[,1:5]),center = T))
all_sub_gene  = as.matrix(all_sub_gene[,1:5])
all_sub_gene[,1:5] <- normalize.quantiles(all_sub_gene[,1:5])
head(all_sub_gene)
all_sub_gene <- as.data.frame(all_sub_gene)
all_sub_gene$term = all_sub_gene_sample$term
dim(all_sub_gene)
table(all_sub_gene$term)

all_sub_gene <- all_sub_gene[which(all_sub_gene$term !="DDR" & all_sub_gene$term !="ECM" ),]
dim(all_sub_gene)
table(all_sub_gene$term)

library(ComplexHeatmap)
set.seed(123)
ha_column = HeatmapAnnotation(df = data.frame(subtype = c(rep("CM1", 1), rep("CM2", 1),rep("CM3", 1),rep("CM4", 1),rep("CM5", 1))),
                              col = list(subtype = c("CM1" =  "#BC3C29FF", "CM2" = "#0072B5FF", "CM3" =  "#E18727FF", "CM4" = "#20854EFF",
                                                     "CM5" =  "#7876B1FF")))
ha_row = rowAnnotation(df = data.frame(BP = c(rep("Antioxidant",12),rep("Gly",21),rep("Mit",11),
                                              rep("ATP",12),rep("PI3K",11),rep("CCR",77))), 
                       col = list(BP = c("ATP" = "#EE4C97FF", "Antioxidant" = "#FFDC91FF", "Gly" = "#7E6148B2",
                                         "PI3K" = "#6F99ADFF", "CCR" = "#60EFDBFF",
                                         "Mit" = "#00008B")))
pdf("./BiocAccessibility/heatmap-Bioc.pdf")
ht <-Heatmap(all_sub_gene[,1:5],clustering_method_rows = "ward.D2",cluster_rows = F,
             , show_row_names = T,  show_column_names = F,cluster_columns = F, show_heatmap_legend = T 
             , col =c("blue","white","red") , column_title="", top_annotation = ha_column)
draw(ht + ha_row, heatmap_legend_side = "right")
dev.off()

#用readcount 来做，取上面的基因
genes = read.table("~/ori/Data/scATAC/cardiacmyocyte_regeneration/IDR/CellSubType/cardio/gene_site/mm10.refseq.tss500bp.bed");
genes.gr = GRanges(genes[,1], IRanges(genes[,2], genes[,3]), name=genes[,4]  );

genes_sel <- genes.gr[which(genes.gr$name %in% rownames(all_sub_gene)),]
genes_sel <- as.data.frame(genes_sel)
write.table(genes_sel,"genes_select_Bioc.txt",quote=F,sep="\t")

#count reads后 画图
table <- read.table("./BiocAccessibility/genes_select_Bioc-ReadsCount.txt",header=T)
head(table)
gly_1 <- table[which(table$gene %in% gly$X),]
head(gly_1)
anti_1 <- table[which(table$gene %in% anti$X),]
ATP_1 <- table[which(table$gene %in% ATP$X),]
ddr_1 <- table[which(table$gene %in% ddr$X),]
mit_1 <- table[which(table$gene %in% mit$X),]
ECM_1 <- table[which(table$gene %in% ECM$X),]
PI3K_1 <-table[which(table$gene %in% PI3K$X),]
CC_1 <- table[which(table$gene %in% CC$X),]

head(mit_1)
all_sub_gene_1 <- rbind(anti_1,gly_1)
all_sub_gene_1 <- rbind(all_sub_gene_1 ,PI3K_1)
all_sub_gene_1 <- rbind(all_sub_gene_1 ,CC_1)
all_sub_gene_1 <- rbind(all_sub_gene_1 ,mit_1)
all_sub_gene_1 <- rbind(all_sub_gene_1 ,ATP_1)

all_sub_gene_1 =all_sub_gene_1[!duplicated(all_sub_gene_1$gene),] #去除相同的基因，因为有iosforms
head(all_sub_gene_1)

all_sub_gene_1[,5:9] =t(scale(t(all_sub_gene_1[,5:9]),center = T))
head(all_sub_gene_1)
all_sub_gene_1 <- as.data.frame(all_sub_gene_1)

all_sub_gene_sample <- all_sub_gene_sample[which(all_sub_gene_sample$term !="DDR" & all_sub_gene_sample$term !="ECM" ),]
all_sub_gene_1$term = all_sub_gene_sample$term

dim(all_sub_gene_1)
table(all_sub_gene_1$term)
head(all_sub_gene_1[5:9])
rownames(all_sub_gene_1) =all_sub_gene_1$gene

names <- data.frame(ori =rownames(all_sub_gene),new =all_sub_gene_1$gene)

write.table(all_sub_gene_1,"./BiocAccessibility/genes_select_Bioc-ReadsCount_term.txt",quote=F,sep="\t")

all_sub_gene_2 <- read.table("./BiocAccessibility/genes_select_Bioc-ReadsCount_term2.txt",header=T)
rownames(all_sub_gene_2) =all_sub_gene_2$gene
head(all_sub_gene_2)

library(pheatmap)
library("RColorBrewer")
set.seed(123)
annotation_row = data.frame(BP = c(rep("Antioxidant",12),rep("Gly",21),rep("PI3K",11),rep("CCR",39),
	rep("Mit",11), rep("ATP",12)))
ann_colors =list(subtype = c(CM1 =  "#BC3C29FF", CM2 = "#0072B5FF", 
	CM3 =  "#E18727FF", CM4 = "#20854EFF", CM5 =  "#7876B1FF"),
	BP = c(ATP = "#EE4C97FF", Antioxidant = "#FFDC91FF", 
		Gly= "#7E6148B2",PI3K = "#6F99ADFF", CCR = "#60EFDBFF",Mit = "#00008B"))
pdf("./BiocAccessibility/heatmap-Bioc-term.pdf",height =12)
ped =pheatmap(all_sub_gene_2[,5:9], cluster_rows=FALSE, 
              clustering_method = "ward.D2" ,
              annotation_row=annotation_row, annotation_colors = ann_colors,angle_col = "45",
              show_rownames=TRUE,treeheight_row = 0,
              color = colorRampPalette((rev(brewer.pal(n = 11, name="RdBu"))))(1000), 
              cluster_cols=F,main = "Cardiomyocyte Subtype Chromatin Accessibility",
              show_colnames = T,fontsize_row = 7.5)
dev.off()


#marker gene in CM
markergene <- c(
"Cacnb2","Camk2d","Cacna1e","Cacnb2","Camk2d","Irx5","Thra",
"Col3a1","Dpt","Postn","Col27a1","Col4a2",
"Sulf1","Adamts9","Tcf4","Smad1",
"S1pr1","Nppb","Pdgfb","Zfp143","Phip",
"Cdh7", "Sdk1"
)

genes = read.table("/md01/Chenzh275/ori/Data/scATAC/cardiacmyocyte_regeneration/IDR/CellSubType/cardio/gene_site/gencode.vM16.gene.bed");
genes.gr = GRanges(genes[,1], IRanges(genes[,2], genes[,3]), name=genes[,4]  );

genes.sel.gr <- genes.gr[which(genes.gr$name %in% markergene)];
genes.sel.gr

markergene <- markergene[which(markergene %in% genes.sel.gr$name)]
x.sp = createGmatFromMat(
    obj=x.sp, 
    input.mat="bmat",
    genes=genes.sel.gr,
    do.par=TRUE,
    num.cores=5
  )

x.sp = scaleCountMatrix(
    obj=x.sp, 
    cov=x.sp@metaData$passed_filters + 1,
    mat="gmat",
    method = "RPM"
  )

x.sp = runMagic(
    obj=x.sp,
    input.mat="gmat",
    step.size=3
  )

x.sp

pdf("Cardiomyocyte-CM1-CM5-specific-gene-in-umap.pdf")
for(i in 1:length(markergene)){
    plotFeatureSingle(
        obj=x.sp,
        feature.value=x.sp@gmat[,markergene[i]],
        method="umap", 
        main=markergene[i],
        point.size=0.8, 
        point.shape=19, 
        down.sample=10000,
#        quantiles=c(0, 1)
  )
}
dev.off()

pdf("Cardiomyocyte-CM1-CM5-specific-gene-in-vlnplot.pdf")
for(i in 1:length(markergene)){
dat = data.frame(x=x.sp@metaData[,"subtype"], y=x.sp@gmat[,markergene[i]]);
p1 <- ggplot(dat, aes(x=x, y=y, fill=x)) + 
  theme_classic() +
  geom_violin() + 
  xlab("subtype") +
  ylab("Accessibility") + 
  ggtitle(markergene[i]) +
  scale_fill_manual(values = c("#BC3C29FF","#0072B5FF", "#E18727FF","#20854EFF","#7876B1FF")) +
  theme(
      plot.margin = margin(5,1,5,1, "cm"),
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.ticks.x=element_blank(),
      legend.position = "none"
   );
print(p1)
}
dev.off()

#IC not work cos the lower cell number
IC<- total[which(total@cluster == 9 |total@cluster == 19 ),]

IC = runDiffusionMaps(
    obj=IC,
    input.mat="bmat", 
    num.eigs=50
  )

plotDimReductPW(
    obj=IC, 
    eigs.dims=1:50,
    point.size=0.3,
    point.color="grey",
    point.shape=19,
    point.alpha=0.6,
    down.sample=,
    pdf.file.name=NULL, 
    pdf.height=7, 
    pdf.width=7
  );

library(leiden);
library(umap)
IC = runKNN(
  obj=IC,
  eigs.dims=1:10,
 k=10
)

IC=runCluster(
  obj=IC,
  tmp.folder=tempdir(),
  louvain.lib="leiden",
  seed.use=18
)

IC = runViz(
  obj=IC,
  tmp.folder=tempdir(),
  dims=2,
  eigs.dims=1:10,
  method="umap",
  seed.use=18
)

markergene <- c("Cd74","Cd22","Cd79b","Cxcr4", #B cell
	"Samsn1","Clec4d","S100a9", "S100a8",#Granulocy
	"Cd86","Mpeg1","Cd68","Ccl3","Cd44",# MP
	"Cd3e","Cd3d", "Cd8a", "Cd8b1","nkg7","Igfbp4","Lat","Cd53","Itk" #T cell
	)

genes = read.table("/md01/Chenzh275/ori/Data/scATAC/cardiacmyocyte_regeneration/IDR/CellSubType/cardio/gene_site/gencode.vM16.gene.bed");
genes.gr = GRanges(genes[,1], IRanges(genes[,2], genes[,3]), name=genes[,4]  );

genes.sel.gr <- genes.gr[which(genes.gr$name %in% markergene)];
genes.sel.gr

markergene <- markergene[which(markergene %in% genes.sel.gr$name)]

IC = createGmatFromMat(
    obj=IC, 
    input.mat="bmat",
    genes=genes.sel.gr,
    do.par=TRUE,
    num.cores=5
  )

IC = scaleCountMatrix(
    obj=IC, 
    cov=IC@metaData$passed_filters + 1,
    mat="gmat",
    method = "RPM"
  )

IC = runMagic(
    obj=IC,
    input.mat="gmat",
    step.size=3
  )

IC

pdf("IC-marker-gene-annotation-in-umap-splitIC.pdf")
for(i in 1:21){
    plotFeatureSingle(
        obj=IC,
        feature.value=IC@gmat[,markergene[i]],
        method="umap", 
        main=markergene[i],
        point.size=0.4, 
        point.shape=19, 
        down.sample=10000,
#        quantiles=c(0, 1)
  )
}
dev.off()

