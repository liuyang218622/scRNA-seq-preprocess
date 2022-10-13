rm(list=ls())
options(stringsAsFactors = F)
##
setwd("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/04.output/01.数据搜集和预处理/")
##载入TCGA_LIHC数据
TCGALIHC_fpkm<- read.table('../../01.data/13.TCGA-LIHC.htseq_fpkm.tsv',header = T)
TCGALIHC_counts<-read.table('../../01.data/12.TCGA-LIHC.htseq_counts.tsv',header = T)
TCGALIHC_counts[,2:ncol(TCGALIHC_counts)]=round((2**TCGALIHC_counts[,2:ncol(TCGALIHC_counts)])-1)
TCGALIHC_fpkm[,2:ncol(TCGALIHC_fpkm)]=round(2**TCGALIHC_fpkm[,2:ncol(TCGALIHC_fpkm)]-1,3)
TCGALIHC_counts$Ensembl_ID<-sapply(strsplit(TCGALIHC_counts$Ensembl_ID,'.',fixed=T),function(x)x[1])
TCGALIHC_fpkm$Ensembl_ID<-sapply(strsplit(TCGALIHC_fpkm$Ensembl_ID,'.',fixed=T),function(x)x[1])
##注释文件##
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75
columns(edb)
keytypes(edb)
grids <- keys(edb, keytype="GENEID")
## Get the data
gene2sym<-select (edb, keys=grids, 
                 columns=c("SYMBOL","GENEBIOTYPE",'GENENAME'),
                 keytype="GENEID")
save(gene2sym,file="/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/04.output/01.数据搜集和预处理/gene2sym.Rdata")
load("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/04.output/01.数据搜集和预处理/gene2sym.Rdata")
table(gene2sym$GENEBIOTYPE)
##mRNA的注释
mRNA<-gene2sym[which(gene2sym$GENEBIOTYPE=="protein_coding"),]
###子函数，计算多个序列号比对到同一个基因的矩阵的表达量（取平均数最大的那个）
merge_genes_matrix<-function(exprSet,ids){
  exprSet=exprSet[rownames(exprSet) %in% ids$GENEID,]
  ids=ids[match(rownames(exprSet),ids$GENEID),]
  tmp=by(exprSet,ids$SYMBOL,
         function(x) rownames(x)[which.max(rowMeans(x))] )
  probes = as.character(tmp)
  exprSet=exprSet[rownames(exprSet) %in% probes ,]
  rownames(exprSet)=ids[match(rownames(exprSet),ids$GENEID),2]
  return (exprSet)
}
##提取mRNA的counts数据
TCGALIHC_mRNA_counts=na.omit(TCGALIHC_counts[match(mRNA$GENEID,TCGALIHC_counts$Ensembl_ID),])
rownames(TCGALIHC_mRNA_counts)=TCGALIHC_mRNA_counts$Ensembl_ID
TCGALIHC_mRNA_counts=TCGALIHC_mRNA_counts[,-1]
#同一个基因不同ensmbl号，取平均值最大的那个基因的情况。
TCGALIHC_mRNA_gene_counts<-merge_genes_matrix(TCGALIHC_mRNA_counts,mRNA)
##只取质量为A的样本
TCGALIHC_mRNA_gene_counts=TCGALIHC_mRNA_gene_counts[,substr(colnames(TCGALIHC_mRNA_gene_counts),16,16)%in%c("A")]
write.table(TCGALIHC_mRNA_gene_counts,file = '../01.数据搜集和预处理/02.TCGALIHC_mRNA_gene_counts.txt',sep='\t',quote = F)
##转换mRNA fpkm的ID
rownames(TCGALIHC_fpkm)=TCGALIHC_fpkm$Ensembl_ID
TCGALIHC_fpkm=TCGALIHC_fpkm[,-1]
TCGALIHC_gene_fpkm<-merge_genes_matrix(TCGALIHC_fpkm,gene2sym)
TCGALIHC_gene_fpkm=TCGALIHC_gene_fpkm[,substr(colnames(TCGALIHC_gene_fpkm),16,16)%in%c("A")]
write.table(TCGALIHC_gene_fpkm,file = '/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/04.output/01.数据搜集和预处理/03.TCGALIHC_gene_fpkm.txt',sep = '\t',quote = F)
###
TCGALIHC_gene_fpkm<-read.table('/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/04.output/01.数据搜集和预处理/03.TCGALIHC_gene_fpkm.txt',header = T)
TCGALIHC_gene_tumor_fpkm=TCGALIHC_gene_fpkm[ ,substr(colnames(TCGALIHC_gene_fpkm),14,14) %in% c("0")]
library(limma)
rt2 = TCGALIHC_gene_tumor_fpkm
rt2 = as.matrix(rt2)
mrna_datexpr_LIHC = avereps(rt2)
type<-unlist(lapply(colnames(mrna_datexpr_LIHC),function(b){substr(b,14,16)}))
index<-which(type %in% c("01A","11A"))
mrna_datexpr_LIHC<-mrna_datexpr_LIHC[,index] ###363   6  19 
colnames(mrna_datexpr_LIHC)<-unlist(lapply(colnames(mrna_datexpr_LIHC),function(b){substr(b,0,15)}))
#############临床信息
tcga_os<-read.table('/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/01.data/15.TCGA-LIHC.survival.tsv',header = T)
tcga_os<-tcga_os[,-3]
tcga_os$sample<-gsub('A$','',tcga_os$sample)
tcga_os$sample<-gsub('-','.',tcga_os$sample)
colnames(tcga_os)=c("sample","OS_STATUS","OS_TIME")
length(intersect(tcga_os[,1],colnames(mrna_datexpr_LIHC)))
cancersample<-intersect(tcga_os[,1],colnames(mrna_datexpr_LIHC))
tcga_os<-tcga_os[match(cancersample,tcga_os[,1]),] ####finally 419 samples
tcga_os$OS_STATUS<-plyr::revalue(as.factor(tcga_os$OS_STATUS),c("1"="DECEASED","0"="LIVING"))
tcga_os=subset(tcga_os,tcga_os$sample %in% colnames(mrna_datexpr_LIHC) )
mrna_datexpr_cancer<-mrna_datexpr_LIHC[,match(cancersample,colnames(mrna_datexpr_LIHC))] ###53317个基因 363个肿瘤样本
table(tcga_os$sample %in% colnames(mrna_datexpr_cancer))
## step5###下载GEO数据库的数据###
##################GEO1 
library("GEOquery")
library(gdata)
gse1 <- getGEO("GSE14520",GSEMatrix=TRUE)
geo_pdata_raw <-pData(gse1[[1]])
geo_exp_raw1 <- exprs(gse1[[1]])
library(dplyr)
probe <- getGEO(filename='/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/04.output/01.数据搜集和预处理/GEO/GSE14520_family.soft')
#############
id_probe1=probe@gpls$GPL571@dataTable@table
probe2gene1 <- id_probe1[,c(1,11)] 
library(stringr)  
probe2gene1$`Gene Symbol`=trimws(str_split(probe2gene1$`Gene Symbol`,'///',simplify = T)[,1])
###############
id_probe2=probe@gpls$GPL3921@dataTable@table
probe2gene2 <- id_probe2[,c(1,11)] 
library(stringr)  
probe2gene2$`Gene Symbol`=trimws(str_split(probe2gene2$`Gene Symbol`,'///',simplify = T)[,1])
######
probe2id<-rbind(probe2gene1,probe2gene2)
colnames(probe2id)=c("ID","symbol")
a<-intersect(rownames(geo_exp_raw1),probe2id[,1])
probe2id<-probe2id[match(a,probe2id[,1]),]
geo_exp_raw1<-geo_exp_raw1[match(a,rownames(geo_exp_raw1)),]
datExpr1_v1<-t(sapply(split(as.data.frame(geo_exp_raw1),probe2id[,2]),function(x) colMeans(x[,1:ncol(x)])))
datExpr1_v1<-datExpr1_v1[-1,] ### 13237个基因 445个样本(包含正常和肝癌组织)
geo_pdata_raw=geo_pdata_raw[geo_pdata_raw$`Tissue:ch1`%in% "Liver Tumor Tissue",] ## 222 个肝癌组织
datExpr1=datExpr1_v1[,rownames(geo_pdata_raw)]
####预后信息###
geo_pdata=read.table('/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/04.output/01.数据搜集和预处理/GEO/GSE14520_Extra_Supplement.txt',header = T,fill = T,sep = '\t')
geo_pdata1<- geo_pdata[,c(3,19,20)]
geo_pdata1_os<-geo_pdata1[,c(1,2,3)]
geo_pdata1_os$Survival.months= geo_pdata1_os$Survival.months*30
geo_pdata1_os=geo_pdata1_os[geo_pdata1$Affy_GSM %in% rownames(geo_pdata_raw),] ###222个样本
colnames(geo_pdata1_os)<-c("sample","OS_STATUS","OS_TIME")
rownames(geo_pdata1_os)=NULL
geo_pdata1_os$OS_STATUS<-plyr::revalue(as.factor(geo_pdata1_os$OS_STATUS),c("1"="DECEASED","0"="LIVING"))
datExpr1<-datExpr1[,which(colnames(datExpr1) %in% geo_pdata1_os[,1])] ###
datExpr1<-na.omit(datExpr1) ##### 13237个基因 222个肿瘤样本
datExpr1<-2**datExpr1 ##这个GEO数据库经历了log2 
table(colnames(datExpr1) %in% geo_pdata1_os$sample)

#############获得GEO2的数据################
gse2 <- getGEO("GSE10186",GSEMatrix=TRUE)
geo_pdata_raw <-pData(gse2[[1]])
geo_exp_raw1 <- exprs(gse2[[1]])
probe <- getGEO(filename='/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/04.output/01.数据搜集和预处理/GEO/GSE10186_family.soft')
id_probe=probe@gpls$GPL5474@dataTable@table
probe2id<-id_probe[,c(1,5)]
a<-intersect(rownames(geo_exp_raw1),probe2id[,1])
probe2id<-probe2id[match(a,probe2id[,1]),]
geo_exp_raw1<-geo_exp_raw1[match(a,rownames(geo_exp_raw1)),]
datExpr2<-t(sapply(split(as.data.frame(geo_exp_raw1),probe2id[,2]),function(x) colMeans(x[,1:ncol(x)])))
#geo_pdata_raw<-geo_pdata_raw[which(geo_pdata_raw[,8] %in% "Bladder tumor"),]
geo_pdata2<- geo_pdata_raw[,c(2,51,52)]
geo_pdata2_os<-geo_pdata2[,c(1,2,3)]
colnames(geo_pdata2_os)<-c("sample","OS_STATUS","OS_TIME")
rownames(geo_pdata2_os)=NULL
geo_pdata2_os$OS_STATUS<-plyr::revalue(as.factor(geo_pdata2_os$OS_STATUS),c("alive or censored, 1: dead) : 1"="DECEASED","alive or censored, 1: dead) : 0"="LIVING"))
datExpr2<-datExpr2[,which(colnames(datExpr2) %in% geo_pdata2_os[,1])]
datExpr2<-na.omit(datExpr2) 
table(colnames(datExpr2) %in% geo_pdata2_os$sample) ### 6100个基因 118个样本
####ICGC##数据##
############ICGC1#############
ICGC1<-read.table('/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/04.output/01.数据搜集和预处理/ICGC/LICA-FR/LICA-FR-express.txt',header = T)
ICGC1$gene<-sapply(strsplit(ICGC1$gene,'.',fixed=T),function(x)x[1])
rownames(ICGC1)=ICGC1$gene
ICGC1=ICGC1[,-1]
ICGC1_gene_fpkm<-merge_genes_matrix(ICGC1,gene2sym)
colnames(ICGC1_gene_fpkm)=gsub('T$','',colnames(ICGC1_gene_fpkm)) ##
###预后信息##
icgc1<-read.table('/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/04.output/01.数据搜集和预处理/ICGC/LICA-FR/donor.LICA-FR.tsv',header = T,sep = '\t',fill = T)
icgc1_os<-icgc1[,c(4,6,17)]
colnames(icgc1_os)<-c("sample","OS_STATUS","OS_TIME")
icgc1_os$OS_STATUS<-plyr::revalue(as.factor(icgc1_os$OS_STATUS),c("deceased"="DECEASED","alive"="LIVING"))
icgc1_os=icgc1_os[icgc1_os$sample %in% colnames(ICGC1_gene_fpkm),] # 153个样本
icgc1_datExpr1<-ICGC1_gene_fpkm[,which(colnames(ICGC1_gene_fpkm) %in% icgc1_os[,1])] ##55773 153样本
############ICGC2##############
ICGC2<-read.table('/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/04.output/01.数据搜集和预处理/ICGC/LINC-RIKEN-JP/LICA-LICA-FR-express.txtFR-express.txt',header = T)
rownames(ICGC2)=ICGC2$gene
ICGC2=ICGC2[,-1]
ICGC2=ICGC2[,substr(colnames(ICGC2),7,12) %in% c("Cancer")] ##22913 243
colnames(ICGC2)=gsub('_Cancer','',colnames(ICGC2))
#####预后信息###
icgc2<-read.table('../01.数据搜集和预处理/ICGC/LINC-RIKEN-JP/donor.LIRI-JP.tsv',header = T,sep="\t",fill =T)
icgc2_os<-icgc2[,c(4,6,17)]
colnames(icgc2_os)<-c("sample","OS_STATUS","OS_TIME")
icgc2_os$OS_STATUS<-plyr::revalue(as.factor(icgc2_os$OS_STATUS),c("deceased"="DECEASED","alive"="LIVING"))
icgc2_datExpr2<-ICGC2[,which(colnames(ICGC2) %in% icgc2_os[,1])] ##22913 231样本
icgc2_os= unique(icgc2_os[colnames(ICGC2) %in% icgc2_os[,1], ]) ## 221 个样本
#############################合并表达谱
tmp1<-intersect(rownames(icgc1_datExpr1),rownames(icgc2_datExpr2))
tmp2<-intersect(tmp1,intersect(rownames(datExpr1),rownames(datExpr2)))
genes<-intersect(rownames(mrna_datexpr_cancer),tmp2)
tcga<-mrna_datexpr_cancer[match(genes,rownames(mrna_datexpr_cancer)),]
geo1<-datExpr1[match(genes,rownames(datExpr1)),]
geo2<-datExpr2[match(genes,rownames(datExpr2)),]
icgc1<-icgc1_datExpr1[match(genes,rownames(icgc1_datExpr1)),]
icgc2<-icgc2_datExpr2[match(genes,rownames(icgc2_datExpr2)),]
datexpr<-cbind(tcga,geo1,geo2,icgc1,icgc2)
datexpr<-na.omit(datexpr) ## 4188个基因  1087个样本
tcga_os=unique(tcga_os)
geo_pdata1_os=unique(geo_pdata1_os)
geo_pdata2_os=unique(geo_pdata2_os)
icgc1_os=unique(icgc1_os)
icgc2_os=unique(icgc2_os)
clin<-rbind(tcga_os,geo_pdata1_os,geo_pdata2_os,icgc1_os,icgc2_os)
clin<-unique(clin) ## 1103 个样本
write.table(clin,file = '/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Table/SupplementaryTable1_clin.txt',row.names = F,sep = '\t')
table(colnames(datexpr) %in% clin[,1]) ##False 12 # True 1075 
#############################矫正表达谱
sampletype<-data.frame(sample=c(colnames(tcga),colnames(geo1),colnames(geo2),colnames(icgc1),colnames(icgc2)),
                       type=c(rep("LIHC",ncol(tcga)),rep("GSE1",ncol(geo1)),rep("GSE2",ncol(geo2)),rep("ICGC1",ncol(icgc1)),rep("ICGC2",ncol(icgc2)))
)
library("sva")
colnames(datexpr)<-gsub("[.]","-",colnames(datexpr))
exp<-as.matrix(datexpr)
dimnames<-list(rownames(exp),colnames(exp))
exp<-matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
batchtype<-sampletype[,2]
exp<-data.frame(ComBat(exp,batchtype))
#######用PCA检查数据情况####
if(FALSE){
  library(ggplot2)
  data<-t(datexpr)  ###data行是样本列是基因
  data.pca <- prcomp(data,retx=T,scale=T,center=T) 
  a <- summary(data.pca)
  tmp <- a$importance
  pro1 <- as.numeric(sprintf("%.3f",tmp[2,1]))*100
  pro2 <- as.numeric(sprintf("%.3f",tmp[2,2]))*100
  pc = as.data.frame(a$x)
  pc1<-pc$PC1
  pc2<-pc$PC2
  #type<-c(rep("TCGA",ncol(datexpr_tcga)),rep("CGGA_mRNAseq_325",ncol(datexpr_cgga_part1)),rep("CGGA_mRNAseq_693",ncol(datexpr_cgga_part2)))
  type<-sampletype[,2]
  gg1<-data.frame(x=pc1,y=pc2,Batch=type,t=rep("Before batch correction",length(type)))
  rownames(gg1)<-rownames(pc)
  
  data<-t(exp)  ###data行是样本列是基因
  data.pca <- prcomp(data,retx=T,scale=T,center=T) 
  a <- summary(data.pca)
  tmp <- a$importance
  pro1 <- as.numeric(sprintf("%.3f",tmp[2,1]))*100
  pro2 <- as.numeric(sprintf("%.3f",tmp[2,2]))*100
  pc = as.data.frame(a$x)
  pc1<-pc$PC1
  pc2<-pc$PC2
  type<-sampletype[,2]
  gg2<-data.frame(x=pc1,y=pc2,Batch=type,t=rep("After batch correction",length(type)))
  rownames(gg2)<-rownames(pc)
  
  gg<-as.data.frame(rbind(gg1,gg2))
  
  ###################plot figureA
  p<-ggplot(gg,aes(x=x,y=y,colour=Batch))+geom_point()
  a<-p+facet_grid(.~t)
  b<-a+theme (panel.background=element_blank(), panel.border=element_rect(colour="black", fill=NA))
  c<-b+labs(x = "PC 1", y = "PC 2")
  d<-c+theme(strip.background=element_rect(fill=NA),strip.text=element_text(size=rel(1.2)))
  e<-d+labs(colour="Batch")
  f<-e+theme_bw()+theme(text=element_text(size=15))
}
###########提取RBP基因####
RBP<-read.table("../01.数据搜集和预处理/01.RBP_genelist.txt",header = T)
RBP<-RBP[,1]
length(intersect(RBP,rownames(exp)))
RBP_genes<-intersect(RBP,rownames(exp)) ### 1374个基因
save(list=ls(),file="../step1.RData")
########################一致性聚类
library(ConsensusClusterPlus)
library(survival, quietly = TRUE)
library(survminer, quietly = TRUE)
#load("./runtime/output/data.RData")
kDefaults <- list(
  clinTsv = "../../../Reuslts/Table/SupplementaryTable1_clin.txt",
  ConClustDir = "../02.RBP基因无监督聚类/cc1/",
  RBPgroupfile = "/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Table/SupplementaryTable2_RBP_cluster.txt",
  survRBPclusterPdf = "../../../Reuslts/Figure/Figure1B_surv_metabolic_cluster.pdf",  # 根据RBP分类的生存曲线
  survRBPclusterTif72 = "../../../Reuslts/Figure/Figure1B_surv_metabolic_cluster-72.tif",  # 根据RBP分类的生存曲线
  survRBPclusterTif300 = "../../../Reuslts/Figure/Figure1B_surv_metabolic_cluster-300.tif"  # 根据RBP分类的生存曲线
)

getOptK <- function(conClust, minCls = 2, maxCls = 10) {
  #   最佳分类数
  Kvec = minCls: maxCls
  x1 = 0.1
  x2 = 0.9 # threshold defining the intermediate sub-interval
  
  PAC = rep(NA, length(Kvec))
  names(PAC) = paste("K=", Kvec, sep = "") # from 2 to maxK
  for(i in Kvec){
    M = conClust[[i]][["consensusMatrix"]]
    Fn = ecdf(M[lower.tri(M)])
    PAC[i-1] = Fn(x2) - Fn(x1)
  }
  optK = Kvec[which.min(PAC)]
  return(optK)
}
####step1: Unsupervised clustering for 21 m6A regulators（m6Acluster）
#d<-sweep(exp,1,apply(exp,1,median,na.rm=TRUE))
emt<-exp
conClust <- ConsensusClusterPlus(
  as.matrix(emt[which(rownames(emt) %in% RBP_genes),]),
  maxK = 5,
  reps = 500,
  pItem = 0.8,
  pFeature = 1,
  clusterAlg = "km",
  distance = 'euclidean',
  corUse = "complete.obs",
  seed = 123456,
  plot = "pdf",
  #title = kDefaults$ConClustDir,
  finalLinkage = 'ward.D',
  innerLinkage = 'ward.D',
  writeTable = FALSE
)
save(conClust,file="../../04.output/conClust.RData")
k<-getOptK(conClust, minCls = 2, maxCls = 5)
message(sprintf("最佳分类数：%d", k))

cluster1<-conClust[[k]]$consensusClass
cluster<-plyr::revalue(as.character(cluster1),c("1"="RBP_clusterA","2"="RBP_clusterB"))
RBP_cluster=data.frame(sample=names(cluster1),cluster=cluster)
RBP_cluster$sample=gsub('-','.',RBP_cluster$sample)
write.table(RBP_cluster,file=kDefaults$RBPgroupfile,quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)

###############
clin1<-read.table(file=kDefaults$clinTsv, header=TRUE,sep="\t",as.is=TRUE)
clin1<-clin1[which(clin1[,3]>0),]
samples<-intersect(RBP_cluster[,1],clin1[,1])
group<-RBP_cluster[match(samples,RBP_cluster[,1]),] #####筛选有临床信息的样本
clin1<-clin1[match(samples,clin1[,1]),]
table(group[,1]==clin1[,1])
survRBPclusterData <- data.frame(
  time = as.numeric(clin1$OS_TIME),
  status = as.numeric(clin1$OS_STATUS == "DECEASED"),
  group = group[,2]
)
# 绘制生存曲线"#DF342A","#DCAB3A","#559ECA"(红，黄，蓝) 黄色"#FC7F23" 蓝色"#3B7FB6"
survRBPcluster <- ggsurvplot(
  survfit(Surv(time, status) ~ group, data = survRBPclusterData),
  font.legend=14,xlab="OS Time (days)",legend.title="cluster",pval=TRUE,
  font.y=c(20),font.x=c(20),font.tickslab =c(18),pval.size=7,palette = c("#FF7F00","#377EB8","#DCAB3A"),risk.table = TRUE,
  ggtheme = theme_survminer()
)
#survm6Acluster
#write.table(group,file=kDefaults$m6Agroupfile,quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)
pdf(kDefaults$survRBPclusterPdf, width = 10, height = 8)
print(survRBPcluster)
dev.off()
tiff(kDefaults$survRBPclusterTif72, width = 10, height = 8, units = "in", res = 72)
print(survRBPcluster)
dev.off()
tiff(kDefaults$survRBPclusterTif300, width = 10, height = 8, units = "in", res = 300)
print(survRBPcluster)
dev.off()

###############GSVA分析
library("GSEABase")
library(GSVA)
library(pheatmap)
###创建文件#########
kDefaults <- list(
  gmt = "../../01.data/11.c2.cp.kegg.v7.1.symbols.gmt",
  gsvaFile_output = "../../../Reuslts/Table/SupplementaryTable3_GSVA.txt",
  wilcox_resultAB = "../../../Reuslts/Table/SupplementaryTable3_GSVA_wilcoxAB.txt",
  wilcox_resultBC = "../../../Reuslts/Table/SupplementaryTable3_GSVA_wilcoxBC.txt",
  RBPgroupfile = "../../../Reuslts/Table/SupplementaryTable2_RBP_cluster.txt",
  GSVA1Pdf = "../../../Reuslts/Figure/Figure2_GSVA1.pdf", 
  GSVA1Tif72 = "../../../Reuslts/Figure/Figure2_GSVA1-72.tif",  
  GSVA1Tif300 = "../../../Reuslts/Figure/Figure2_GSVA1-300.tif",
  GSVA2Pdf = "../../../Reuslts/Figure/Figure2_GSVA2.pdf", 
  GSVA2Tif72 = "../../../Reuslts/Figure/Figure2_GSVA2-72.tif",  
  GSVA2Tif300 = "../../../Reuslts/Figure/Figure2_GSVA2-300.tif",
  ssGSEA_file = "../../../Reuslts/Table/SupplementaryTable4_ssGSVA.txt",
  p1Pdf = "../../../Reuslts/Figure/Figure3A_ssGSVA.pdf", 
  p1Tif72 = "../../../Reuslts/Figure/Figure3A_ssGSVA-72.tif",  
  p1Tif300 = "../../../Reuslts/Figure/Figure3A_ssGSVA-300.tif",
  p2pdf="../../../Reuslts/Figure/Figure3B_ssGSVA.pdf",
  p2Tif72 = "../../../Reuslts/Figure/Figure3B_ssGSVA-72.tif", 
  p2Tif300 = "../../../Reuslts/Figure/Figure3B_ssGSVA-300.tif"
  
)
############计算GSVA
geneSets <- getGmt(kDefaults$gmt)
colnames(exp)=gsub('-','.',colnames(exp))
mydata = as.matrix(exp)
res_es_TCGA_cancer <- gsva(mydata, geneSets, mx.diff=FALSE, verbose=FALSE,method='gsva')
write.table(res_es_TCGA_cancer,file=kDefaults$gsvaFile_output,quote=FALSE,sep="\t",col.names=TRUE,row.names=TRUE)

data<-read.table(kDefaults$RBPgroupfile,header=TRUE,sep="\t",as.is=TRUE)
data$sample=gsub('-','.',data$sample)
sampletype<-data.frame(samples=c(colnames(tcga),colnames(geo1),colnames(geo2),colnames(icgc1),colnames(icgc2)),
                       project=c(rep("LIHC",ncol(tcga)),rep("GSE1",ncol(geo1)),rep("GSE2",ncol(geo2)),
                                 rep("ICGC1",ncol(icgc1)),rep("ICGC2",ncol(icgc2)))
)

gg<-merge(data,sampletype,by.x="sample",by.y="samples")
colnames(res_es_TCGA_cancer)=gsub('-','.',colnames(res_es_TCGA_cancer))
res_es_TCGA_cancer<-res_es_TCGA_cancer[,match(gg[,1],colnames(res_es_TCGA_cancer))]
wilcox_result_AB <- apply(res_es_TCGA_cancer, 1, function(geneVec){
  resultAB <- wilcox.test(geneVec[which(gg[,2] == "RBP_clusterA")], geneVec[which(gg[,2] == "RBP_clusterB")], paired = FALSE)
  return(c(resultAB$p.value))
})
#p值矫正
wilcox_result_adjust_AB <- p.adjust(wilcox_result_AB, method = "bonferroni")
pathway_AB<-names(wilcox_result_adjust_AB)[which(wilcox_result_adjust_AB<0.05)]
pathway_top20_AB<-names(sort(wilcox_result_adjust_AB))[1:15]
write.table(cbind(wilcox_result_AB=wilcox_result_AB,wilcox_result_adjust_AB=wilcox_result_adjust_AB),
            file=kDefaults$wilcox_resultAB,quote=FALSE,sep="\t",col.names=TRUE,row.names=TRUE)
###########比较A和B
gg_AB=subset(gg,gg$cluster%in%c("RBP_clusterA","RBP_clusterB"))
annotation_col = data.frame(gg_AB[,c(3,2)]) 
rownames(annotation_col)<-gg_AB[,1]
annotation_col = annotation_col[order(annotation_col$cluster),]
ann_colors = list(
  cluster =c("RBP_clusterA"="#FF7F00","RBP_clusterB"="#377EB8")
)
sample_order<-rownames(annotation_col)
#sample_order<-subset(sample_order,colnames(res_es_TCGA_cancer) %in% sample_order)
pdf(kDefaults$GSVA1Pdf, width = 12, height = 8)
pheatmap(res_es_TCGA_cancer[pathway_top20_AB,sample_order],clustering_method = "ward.D",show_rownames =T,show_colnames =F,fontsize_row=10,
         annotation_col = annotation_col,cluster_cols = F,annotation_colors = ann_colors,color = colorRampPalette(c("navy", "white", "firebrick3"))(100)) 
dev.off()
tiff(kDefaults$GSVA1Tif72, width = 12, height = 8, units = "in", res = 72)
pheatmap(res_es_TCGA_cancer[pathway_top20_AB,sample_order],clustering_method = "ward.D",show_rownames =T,show_colnames =F,fontsize_row=10,
         annotation_col = annotation_col,cluster_cols = F,annotation_colors = ann_colors,color = colorRampPalette(c("navy", "white", "firebrick3"))(100)) 
dev.off()
tiff(kDefaults$GSVA1Tif300, width = 12, height = 8, units = "in", res = 300)
pheatmap(res_es_TCGA_cancer[pathway_top20_AB,sample_order],clustering_method = "ward.D",show_rownames =T,show_colnames =F,fontsize_row=10,
         annotation_col = annotation_col,cluster_cols = F,annotation_colors = ann_colors,color = colorRampPalette(c("navy", "white", "firebrick3"))(100)) 
dev.off()
##########################
library(openxlsx)
group<-read.table(kDefaults$RBPgroupfile,header=TRUE,sep="\t",as.is=TRUE)
group$sample=gsub("-",".",group$sample)
gene_set<-read.csv("../01.数据搜集和预处理/04.mmc3.csv")
head(gene_set)
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
ssGSEA_result<- gsva(as.matrix(exp),list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
write.table(ssGSEA_result,file=kDefaults$ssGSEA_file,quote=FALSE,sep="\t",col.names=TRUE,row.names=TRUE)

w_df = reshape2::melt(t(ssGSEA_result), varnames = c("sample", "ssGSEA"))
gg<-merge(w_df,group,by.x="sample",by.y="sample")
p<-ggboxplot(gg,x='ssGSEA',y='value',fill='cluster')+
  stat_compare_means(aes(group = cluster),label = "p.signif")+xlab('')+ylab('Immune score')
p2<-p+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +scale_fill_manual(values=c("#FF7F00", "#377EB8","#DF342A"))
print(p2)
pdf(kDefaults$p1Pdf,width=15,height=8)
print (p2)
dev.off()
tiff(kDefaults$p1Tif72,width=15,height=8,units="in",res=72)
print (p2)
dev.off()
tiff(kDefaults$p1Tif300,width=15,height=8,units="in",res=300)
print (p2)
dev.off()
#############################
survival_data<-clin
gsva_matrix=t(ssGSEA_result)
library("survival")
library("survminer")
library("plyr")
colnames(gsva_matrix)=c("Activated.B.cell","Activated.CD4.T.cell","Activated.CD8.T.cell","Activated.dendritic.cell", "CD56bright.natural.killer.cell",
                        "CD56dim.natural.killer.cell","Central.memory.CD4.T.cell","Central.memory.CD8.T.cell","Effector.memeory.CD4.T.cell","Effector.memeory.CD8.T.cell","Eosinophil","Gamma.delta.T.cell",            
                        "Immature.B.cell","Immature.dendritic.cell","Macrophage","Mast.cell","MDSC","Memory.B.cell","Monocyte","Natural.killer.cell","Natural.killer.T.cell","Neutrophil","Plasmacytoid.dendritic.cell","Regulatory.T.cell","T.follicular.helper.cell",
                        "Type.1.T.helper.cell","Type.17.T.helper.cell","Type.2.T.helper.cell" )
survival_all_data<-merge(survival_data,gsva_matrix,by.x="sample",by.y="row.names",all=F)
survival_all_data<-na.omit(survival_all_data)
formula_for_multivariate<-as.formula(paste0('Surv(OS.time,OS)~',
                                            paste(as.character(colnames(gsva_matrix)[1:15]),sep = '',collapse = '+')))
survival_all_data$OS_TIME=as.numeric(survival_all_data$OS_TIME)
survival_all_data$OS_STATUS=gsub("DECEASED",'1',survival_all_data$OS_STATUS)
survival_all_data$OS_STATUS=gsub("LIVING",'0',survival_all_data$OS_STATUS)
survival_all_data$OS_STATUS=as.numeric(survival_all_data$OS_STATUS)
colnames(survival_all_data)[2:3]=c("OS","OS.time")
model<-coxph(formula_for_multivariate,data=survival_all_data)
p3<-ggforest(model,data=survival_all_data,main='Hazard rations of candidate genes',fontsize=1)
pdf(kDefaults$p2pdf,w=15,h=12)
print(p3)
dev.off()
tiff(kDefaults$p2Tif72,w=15,h=12,units="in",res=72)
print(p3)
dev.off()
tiff(kDefaults$p2Tif300,w=15,h=12,units="in",res=300)
print(p3)
dev.off()
#######################针对TCGA数据
###############step5:临床特征提取
clin_down<-read.csv('../../01.data/14.TCGA-LIHC.GDC_phenotype.csv',header = T)
feature<-data.frame(sample=clin_down$submitter_id.samples,clin_down)
feature$sample=gsub('A$','',feature$sample)
feature$sample=gsub('-','.',feature$sample)
aa<-RBP_cluster
infor<-merge(feature,aa,by.x="sample",by.y="sample")

###############热图##########这个针对TCGA数据因为TCGA数据有完整的临床信息，同时RBP基因用全部的基因#############
####################
RBPgenelist<-read.table('../01.数据搜集和预处理/01.RBP_genelist.txt',header = T)
RBPgenelist<-RBPgenelist[,1]
#TCGALIHC_mRNA_counts
feature<-data.frame(sample=infor$sample,cluster=infor$cluster,AGE=infor$age_at_index.demographic,
                    SEX=infor$gender.demographic)
feature$AGE<-ifelse(feature$AGE>=60,">=60","<60")
feature$sample=gsub('-','.',feature$sample)
ggdata<-tcga[intersect(RBPgenelist,rownames(tcga)),intersect(feature[,1],colnames(tcga))]
ggdata<-na.omit(ggdata) ## 1481 个基因 363 个样本 （TCGA和RBP基因交集，以及临床样本交集）
colnames(ggdata)=gsub('A$','',colnames(ggdata))
annotation_col=feature
colnames(TCGALIHC_mRNA_gene_counts)=gsub('A$','',colnames(TCGALIHC_mRNA_gene_counts))
TCGALIHC_mRNA_counts_insect=TCGALIHC_mRNA_gene_counts[intersect(RBPgenelist,rownames(tcga)),intersect(feature[,1],colnames(tcga))]##1481 个基因 363 个样本 （TCGA和RBP基因交集，以及临床样本交集）
TCGALIHC_mRNA_counts_insect=na.omit(TCGALIHC_mRNA_counts_insect)
feature=feature[match(colnames(TCGALIHC_mRNA_counts_insect),feature$sample),]
colData<-data.frame(row.names = as.factor(feature$sample),condition=as.factor(feature$cluster))
###差异表达分析##
DEseq2function<-function(datamatrix,sampleinformation,logFC_cutoff=1){
  library(DESeq2)
  dds<-DESeqDataSetFromMatrix(countData = datamatrix,
                              colData = sampleinformation,
                              design = ~ condition)
  dds<-DESeq(dds)
  ### 设置condition 两两比较
  res<- results(dds,contrast = c('condition',levels(sampleinformation$condition)))
  resOrdered <- res[order(res$pvalue),]
  DEG <- as.data.frame(resOrdered)
  ##去掉NA值
  DEG <- na.omit(DEG)
  return(DEG)
}
DEG<-DEseq2function(TCGALIHC_mRNA_counts_insect,colData)
###
logFC_cutoff<-1
DEG$change=as.factor(
  ifelse(DEG$padj<0.05 & abs(DEG$log2FoldChange)>logFC_cutoff, ifelse(DEG$log2FoldChange>logFC_cutoff,'ClusterA-UP','ClusterB-UP'),'NOSignificant')
)
DEG <- na.omit(DEG)
deseq_DEG=subset(DEG,!DEG$change %in% ('NOSignificant')) 
deseq_gene=rownames(deseq_DEG)
#########
rownames(annotation_col)<-feature[,1]
annotation_col=annotation_col[,-1]
annotation_col=annotation_col[order(annotation_col$cluster),]
sample_order<-rownames(annotation_col)

color = c(color = colorRampPalette(c("navy", "white","firebrick3"))(50))
bk<- unique(c(seq(0,5,by=0.1)))
data<-ggdata[deseq_gene,sample_order]
data<-log2(data+1)
pheatmap(data,show_rownames =F,show_colnames =F,fontsize_row=10,
         annotation_col = annotation_col,cluster_cols = F,color = color,breaks=bk) 

pdf("../../../Reuslts/Figure/Figure4_pheatmap.pdf", width=10,height=10)
pheatmap(data,show_rownames =F,show_colnames =F,fontsize_row=10,
         annotation_col = annotation_col,cluster_cols = F,color = color,breaks=bk) 
dev.off()
tiff("../../../Reuslts/Figure/Figure4_pheatmap-72.tif", width=10,height=10, units = "in", res = 72)
pheatmap(data,show_rownames =F,show_colnames =F,fontsize_row=10,
         annotation_col = annotation_col,cluster_cols = F,color = color,breaks=bk) 
dev.off()
tiff("../../../Reuslts/Figure/Figure4_pheatmap-300.tif", width=10,height=10, units = "in", res = 300)
pheatmap(data,show_rownames =F,show_colnames =F,fontsize_row=10,
         annotation_col = annotation_col,cluster_cols = F,color = color,breaks=bk) 
dev.off()
save(list=ls(),file="/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/04.output/step1-大样本.RData")







