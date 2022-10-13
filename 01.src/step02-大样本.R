#######################
rm(list = ls())
options(stringsAsFactors = F)
set.seed(123456)
library(survival)
library(maftools)
library("ComplexHeatmap")
source("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/05.src/sub/createOncoMatrix.r")
source("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/05.src/sub/oncoplot1.r")
load("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/04.output/conClust.RData")
#####
##输入数据
load("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/04.output/step1-大样本.RData")
colnames(exp)=gsub('-','.',colnames(exp))
RBP_gene_exp=data.frame(t(subset(exp,rownames(exp) %in% deseq_gene)))
RBPs_survival_genes<-merge(clin,RBP_gene_exp,by.x="sample",by.y="row.names")
RBPs_survival_genes<-na.omit(RBPs_survival_genes)
##单因素分析
library("survival")
library("survminer")
library("plyr")
RBPs_survival_genes$OS_TIME<-as.numeric(RBPs_survival_genes$OS_TIME)
RBPs_survival_genes$OS_STATUS=gsub('DECEASED',1,RBPs_survival_genes$OS_STATUS)
RBPs_survival_genes$OS_STATUS=gsub('LIVING',0,RBPs_survival_genes$OS_STATUS)
RBPs_survival_genes$OS_STATUS<-as.numeric(RBPs_survival_genes$OS_STATUS)
BaseLine<-RBPs_survival_genes
BaSurv<-Surv(time = BaseLine$OS_TIME,event = BaseLine$OS_STATUS)
UniCox<-function(x){
  FML<-as.formula(paste0('BaSurv~',x))
  GCox<-coxph(FML,data = BaseLine)
  GSum<-summary(GCox)
  HR<-round(GSum$coefficients[,2],2)
  PValue<-round(GSum$coefficients[,5],3)
  CI95<-paste0(round(GSum$conf.int[,3:4],2),collapse = '-')
  Unicox<-data.frame('Characteristics'=x,
                     'Hazard.Ratio'=HR,
                     'CI95'=CI95,
                     'P.value'=PValue)
  return(Unicox)
}
varNames<-colnames(BaseLine)[c(4:ncol(BaseLine))]
UniVar<-lapply(varNames,UniCox)
UniVar<-ldply(UniVar,data.frame)
uni_cox_df=UniVar$Characteristics[UniVar$P.value < 0.05]
colnames(RBPs_survival_genes)[2:3]=c("OS","OS.time")
select_col=c(c("sample","OS","OS.time"),uni_cox_df)
RBPs_survival_genes_cox=RBPs_survival_genes[,colnames(RBPs_survival_genes)%in% select_col]
###随机森林##
library("party")
library("randomForest")
formula_for_multivariate<-as.formula(paste0('OS.time ~',
                              paste(colnames(RBPs_survival_genes_cox)[4:ncol(RBPs_survival_genes_cox)],sep = '',collapse = '+')))
output.forest <- randomForest(formula_for_multivariate, 
                              data = RBPs_survival_genes_cox)
print(output.forest)
scores=importance(output.forest,type = 2)
##取中值以上的得分的基因
scores_select=subset(scores,scores[,1]>2000)
select_genes=rownames(scores_select)
select_genes<-c(select_genes,c("sample"))
RBPs_survival_genes_cox_randomForest=RBPs_survival_genes[,colnames(RBPs_survival_genes)%in% select_genes]
rownames(RBPs_survival_genes_cox_randomForest)=RBPs_survival_genes_cox_randomForest$sample
RBPs_survival_genes_cox_randomForest=RBPs_survival_genes_cox_randomForest[,-1]
###PCA###
library("FactoMineR")
df=RBPs_survival_genes_cox_randomForest
PCA_data <- PCA(df, graph = FALSE)
library("factoextra")
eig.val <- get_eigenvalue(PCA_data)
var <- get_pca_var(PCA_data)
pcadata=data.frame(var$coord)
##PCA分析##
#library(psych)
#df=RBPs_survival_genes_cox_randomForest
#fa.parallel(df, fa = "pc", n.iter = 100,
 #           show.legend = F, main = "Scree plot with parallel analysis")
#pc<-principal(df, nfactors = 24, score = T, rotate = "varimax")
#PCA_data=data.frame(round(unclass(pc$weights),2))
PCA_data_subset1=subset(pcadata,abs(pcadata$Dim.1)>=0.35)
PCA_data_subset2=subset(pcadata,abs(pcadata$Dim.2)>=0.35)
PCA_data_subset=unique(rbind(PCA_data_subset1,PCA_data_subset2))
###
RBP_signature_genes=rownames(PCA_data_subset)
##
RBPs_select_colum<-RBP_signature_genes
RBPs_signature_genes_fpkm=RBP_gene_exp[,colnames(RBP_gene_exp)%in% RBPs_select_colum]
write.table(RBP_signature_genes,file='/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Table/SupplementaryTable6-RBP_signature_gene.txt',sep = '\t',quote = F)
##聚类分析##
library(ConsensusClusterPlus)
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
emt<-t(RBPs_signature_genes_fpkm)
conClust <- ConsensusClusterPlus(
  as.matrix(emt),
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
save(conClust,file="/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/04.output/单因素回归分析-conClust.RData")
k<-getOptK(conClust, minCls = 2, maxCls = 5)
message(sprintf("最佳分类数：%d", k))
kDefaults <- list(
  RBPsiggroupfile = "/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Table/SupplementaryTable5_RBPsignature_clusters.txt",
  survRBPsigclusterPdf = "/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure5A_surv_RBP_cluster.pdf",  # 根据RBP分类的生存曲线
  survRBPsigclusterTif72 = "/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure5A_surv_RBP_cluster-72.tif",  # 根据RBP分类的生存曲线
  survRBPsigclusterTif300 = "/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure5A_surv_RBP_cluster-300.tif",  # 根据RBP分类的生存曲线
  boxsigclusterpdf="/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure5B_surv_RBP_cluster.pdf",
  boxsigclusterTif72="/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure5B_surv_RBP_cluster-72.tif",
  boxsigclusterTif300="/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure5B_surv_RBP_cluster-300.tif",
  alluvialDataTsv = "/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Table/SupplementaryTable7_alluvial.txt",  # 冲积图数据
  p6alluvialPdf = "/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure6A_alluvial.pdf",  # 桑基图
  p6alluvialTif72 = "/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure6A_alluvial-72.tif",  # 冲积图
  p6alluvialTif300 = "/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure6A_alluvial-300.tif",
  p6Pdf = "/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure6B_survival_score.pdf", 
  p6Tif72 = "/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure6B_survival_score-72.tif",  
  p6Tif300 = "/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure6B_survival_score-300.tif"
)
cluster1<-conClust[[k]]$consensusClass
cluster<-plyr::revalue(as.character(cluster1),c("1"="RBP_clusterA","2"="RBP_clusterB"))
RBP_sigcluster<-data.frame(sample=as.character(names(conClust[[k]]$consensusClass)),cluster=cluster)
write.table(RBP_sigcluster,file=kDefaults$RBPsiggroupfile,quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)

clin1<-clin
clin1<-clin1[which(clin1[,3]>0),]
samples<-intersect(RBP_sigcluster[,1],clin1[,1])
group<-RBP_sigcluster[match(samples,RBP_sigcluster[,1]),] #####筛选有临床信息的样本
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
  font.legend=14,xlab="OS Time (months)",legend.title="cluster",pval=TRUE,
  font.y=c(20),font.x=c(20),font.tickslab =c(18),pval.size=7,palette = c("#FF7F00","#377EB8","#DCAB3A"),risk.table = TRUE,
  ggtheme = theme_survminer()
)
#survRBPcluster
#write.table(group,file=kDefaults$m6Agroupfile,quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)
pdf(kDefaults$survRBPsigclusterPdf, width = 10, height = 8)
print(survRBPcluster)
dev.off()
tiff(kDefaults$survRBPsigclusterTif72, width = 10, height = 8, units = "in", res = 72)
print(survRBPcluster)
dev.off()
tiff(kDefaults$survRBPsigclusterTif300, width = 10, height = 8, units = "in", res = 300)
print(survRBPcluster)
dev.off()

###boxplot###
RBP_genes_fpkm=data.frame(RBPs_signature_genes_fpkm)
RBP_genes_fpkm$sample=rownames(RBP_genes_fpkm)
rownames(RBP_genes_fpkm)=NULL
datanames=RBP_sigcluster
RBP_genes_fpkm_box=merge(RBP_genes_fpkm,datanames,by="sample",all=F)
colnames(RBP_genes_fpkm_box)[1]="sample"
names=colnames(RBP_genes_fpkm_box)[2:73]
tmp_sum=data.frame(expression=1,cluster=1,gene=1)
for (name in names){
  tmp1=data.frame(RBP_genes_fpkm_box[,c(name,"cluster")])
  tmp1$gene=name
  colnames(tmp1)[1]="expression"
  tmp_sum=rbind(tmp_sum,tmp1)
  print("ok")
}
tmp_sum=tmp_sum[-1,]
p<-ggboxplot(tmp_sum, x="gene", y="expression", fill= "cluster")+
  stat_compare_means(aes(group=cluster),label = "p.signif")+xlab('')+ylab('scale batch expression')
p2<-p+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +scale_fill_manual(values=c("#FF7F00", "#377EB8","#DF342A"))
pdf(kDefaults$boxsigclusterpdf, width = 30, height = 8)
print(p2)
dev.off()
tiff(kDefaults$boxsigclusterTif72, width = 30, height = 8, units = "in", res = 72)
print(p2)
dev.off()
tiff(kDefaults$boxsigclusterTif300, width = 30, height = 8, units = "in", res = 300)
print(p2)
dev.off()
###############PCA计算得分####################
##PCA分析##
library("FactoMineR")
df=data.frame(t(RBPs_signature_genes_fpkm))
PCA_data <- PCA(df, graph = FALSE)
library("factoextra")
eig.val <- get_eigenvalue(PCA_data)
var <- get_pca_var(PCA_data)
pcadata=data.frame(var$coord)
pcadata$score=pcadata$Dim.1+pcadata$Dim.2
pcadata=pcadata[order(pcadata$score),]
RBP_score=pcadata
RBP_score$class_score=as.factor(
  ifelse(RBP_score$score > median(RBP_score$score),'high','low')
)
write.table(RBP_score,file="/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Table/SupplementaryTable6_RBPscore.txt",sep="\t",quote = F)
############
RBP_clusters<-read.table('/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Table/SupplementaryTable2_RBP_cluster.txt',header = T)
RBP_clusters<-na.omit(RBP_clusters)
RBP_clusters$sample=gsub('-','.',RBP_clusters$sample)
RBP_sig_clusters<-read.table('/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Table/SupplementaryTable5_RBPsignature_clusters.txt',header = T)
RBP_sig_clusters$sample=gsub('-','.',RBP_sig_clusters$sample)
tmp1=merge(RBP_clusters,RBP_sig_clusters,by="sample",all=F)
tmp2=merge(tmp1,RBP_score,by.x="sample",by.y="row.names",all=F)
sample_all=merge(tmp2,clin1,by="sample",all=F)
colnames(sample_all)=c("sample","RBP_subtype","RBP_signature_subtype","Dim.1","Dim.2","Dim.3","Dim.4","Dim.5","score","class_score","OS.status","OS.time")

library(ggalluvial)
library(ggplot2)
library(dplyr)
write.table(sample_all,file = "/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Table/SupplementaryTable7_RBPclusters_pic.txt",sep = '\t',quote = F)
alluvial<-ggplot(
  as.data.frame(sample_all),
  aes(
    axis1 = RBP_subtype, 
    axis2 = RBP_signature_subtype,
    axis3 = class_score,
    axis4 = OS.status,
    y =RBP_subtype
  )
) + scale_x_discrete(
  limits = c("RBP_subtype","RBP_signature_subtype", "class_score","OS.status"),
  expand = c(.1, .05)
) + geom_alluvium(
  aes(fill = as.factor(RBP_subtype)),
  show.legend = FALSE
) + geom_stratum(
) + geom_text(
  stat = "stratum",
  aes(label = after_stat(stratum))
) + theme_minimal()+
ggtitle("Patients in the LIHC cohort",
          "stratified by RBP genes and survival")
pdf(kDefaults$p6alluvialPdf, width = 8, height = 8)
print(alluvial)
dev.off()
tiff(kDefaults$p6alluvialTif72, width = 8, height = 8, units = "in", res = 72)
print(alluvial)
dev.off()
tiff(kDefaults$p6alluvialTif300, width = 8, height = 8, units = "in", res = 300)
print(alluvial)
dev.off()

sample_pic=sample_all
sample_pic$OS.status=gsub('DECEASED','1',sample_pic$OS.status)
sample_pic$OS.status=gsub('LIVING','0',sample_pic$OS.status)
sample_pic$OS.status=as.numeric(sample_pic$OS.status)
sample_pic$OS.time=as.numeric(sample_pic$OS.time)
fit_km<-survfit(Surv(OS.time,OS.status)~class_score,data=sample_pic)
p<-ggsurvplot(fit_km,conf.int=F,pval=T,legend.title='RBP clusters',risk.table=T, palette=c('red','blue',"green","pink"),surv.median.line='hv')
pdf(kDefaults$p6Pdf, width = 8, height = 8)
print(p)
dev.off()
tiff(kDefaults$p6Tif72, width = 8, height = 8, units = "in", res = 72)
print(p)
dev.off()
tiff(kDefaults$p6Tif300, width = 8, height = 8, units = "in", res = 300)
print(p)
dev.off()

################针对TCGA的数据##############
feature<-read.csv('/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/01.data/14.TCGA-LIHC.GDC_phenotype.csv',header = T)
colnames(feature)[1]='sample'
feature$sample=gsub('-','.',feature$sample)
feature$sample=gsub('A$','',feature$sample)
aa<-data.frame(sample=sample_all$sample,sample_all)
infor<-merge(feature,aa,by.x="sample",by.y="sample")
#######################针对TCGA数据
###############step3:临床特征分析
gg2<-data.frame(pathologic_M=infor$pathologic_M,group=infor$class_score,score=infor$score)
p1<-ggboxplot(data = gg2,x="pathologic_M" , y = "score", fill="group")+stat_compare_means(aes(group = group))
print(p1)
pdf("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure7A_pathologic_M.pdf",width=8,height=8)
print (p1)
dev.off()
tiff("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure7A_pathologic_M-72.tif",width=8,height=8,units="in",res=72)
print (p1)
dev.off()
tiff("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure7A_pathologic_M-300.tif",width=8,height=8,units="in",res=300)
print (p1)
dev.off()

###########person_neoplasm_cancer_status########
gg2<-data.frame(person_neoplasm_cancer_status=infor$person_neoplasm_cancer_status,cluster=infor$class_score,score=infor$score)
gg2=subset(gg2,gg2$person_neoplasm_cancer_status%in% c('TUMOR FREE','WITH TUMOR'))
p2<-ggboxplot(data = gg2,x="person_neoplasm_cancer_status" , y = "score", fill="cluster")+stat_compare_means(aes(group = cluster))
print(p2)
pdf("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure7B_person_neoplasm_cancer_status.pdf",width=8,height=8)
print (p2)
dev.off()
tiff("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure7B_person_neoplasm_cancer_status-72.tif",width=8,height=8,units="in",res=72)
print (p2)
dev.off()
tiff("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure7B_person_neoplasm_cancer_status-300.tif",width=8,height=8,units="in",res=300)
print (p2)
dev.off()
##################
gg2<-data.frame(SEX=infor$gender.demographic,cluster=infor$class_score,score=infor$score)
p3<-ggboxplot(data = gg2,x="SEX" , y = "score", fill="cluster")+stat_compare_means(aes(group = cluster))
print(p3)
pdf("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure7C_SEX.pdf",width=8,height=8)
print (p3)
dev.off()
tiff("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure7C_SEX-72.tif",width=8,height=8,units="in",res=72)
print (p3)
dev.off()
tiff("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure7C_SEX-300.tif",width=8,height=8,units="in",res=300)
print (p3)
dev.off()
############race.demographic
gg2<-data.frame(race.demographic=infor$race.demographic,cluster=infor$class_score,score=infor$score)
p4<-ggboxplot(data = gg2,x="race.demographic" , y = "score", fill="cluster")+stat_compare_means(aes(group = cluster))
pdf("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure7D_race.pdf",width=12,height=8)
print (p4)
dev.off()
tiff("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure7D_race-72.tif",width=12,height=8,units="in",res=72)
print (p4)
dev.off()
tiff("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure7D_race-300.tif",width=12,height=8,units="in",res=300)
print (p4)
dev.off()
########neoplasm_histologic_grade
gg2<-data.frame(neoplasm_histologic_grade=infor$neoplasm_histologic_grade,cluster=infor$class_score,score=infor$score)
gg2=subset(gg2,gg2$neoplasm_histologic_grade %in% c('G1','G2','G3','G4'))
p5<-ggboxplot(data = gg2,x="neoplasm_histologic_grade" , y = "score", fill="cluster")+stat_compare_means(aes(group = cluster))
pdf("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure7E_Grade.pdf",width=8,height=8)
print (p5)
dev.off()
tiff("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure7E_Grade-72.tif",width=8,height=8,units="in",res=72)
print (p5)
dev.off()
tiff("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure7E_Grade-300.tif",width=8,height=8,units="in",res=300)
print (p5)
dev.off()
############age_at_initial_pathologic_diagnosis
gg2<-data.frame(AGE=infor$age_at_initial_pathologic_diagnosis,cluster=infor$class_score,score=infor$score)
gg2$AGE<-ifelse(gg2$AGE>=60,">=60","<60")
p6<-ggboxplot(data = gg2,x="AGE" , y = "score", fill="cluster")+stat_compare_means(aes(group = cluster))
print(p6)
pdf("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure7F_AGE.pdf",width=8,height=8)
print (p6)
dev.off()
tiff("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure7F_AGE-72.tif",width=8,height=8,units="in",res=72)
print (p6)
dev.off()
tiff("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure7F_AGE-300.tif",width=8,height=8,units="in",res=300)
print (p6)
dev.off()

############生存
data=sample_all
index<-grep("TCGA",data[,1])
survRBPScoreData1 <- data.frame(
  time = as.numeric(data$OS.time[index]),
  status = as.numeric(data$OS.status[index]),
  RBPscoreGroup = data$class_score[index]
)

survRBPScore1 <- ggsurvplot(
  survfit(Surv(time, status) ~ RBPscoreGroup, data = survRBPScoreData1),
  font.legend=14,xlab="OS Time (months)",legend.title="TCGA cluster",pval=TRUE,
  font.y=c(20),font.x=c(20),font.tickslab =c(18),pval.size=7,palette = c("#FF7F00","#377EB8"),risk.table = TRUE,
  ggtheme = theme_survminer()
)

pdf("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure7G.pdf", width = 12, height = 8)
print(survRBPScore1)
dev.off()
tiff("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure7G-72.tif", width = 12, height = 8, units = "in", res = 72)
print(survRBPScore1)
dev.off()
tiff("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure7G-300.tif", width = 12, height = 8, units = "in", res = 300)
print(survRBPScore1)
dev.off()

##############
########################
###########2,突变信息展示
library("DT")
library(TCGAbiolinks)
query.maf.hg19 <- GDCquery(project = "TCGA-LIHC", 
                           data.category = "Simple nucleotide variation", 
                           data.type = "Simple somatic mutation",
                           access = "open", 
                           file.type = "LIHC_pairs.aggregated.capture.tcga.uuid.automated.somatic.maf",
                           legacy = TRUE)
GDCdownload(query.maf.hg19,files.per.chunk = 200)
maf <- GDCprepare(query.maf.hg19)
SNP_samples_LIHC<-unlist(lapply(unique(maf[,16]$Tumor_Sample_Barcode),function(b){substr(b,1,15)}))
count<-sort(table(SNP_samples_LIHC))
index<-which(count>1)
if(length(index)>1){
  rm_sample<-names(count)[index]
  fullname<-unique(maf[,16]$Tumor_Sample_Barcode)
  index_rm<-c()
  for(i in 1:length(rm_sample)){
    index<-grep(rm_sample[i],fullname)
    index_rm<-c(index_rm,index[2:length(index)])
  }
  retain_samples<-fullname[-index_rm]
  index<-which(maf[,16]$Tumor_Sample_Barcode %in% retain_samples)
  maf<-maf[index,]
  SNP_samples_LIHC<-unlist(lapply(maf[,16]$Tumor_Sample_Barcode,function(b){substr(b,1,15)})) }###413
sort(table(SNP_samples_LIHC))
write.table(maf,file="/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Table/hg19_somatic_maf.maf",quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)


laml = read.maf(maf = "/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Table/hg19_somatic_maf.maf")
SNP_samples<-laml@data$Tumor_Sample_Barcode
SNP_samples<-unlist(lapply(SNP_samples,function(b){substr(b,1,15)}))
SNP_samples=gsub('-','.',SNP_samples)
mutation_count<-getSampleSummary(laml)
mutation<-mutation_count$total
Mut<-data.frame(mutation_count[,1],Mutation=mutation)
Mut[,1]<-unlist(lapply(Mut[,1],function(b){substr(b,1,15)}))
Mut$Tumor_Sample_Barcode=gsub('-','.',Mut$Tumor_Sample_Barcode)
clin_anno<-merge(infor,Mut,by.x="sample",by.y="Tumor_Sample_Barcode")

high_samples<-as.character(clin_anno[clin_anno$class_score=="high",1])
low_samples<-as.character(clin_anno[clin_anno$class_score=="low",1])


maf1<-laml@data[which(SNP_samples %in% high_samples),]
write.table(maf1,file="/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Table/hg19_somatic_highsamples.maf",quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)
laml_high = read.maf(maf = "/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Table/hg19_somatic_highsamples.maf")
genes = getGeneSummary(x = laml_high)[1:nrow(getGeneSummary(x = laml_high)), Hugo_Symbol]
om = createOncoMatrix(m = laml_high, g = genes)
snp_matrix_high = om$oncoMatrix
colnames(snp_matrix_high)<-unlist(lapply(colnames(snp_matrix_high),function(b){substr(b,1,15)}))
top_gene_high <- getGeneSummary(x = laml_high)[1:20, Hugo_Symbol]
#################### low sample
maf2<-laml@data[which(SNP_samples %in% low_samples),]
write.table(maf2,file="/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Table/hg19_somatic_lowsamples.maf",quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)
laml_low = read.maf(maf = "/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Table/hg19_somatic_lowsamples.maf")
genes = getGeneSummary(x = laml_low)[1:nrow(getGeneSummary(x = laml_low)), Hugo_Symbol]
om = createOncoMatrix(m = laml_low, g = genes)
snp_matrix_low = om$oncoMatrix
colnames(snp_matrix_low)<-unlist(lapply(colnames(snp_matrix_low),function(b){substr(b,1,15)}))
top_gene_low <- getGeneSummary(x = laml_low)[1:20, Hugo_Symbol]
#################################################
col<-get_vcColors()
ovlp<-intersect(top_gene_high,top_gene_low)
length(ovlp) ##13
###########################high figure15_1
rownames(clin_anno)<-clin_anno[,1]
colnames(snp_matrix_high)<-unlist(lapply(colnames(snp_matrix_high),function(b){substr(b,0,15)}))
colnames(snp_matrix_high)<-gsub('-','.',colnames(snp_matrix_high))
anno<-clin_anno[match(colnames(snp_matrix_high),rownames(clin_anno)),]
for (i in 1:nrow(anno)){
  if (anno$Mutation[i]>500){
      anno$Mutation[i]=500
  }
}
p1<-oncoPrint(snp_matrix_high[ovlp,],alter_fun = alter_fun, width = 133,
              col = col, remove_empty_columns = FALSE,row_order = ovlp,
              top_annotation = HeatmapAnnotation(gap = unit(0.5, "mm"),height = unit(60, "mm"),
                                                 Mutation = anno_barplot(as.numeric(anno$Mutation),border = F,axis = T,gp = gpar(fill = "honeydew4")),
                                                 #bar4 = anno_barplot(anno$TMB,border = F,axis = T,gp = gpar(fill = "honeydew4")),
                                                 RBPscore = anno_barplot(anno$score,border = F,axis = T,gp = gpar(fill = "honeydew4")),
                                                 RBPcoreGroup=anno$class_score,
                                                 #Stage=anno$stage,
                                                 OS_status=anno$OS.status,
                                                 Gender=anno$gender.demographic,
                                                 pathologic_M=anno$pathologic_M,
                                                 col=list(
                                                   RBPcoreGroup=c('high'='#FF7F00','low'='#377EB8'),
                                                   OS_status=c('DECEASED'='firebrick4','LIVING'='darkolivegreen3'),
                                                   Gender=c('female'='red','male'='blue'),
                                                   pathologic_M=c('M0'='red','MX'='yellow','M1'='blue')
                                                   #Stage=c('Stage I'='#FF7F00FF','Stage II'='firebrick2','Stage III'='firebrick','Stage IV'='firebrick4','Stage X'='SlateGray','Not Available'='grey')
                                                 )                               
              ))

colnames(snp_matrix_low)<-unlist(lapply(colnames(snp_matrix_low),function(b){substr(b,0,16)}))
colnames(snp_matrix_low)<-gsub('-','.',colnames(snp_matrix_low))
anno<-clin_anno[match(colnames(snp_matrix_low),rownames(clin_anno)),]
for (i in 1:nrow(anno)){
  if (anno$Mutation[i]>500){
    anno$Mutation[i]=500
  }
}
p2<-oncoPrint(snp_matrix_low[ovlp,],alter_fun = alter_fun, width = 133,
              col = col, remove_empty_columns = FALSE,row_order = ovlp,
              top_annotation = HeatmapAnnotation(gap = unit(0.5, "mm"),height = unit(60, "mm"),
                                                 Mutation = anno_barplot(as.numeric(anno$Mutation),border = F,axis = T,gp = gpar(fill = "honeydew4")),
                                                 RBPscore = anno_barplot(anno$score,border = F,axis = T,gp = gpar(fill = "honeydew4")),
                                                 RBPcoreGroup=anno$class_score,
                                                 OS_status=anno$OS.status,
                                                 Gender=anno$gender.demographic,
                                                 pathologic_M=anno$pathologic_M,
                                                 col=list(
                                                   RBPcoreGroup=c('high'='#FF7F00','low'='#377EB8'),
                                                   OS_status=c('DECEASED'='firebrick4','LIVING'='darkolivegreen3'),
                                                   Gender=c('female'='red','male'='blue'),
                                                   pathologic_M=c('M0'='red','MX'='yellow','M1'='blue')
                                                   #Stage=c('Stage I'='#FF7F00FF','Stage II'='firebrick2','Stage III'='firebrick','Stage IV'='firebrick4','Stage X'='SlateGray','Not Available'='grey')
                                                 )                               
              ))

kDefaults=list(
  p7Pdf="/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure8A.pdf",
  p7Tif72='/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure8A-72.tif',
  p7Tif300='/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure8A-300.tif',
  p8Pdf="/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure8B.pdf",
  p8Tif72='/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure8B-72.tif',
  p8Tif300='/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure8B-300.tif'
)

pdf(kDefaults$p7Pdf,width=10,height=6)
print (p1)
dev.off()
tiff(kDefaults$p7Tif72,width=10,height=6,units="in",res=72)
print (p1)
dev.off()
tiff(kDefaults$p7Tif300,width=10,height=6,units="in",res=300)
print (p1)
dev.off()


pdf(kDefaults$p8Pdf,width=10,height=6)
print (p2)
dev.off()
tiff(kDefaults$p8Tif72,width=10,height=6,units="in",res=72)
print (p2)
dev.off()
tiff(kDefaults$p8Tif300,width=10,height=6,units="in",res=300)
print (p2)
dev.off()

save(list=ls(),file="/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/04.output/tmp2.RData")


##############需要修改#########
##############3,CNV gistic提示有overlap，是因为有样本重复
seg = read.table("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/01.data/17.TCGA-LIHC.masked_cnv.tsv",header=TRUE,sep="\t",as.is=TRUE)
a<-sort(table(paste(seg[,1],seg[,2],seg[,3],sep="_")),decreasing = TRUE)
tt<-names(a)[which(a==2)]
rm<-names(table(unlist(lapply(tt,function(b){strsplit(b,split="_")[[1]][1]}))))
aa<-data.frame(seg,over=paste(seg[,1],seg[,2],seg[,3],sep="_"))
bb<-data.frame(count=table(paste(seg[,1],seg[,2],seg[,3],sep="_")),over=names(table(paste(seg[,1],seg[,2],seg[,3],sep="_"))))
tt<-merge(aa,bb,by.x="over",by.y="over")
index<-which(tt$count.Freq>1)
seg<-tt[-index,c(2,3,4,5,6)]

index<-grep("01A$",seg[,1])
seg<-seg[index,]
seg[,1]<-unlist(lapply(seg[,1],function(b){substr(b,0,15)}))
sort(table(unlist(lapply(unique(seg[,1]),function(b){substr(b,0,15)}))))


high_samples<-rownames(RBP_score)[which(RBP_score$class_score=="high")]
low_samples<-rownames(RBP_score)[which(RBP_score$class_score=="low")]
seg$sample=gsub('-','.',seg$sample)
high_seg<-seg[which(seg[,1] %in% high_samples),]
low_seg<-seg[which(seg[,1] %in% low_samples),]
write.table(high_seg,file="/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Table/SupplementaryTable8_high_seg-v1.txt",quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)
write.table(low_seg,file="/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Table/SupplementaryTable8_low_seg-v1.txt",quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)



library(maftools)
all.lesions <- '/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/04.output/GISTIC/example_results_TCGA_high/all_lesions.conf_95.txt'
amp.genes <- '/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/04.output/GISTIC/example_results_TCGA_high/amp_genes.conf_95.txt'
del.genes <- '/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/04.output/GISTIC/example_results_TCGA_high/del_genes.conf_95.txt'
scores.gis <- '/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/04.output/GISTIC/example_results_TCGA_high/scores.gistic'
laml.gistic = readGistic(gisticAllLesionsFile = all.lesions, gisticAmpGenesFile = amp.genes, gisticDelGenesFile = del.genes, gisticScoresFile = scores.gis)
pdf("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure8C_RBP_high_cnv.pdf",width=10,height=6)
gisticChromPlot(gistic = laml.gistic)  
dev.off()

tiff("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure8C_RBP_high_cnv-72.tif",width=10,height=6,units="in",res=72)
gisticChromPlot(gistic = laml.gistic)  
dev.off()

tiff("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure8C_RBP_high_cnv-300.tif",width=10,height=6,units="in",res=300)
gisticChromPlot(gistic = laml.gistic)  
dev.off()


all.lesions <- '/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/04.output/GISTIC/example_results_TCGA_low/all_lesions.conf_95.txt'
amp.genes <- '/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/04.output/GISTIC/example_results_TCGA_low/amp_genes.conf_95.txt'
del.genes <- '/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/04.output/GISTIC/example_results_TCGA_low/del_genes.conf_95.txt'
scores.gis <- '/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/04.output/GISTIC/example_results_TCGA_low/scores.gistic'
laml.gistic = readGistic(gisticAllLesionsFile = all.lesions, gisticAmpGenesFile = amp.genes, gisticDelGenesFile = del.genes, gisticScoresFile = scores.gis, isTCGA = TRUE)
#gisticBubblePlot(gistic = laml.gistic,fdrCutOff=0.01)
pdf("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure8D_RBP_low_cnv.pdf",width=10,height=6)
gisticChromPlot(gistic = laml.gistic)  
dev.off()

tiff("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure8D_RBP_low_cnv-72.tif",width=10,height=6,units="in",res=72)
gisticChromPlot(gistic = laml.gistic)  
dev.off()

tiff("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure8D_RBP_low_cnv-300.tif",width=10,height=6,units="in",res=300)
gisticChromPlot(gistic = laml.gistic)  
dev.off()


library(survival, quietly = TRUE)
library(survminer, quietly = TRUE)
library(ggplot2)
library(ggpubr)
kDefaults <- list(
  RBPScore = "/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Table/SupplementaryTable6_RBPscore.txt",
  p1Pdf = "/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure9A_ic50.pdf",  ###TCGA KM 
  p1Tif72 = "/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure9A_ic50-72.tif",  ###TCGA KM 
  p1Tif300 = "/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure9A_ic50-300.tif",
  p2Pdf = "/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure9B_ic50_2.pdf",  ###TCGA KM 
  p2Tif72 = "/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure9B_ic50_2-72.tif",  ###TCGA KM 
  p2Tif300 = "/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure9B_ic50_2-300.tif",
  
  ic50="/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Table/SupplementaryTable11_ic50.csv",###TCGA KM
  #tide="./result/Table/SupplementaryTable19_tide.txt",
  p3Pdf = "/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure9C_tide.pdf",  ###ROC
  p3Tif72 = "/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure9C_tide-72.tif",  ###ROC
  p3Tif300 = "/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure9C_tide-300.tif"
)

library("pROC")
RBPscore<-read.table(file=kDefaults$RBPScore,header=TRUE,sep="\t",as.is=TRUE,row.names=1)
library(pRRophetic)
library(ggplot2)
set.seed(123)

da<-exp
pRRopheticQQplot("Gemcitabine")
cvOut <- pRRopheticCV("Gemcitabine", cvFold=10,testExprData=as.matrix(da))
summary(cvOut)
predictedPtype1 <- pRRopheticPredict(as.matrix(da), "Cisplatin", selection = 1)
predictedPtype2 <- pRRopheticPredict(as.matrix(da), "Gemcitabine", selection = 1)


library(ggpubr)
predictedPtype1<-predictedPtype1[rownames(RBPscore)]
predictedPtype2<-predictedPtype2[rownames(RBPscore)]
dat<-data.frame('RBPscoreGroup'=RBPscore$class_score,'IC50'=predictedPtype1) #Cisplatin
dat1<-data.frame('RBPscoreGroup'=RBPscore$class_score,'IC50'=predictedPtype2) #Gemcitabine
#ggscatter(dat,x='Class',y='IC50')
p1<-ggviolin(dat,x='RBPscoreGroup',y='IC50',fill='RBPscoreGroup',title = 'Cisplatin',
             palette=c("indianred2","turquoise3"),add.params=list(fill="white"),add="boxplot")+stat_compare_means()

p2<-ggviolin(dat1,x='RBPscoreGroup',y='IC50',fill='RBPscoreGroup',title = 'Gemcitabine',
             palette=c("indianred2","turquoise3"),add.params=list(fill="white"),add="boxplot")+stat_compare_means()

write.csv(cbind(RBPscore,'Cisplatin'=dat$IC50,'Gemcitabine'=dat1$IC50),kDefaults$ic50,quote = F)

pdf(kDefaults$p1Pdf, width = 8, height = 8)
print(p1)
dev.off()
tiff(kDefaults$p1Tif72, width = 8, height = 8, units = "in", res = 72)
print(p1)
dev.off()
tiff(kDefaults$p1Tif300, width = 8, height = 8, units = "in", res = 300)
print(p1)
dev.off()

pdf(kDefaults$p2Pdf, width = 8, height = 8)
print(p2)
dev.off()
tiff(kDefaults$p2Tif72, width = 8, height = 8, units = "in", res = 72)
print(p2)
dev.off()
tiff(kDefaults$p2Tif300, width = 8, height = 8, units = "in", res = 300)
print(p2)
dev.off()













