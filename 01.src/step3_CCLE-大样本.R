#######################
rm(list = ls())
options(stringsAsFactors = F)
set.seed(123456)
library(survival)
#library(Boruta, quietly = TRUE)
library(maftools)
library("ComplexHeatmap")
source("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/05.src/sub/createOncoMatrix.r")
source("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/05.src/sub/oncoplot1.r")
load("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/04.output/conClust.RData")
#####################
data<-read.table(file="/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/01.data/18.CCLE_RNAseq_genes_rpkm_20180929.gct",skip = 2, header = FALSE, sep = "\t",as.is=TRUE)
RBP_signature_genes<-read.table('/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Table/SupplementaryTable6-RBP_signature_gene.txt',header=T)
pp<-data[-1,-c(1,2)]
colnames(pp)<-data[1,][-c(1,2)]
ga<-data[,2][-1]
index<-which(ga %in% as.character(RBP_signature_genes[,1]))
pp<-pp[index,]
rownames(pp)<-ga[index]
infor<-read.table(file="/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/01.data/19.CCLE_sample_info_file_2012-10-18.txt",header=T,sep="\t",as.is=TRUE) #1046
infor<-infor[which(infor[,6] %in% "carcinoma"),]
celline<-intersect(colnames(pp),infor[,1]) #592
index<-which(colnames(pp) %in% celline)
ccle_exp<-as.matrix(pp[,index])
write.table(ccle_exp,file="/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/01.data/CCLE/aa.txt",quote=FALSE,sep="\t",col.names=TRUE,row.names=TRUE)
ccle_exp<-read.table(file="/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/01.data/CCLE/aa.txt",header=TRUE,sep="\t",as.is=TRUE)
##PCA分析##
library("FactoMineR")
df=data.frame(ccle_exp)
PCA_data <- PCA(df, graph = FALSE)
library("factoextra")
eig.val <- get_eigenvalue(PCA_data)
var <- get_pca_var(PCA_data)
pcadata=data.frame(var$coord)
pcadata$score=pcadata$Dim.1+pcadata$Dim.2
pcadata=pcadata[order(pcadata$score),]
score_ccle=pcadata
score_ccle$class_score=as.factor(
    ifelse(score_ccle$score > median(score_ccle$score),'high','low')
)

high<-rownames(score_ccle)[score_ccle$class_score=='high']
low<-rownames(score_ccle)[score_ccle$class_score=='low']
write.table(score_ccle,file="/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Table/SupplementaryTable10_score_ccle.txt",quote=FALSE,sep="\t",col.names=FALSE,row.names=TRUE)

################ 突变信息
library("maftools")
laml<-read.table("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/01.data/20.CCLE_DepMap_18q3_maf_20180718.txt",header=TRUE,as.is=TRUE,sep="\t")
laml2<-cbind(laml[,c(1:11)],Tumor_Seq_Allele2="NA",laml[,c(12:33)])
write.table(laml2,file="/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/01.data/CCLE/CCLE.maf",quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)
laml3<-read.maf("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/01.data/CCLE/CCLE.maf")
genes = getGeneSummary(x = laml3)[1:nrow(getGeneSummary(x = laml3)), Hugo_Symbol]
om_all = createOncoMatrix(m = laml3, g = genes)
mut<-om_all$numericMatrix
mut[mut>0]<-1
mut<-apply(mut,2,sum)
mut<-data.frame(sample=names(mut),mut)
type<-data.frame(sample=c(high,low),group=c(rep("high",length(high)),rep("low",length(low))))
gg<-merge(mut,type,by="sample")
write.table(gg,file="/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Table/SupplementaryTable11_mutCount.txt",quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)
for (i in 1: nrow(gg)){
    if (gg$mut[i]>2000){
        gg$mut[i]=2000
    }
}
p<-ggboxplot(gg,x="group",y="mut",color="group",add="jitter")
#CR  PR  SD  PD
my_comparisons<-list(c("high","low"))
p2<-p+stat_compare_means(comparisons=my_comparisons,method="wilcox.test")+ylab("Mutation count")

pdf("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure9C_Mutation.pdf",width=8,height=8)
print (p2)
dev.off()
tiff("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure9C_Mutation-72.tif",width=8,height=8,units="in",res=72)
print (p2)
dev.off()
tiff("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure9C_Mutation-300.tif",width=8,height=8,units="in",res=300)
print (p2)
dev.off()



SNP_samples<-laml3@data$Tumor_Sample_Barcode
maf1<-laml3@data[which(SNP_samples %in% high),]
write.table(maf1,file="/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Table/hg19_somatic_highsamples_ccle.maf",quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)
laml_high = read.maf(maf = "/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Table/hg19_somatic_highsamples_ccle.maf")
genes = getGeneSummary(x = laml_high)[1:nrow(getGeneSummary(x = laml_high)), Hugo_Symbol]
om = createOncoMatrix(m = laml_high, g = genes)
snp_matrix_high = om$oncoMatrix
#colnames(snp_matrix_high)<-unlist(lapply(colnames(snp_matrix_high),function(b){substr(b,1,15)}))
top_gene_high <- getGeneSummary(x = laml_high)[1:20, Hugo_Symbol]
#################### low sample
maf2<-laml3@data[which(SNP_samples %in% low),]
write.table(maf2,file="/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Table/hg19_somatic_lowsamples_ccle.maf",quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)
laml_low = read.maf(maf = "/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Table/hg19_somatic_lowsamples_ccle.maf")
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
oncoplot(maf = laml_high, top = 20)  

pdf("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure9A_ccle_high.pdf",width=10,height=6)
oncoplot(maf = laml_high, top = 20)  
dev.off()
tiff("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure9A_ccle_high-72.tif",width=10,height=6,units="in",res=72)
oncoplot(maf = laml_high, top = 20)  
dev.off()
tiff("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure9A_ccle_high-300.tif",width=10,height=6,units="in",res=300)
oncoplot(maf = laml_high, top = 20)  
dev.off()


pdf("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure9B_ccle_low.pdf",width=10,height=6)
oncoplot(maf = laml_low, top = 20)  
dev.off()
tiff("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure9B_ccle_low-72.tif",width=10,height=6,units="in",res=72)
oncoplot(maf = laml_low, top = 20)  
dev.off()
tiff("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure9B_ccle_low-300.tif",width=10,height=6,units="in",res=300)
oncoplot(maf = laml_low, top = 20)  
dev.off()

##########
#####################IMvigor210CoreBiologies. install.packages("/Users/lily/Desktop/IMvigor210CoreBiologies_1.0.0.tar.gz", repos=NULL)
library("IMvigor210CoreBiologies")
data(cds)
counts<-counts(cds)
fdata<-fData(cds)
pdata<-pData(cds)
countToFpkm <- function(counts, effLen)
{
    N <- sum(counts)
    exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}
fpkm<-countToFpkm(counts,fdata[,4])
fpkm_log<-log(fpkm+1,base=2)
fpkm_log<-t(sapply(split(as.data.frame(fpkm_log),fdata[,2]),function(x) colMeans(x[,1:ncol(x)])))
RBP_signature_genes<-read.table('/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Table/SupplementaryTable6-RBP_signature_gene.txt',header = T)
RBP_signature_dataset<-subset(fpkm_log,rownames(fpkm_log)%in%RBP_signature_genes[,1])
##PCA分析##
library("FactoMineR")
df=data.frame(RBP_signature_dataset)
PCA_data <- PCA(df, graph = FALSE)
library("factoextra")
eig.val <- get_eigenvalue(PCA_data)
var <- get_pca_var(PCA_data)
pcadata=data.frame(var$coord)
pcadata$score=pcadata$Dim.1+pcadata$Dim.2
pcadata=pcadata[order(pcadata$score),]
score_data=pcadata
score_data$class_score=as.factor(
    ifelse(score_data$score > median(score_data$score),'high','low')
)

high<-rownames(score_data)[score_data$class_score=='high']
low<-rownames(score_data)[score_data$class_score=='low']
RBPScore_vigor=score_data
################生存曲线
rownames(pdata)==colnames(fpkm_log)
clin_vigor<-pdata[,c(21:22)]
survm6AScoreData_vigor <- data.frame(
    time = clin_vigor$os,
    status = as.numeric(clin_vigor$censOS),
    RBPscore = RBPScore_vigor
)
res_cut <- surv_cutpoint(survm6AScoreData_vigor, time = "time", event = "status",variables = "RBPscore.score")      
cutoff1<-summary(res_cut)$cutpoint  ## 1.14
#cutoff<-median(survm6AScoreData$m6Ascore)
#survm6AScoreData_vigor$metabolicscoreGroup = ifelse(
#    survm6AScoreData_vigor$metabolicscore > cutoff1,
 #   "high", "low"
#)
write.table(survm6AScoreData_vigor, file="/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Table/SupplementaryTable12_RBPScore_IMvigor210CoreBiologies.txt", sep = "\t", col.names = NA)

# 绘制生存曲线
survm6AScore <- ggsurvplot(
    survfit(Surv(time, status) ~ RBPscore.class_score, data = survm6AScoreData_vigor),
    font.legend=14,xlab="OS Time (years)",legend.title="cluster",pval=TRUE,
    font.y=c(20),font.x=c(20),font.tickslab =c(18),pval.size=7,palette = c("#FF7F00","#377EB8"),risk.table = TRUE,
    ggtheme = theme_survminer()
)
pdf("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure11A_survRBPScorecluster_IMvigor210CoreBiologies.pdf", width = 12, height = 8)
print(survm6AScore)
dev.off()
tiff("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure11A_survRBPScorecluster_IMvigor210CoreBiologies-72.tif", width = 12, height = 8, units = "in", res = 72)
print(survm6AScore)
dev.off()
tiff("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure11A_survmRBPScorecluster_IMvigor210CoreBiologies-300.tif", width = 12, height = 8, units = "in", res = 300)
print(survm6AScore)
dev.off()

table(rownames(survm6AScoreData_vigor)==rownames(pdata))
aa<-cbind(survm6AScoreData_vigor,pdata)
write.table(aa, file="/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Table/SupplementaryTable12_infor_IMvigor210CoreBiologies.txt", sep = "\t", col.names = NA)

gg2<-aa[,c(10,9)]
count<-gg2 %>% table()
ss<-apply(count,2,function(b){b/sum(b)})
chisq.test(count)[[3]][1] ##0.0353723

ggdata<-reshape2::melt(ss, varnames = c("overall_response", "cluster"))
ggdata[,1]<-as.factor(ggdata[,1])
p2<-ggplot(data = ggdata, mapping = aes(x = cluster, y = value, fill=overall_response)) + geom_bar(stat = 'identity', position = 'stack') +
       ylab("percentage of samples(%)")+xlab("")+theme_bw() + theme(text=element_text(size=15))+
       theme(axis.line.x = element_line(colour = "black"))#+scale_fill_manual(values=c("#FF7F00", "#377EB8"))
pdf("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure11B_response.pdf",width=8,height=8)
print (p2)
dev.off()
tiff("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure11B_response-72.tif",width=8,height=8,units="in",res=72)
print (p2)
dev.off()
tiff("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure11B_response-300.tif",width=8,height=8,units="in",res=300)
print (p2)
dev.off()

library("ggpubr")
df<-aa[-which(aa[,10] %in% "NE"),]
colnames(df)[10]<-"Overall_Response"
df[,10]<-as.factor(df[,10])
#df$Overall_Response<-factor(df$Overall_Response) 
p<-ggboxplot(df,x="Overall_Response",y="RBPscore.score",color="Overall_Response",add="jitter")
#CR  PR  SD  PD
my_comparisons<-list(c("CR","PR"),c("CR","SD"),c("CR","PD"),c("PR","SD"),c("PR","PD"),c("SD","PD"))
p2<-p+stat_compare_means(comparisons=my_comparisons,method="t.test")

pdf("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure11C_response_2.pdf",width=8,height=8)
print (p2)
dev.off()
tiff("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure11C_response_2-72.tif",width=8,height=8,units="in",res=72)
print (p2)
dev.off()
tiff("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure11C_response_2-300.tif",width=8,height=8,units="in",res=300)
print (p2)
dev.off()


gg<-aa[,c(9,11,8)]
gg<-gg[which(gg[,2] %in% c("CR/PR","SD/PD")),]
library("pROC")
data<-data.frame(status=gg$binaryResponse,risk=as.numeric(gg$RBPscore.class_score=="high"),score=gg$RBPscore.score)
#data$status<-ifelse(data$status=="CR/PR",1,0)
a<-plot.roc(status~risk,data,col="#66C2A5",xlim=c(1,0))
a$auc
pdf("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure11D_roc.pdf", width = 8, height = 8)
plot.roc(status~risk,data,col="#66C2A5",xlim=c(1,0))
legend("bottomright", legend=c("AUC: 0.53"),col=c("#66C2A5"), lwd=2)
dev.off()
tiff("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure11D_roc-72.tif", width = 8, height = 8, units = "in", res = 72)
plot.roc(status~risk,data,col="#66C2A5",xlim=c(1,0))
legend("bottomright", legend=c("AUC: 0.53"),col=c("#66C2A5"), lwd=2)
dev.off()
tiff("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Figure/Figure11D_roc-300.tif", width = 8, height = 8, units = "in", res = 300)
plot.roc(status~risk,data,col="#66C2A5",xlim=c(1,0))
legend("bottomright", legend=c("AUC: 0.53"),col=c("#66C2A5"), lwd=2)

dev.off()

save(list=ls(),file="/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/04.output/data20200811.RData")
































