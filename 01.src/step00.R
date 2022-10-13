rm(list = ls())
options(stringsAsFactors = F)
setwd("/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Runtime/01.data/")
##
RBP_gene1=read.csv("./01.LIHC_88RBP_gene.csv",header = T)
RBP_genelist1=unique(RBP_gene1[,1])
##
RBP_gene2=read.csv("02.A census of human RNA-binding proteins-Gerstbergers-1542.csv",header = T)
RBP_genelist2=unique(RBP_gene2[,1])
###
RBP_gene3=read.csv("03.SONAR Discovers RNA-Binding Proteins from Analysis of Large-Scale Protein-Protein Interactomes-1916.csv",header = T)
RBP_gene3=subset(RBP_gene3,RBP_gene3$score_SONAR>0.79,)
RBP_genelist3=unique(RBP_gene3[,1])
###
RBP_gene4=read.csv("04.Insights into RNA biology from an Atlas of Mammalian mRNA bing proteins.csv",header = T)
RBP_genelist4=unique(RBP_gene4$Symbol)
##
RBP_gene5=read.csv("05.The mRNA-bound proteome and its global occupancy profile on protein-coding transcripts.csv",header = T)
RBP_gene5=subset(RBP_gene5,RBP_gene5$Functional.Group %in% c("RNAbinding"))
RBP_genelist5=unique(RBP_gene5$Offical.gene.symbol)
###
RBP_gene6=read.csv("06.The RNA-binding proteomes from yeast to man harbour conserved enigmRBPs.csv",header = T)
RBP_gene6_1=subset(RBP_gene6,RBP_gene6$Conserved..core..RBP %in% c("core"))
RBP_gene6=RBP_gene6_1
RBP_genelist6=unique(RBP_gene6$HuH7.mRNA.interactome_ENSEMBL.Gene.Name)
###
RBP_gene7=read.csv("07.Comprehensive Identification of RNA-Binding Domains in Human Cells.csv",header = T)
RBP_gene7=subset(RBP_gene7,!RBP_gene7$domain %in% c("other"))
RBP_genelist7=unique(RBP_gene7$Symbol)
###
RBP_gene8=read.csv("08.Serial interactome capture of the human cell nucleus.csv",header = T)
RBP_genelist8=unique(RBP_gene8$symbol)
###
RBP_gene9=read.csv("09.The Human RNA-Binding Proteome and Its Dynamics during Translational -2483.csv",header = T)
RBP_genelist9=unique(RBP_gene9$Gene.name)
###
RBP_gene10=read.csv("10.Transcriptome-wide discovery of coding and noncoding RNA-binding proteinsxlsx-597.csv",header = T)
RBP_genelist10=unique(RBP_gene10$Gene)
###
RBP_genelist=unique(c(RBP_genelist1,RBP_genelist2,RBP_genelist3,RBP_genelist4,RBP_genelist5,RBP_genelist6,RBP_genelist7,RBP_genelist8,RBP_genelist9,RBP_genelist10))
write.table(RBP_genelist,file="../04.output/01.数据搜集和预处理/01.RBP_genelist.txt",sep="\t",quote = F,row.names = F)

















