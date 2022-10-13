setwd('/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Table/')
#下载Copy Number Variation数据
rm(list = ls())
options(stringsAsFactors = F)
set.seed(123456)
library(tidyverse)
library(DT)
TCGAKIRC_segment_cluster_high = read.table('/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Table/SupplementaryTable8_high_seg-v1.txt',header = T)
TCGAKIRC_segment_cluster_low  = read.table('/Users/liuyang/Desktop/资料/项目报告/GAP-F200610001-基于RBP/Reuslts/Table/SupplementaryTable8_low_seg-v1.txt' ,header = T) 
#
delete_sample_segment<-function(TCGAKIRC_segment) {
  delete_sample<-c()
  length=nrow(TCGAKIRC_segment)
  if (length>1){
    for ( i in 1:(length-1)){
      if (TCGAKIRC_segment[(i+1),3]<=TCGAKIRC_segment[i,3] & TCGAKIRC_segment[i+1,4]>=TCGAKIRC_segment[i,4]){
        delete_sample<-append(delete_sample,i)
      }
      if ( TCGAKIRC_segment[(i+1),3]>=TCGAKIRC_segment[i,3] & TCGAKIRC_segment[i+1,4]<=TCGAKIRC_segment[i,4]){
        delete_sample<-append(delete_sample,i)
        TCGAKIRC_segment[i+1,3]=TCGAKIRC_segment[i,3]
        TCGAKIRC_segment[i+1,4]=TCGAKIRC_segment[i,4]
        TCGAKIRC_segment[i+1,5]=(TCGAKIRC_segment[i,5]+TCGAKIRC_segment[i+1,5])/2
      }
      if (TCGAKIRC_segment[(i+1),3]>=TCGAKIRC_segment[i,3] & TCGAKIRC_segment[(i+1),3]<TCGAKIRC_segment[i,4] &TCGAKIRC_segment[i+1,4]>=TCGAKIRC_segment[i,4]){
        delete_sample<-  append(delete_sample,i)
        TCGAKIRC_segment[i+1,3]=TCGAKIRC_segment[i,3]
        TCGAKIRC_segment[i+1,5]=(TCGAKIRC_segment[i,5]+TCGAKIRC_segment[i+1,5])/2
      }
      if (TCGAKIRC_segment[(i+1),3]<TCGAKIRC_segment[i,3] & TCGAKIRC_segment[i+1,4]>TCGAKIRC_segment[i,3] & TCGAKIRC_segment[i,4]>TCGAKIRC_segment[i+1,4]){
        delete_sample<-  append(delete_sample,i)
        TCGAKIRC_segment[i+1,4]=TCGAKIRC_segment[i,4]
        TCGAKIRC_segment[i+1,5]=(TCGAKIRC_segment[i,5]+TCGAKIRC_segment[i+1,5])/2
      }
    }
  }
  if (length(delete_sample)>0){
    TCGAKIRC_segment=TCGAKIRC_segment[-delete_sample,]
  }
  return (TCGAKIRC_segment)
}

####
sps<-split(TCGAKIRC_segment_cluster_high,TCGAKIRC_segment_cluster_high[,c('sample','Chrom')],drop=T)
data_matrix=matrix(1:5,ncol=5,nrow=1,byrow = T)
colnames(data_matrix)=c('sample','Chrom','Start','End','value')
for (i in 1:length(sps)){
  sps[[i]]=sps[[i]] %>% arrange(Start)
  sps[[i]]=delete_sample_segment(sps[[i]])
  data_matrix=rbind(data_matrix,sps[[i]])
}
data_matrix=data_matrix[-1,]
TCGAKIRC_segment_cluster_high<-data_matrix
####
sps<-split(TCGAKIRC_segment_cluster_low,TCGAKIRC_segment_cluster_low[,c('sample','Chrom')],drop=T)
data_matrix=matrix(1:5,ncol=5,nrow=1,byrow = T)
colnames(data_matrix)=c('sample','Chrom','Start','End','value')
for (i in 1:length(sps)){
   sps[[i]]=sps[[i]] %>% arrange(Start)
   sps[[i]]=delete_sample_segment(sps[[i]])
   data_matrix=rbind(data_matrix,sps[[i]])
 }
data_matrix=data_matrix[-1,]
TCGAKIRC_segment_cluster_low<-data_matrix
write.table(TCGAKIRC_segment_cluster_high,file = "./TCGALIHC-cluster_high-segment.seg",sep = '\t',quote = F,row.names = F)
write.table(TCGAKIRC_segment_cluster_low,file = "./TCGALIHC-cluster_low-segment.seg",sep = '\t',quote = F,row.names = F)

#setwd('./')

