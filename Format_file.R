

otu_from_explicet=function(data){
  colnames=data[,1]
  Id=colnames(data)[-1]
  data=t(data)
  colnames(data)=colnames
  data=data[-1,]
  data=apply(apply(data,2,as.character),2,as.numeric)
  data=cbind.data.frame(Id,data)
  rownames(data)=seq(nrow(data))
  return(data)}

name.split<-function(names){
  names=names
  save=strsplit(names,"/")
  
  h=0
  for(i in 1:length(names)){
    h[i]=length(save[[i]])}
  
  i=0
  name.list=NULL
  for(i in 1:length(save)){
    name.list[i]=save[[i]][h[i]]}
  return(name.list)
}


#Creates a block diagonal matrix
rblock<- function(nb) {
  .bdiag(replicate(nb, {
    Matrix(c(0,1,1,0), 2,2) }))
}



library("vegan")
library("reshape2")
library("ggplot2")
library("scales")
library(diagonals)
library(Matrix)
library(sas7bdat)
library(exactRankTests)
library(dplyr)

#Read in files
Metadata = read.table("C:/Users/wkayla/OneDrive for Business/VAP/Methods/Data/VAP_Methods_Corrected2_Dec2015_Metadata.txt",header=T)

data= read.table("C:/Users/wkayla/OneDrive for Business/VAP/Methods/Data/VAP_Methods_Corrected2_Dec2015_Workspace_1_OTU.txt",header=T)
otu=otu_from_explicet(data)

#Create Second subject for S1004 (two samples)
# Metadata$Subject=replace(as.character(Metadata$Subject),grepl("T1639",Metadata$Lib)==T,"S10042")

#reformat otu table
otu=otu_from_explicet(data)
colnames(otu)=gsub(">","",colnames(otu))

otu$Pseudomonaplus=otu$"Bacteria/Proteobacteria/Gammaproteobacteria/Pseudomonadales"+otu$"Bacteria/Proteobacteria/Gammaproteobacteria/Pseudomonadales/Pseudomonadaceae/Pseudomonas"
otu$"Bacteria/Proteobacteria/Gammaproteobacteria/Pseudomonadales"=NULL
otu$"Bacteria/Proteobacteria/Gammaproteobacteria/Pseudomonadales/Pseudomonadaceae/Pseudomonas"=NULL

#Create a variable to denote which samples were modified
modified=ifelse(grepl("X",Metadata$Lib)==T,0,1)
Metadata$modified=factor(modified,levels=c(0,1),labels=c("Standard","Modified"))

