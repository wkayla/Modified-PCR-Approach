library(stringr)
library(dplyr)
library("vegan")
library("reshape2")
library("ggplot2")
library("scales")
library(exactRankTests)
library(diagonals)
library(Matrix)

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

name_split<-function(names){lapply(strsplit(names,"/"), tail, n = 1L)}



###Beta Diversity Functions
MH_values=function(data,n){
  
#Creates a block diagonal matrix
rblock<- function(nb) {
  .bdiag(replicate(nb, {
    Matrix(c(0,1,1,0), 2,2) }))
}
  
values=1-as.matrix(vegdist(data,method="morisita"))

a<-rblock(n)
a=as.matrix(a)

#overlay diagonal matrix onto distance matrix
k=cbind.data.frame(values*a)
f=k[k!=0]

#Delete duplicate MH values; only need one per pair
mh=as.data.frame(unique(f))
return(mh)}



###PcOA Plot
PCOA_plot<-function(data,root,color,group=NULL,group_name,paired){
dist=vegdist(data,method="morisita")
pcoa=cmdscale(dist)

otu_rel=data/root
names=as.character(colnames(data))
name_list=name_split(names)


ds=cbind.data.frame(otu_rel,pcoa)
colnames(ds)[(dim(ds)[2]-1):dim(ds)[2]]=c("m1","m2")

# sums=colSums(otu_rel)
# sums=sort(sums,decreasing = T)
# topRA=name_split(names(sums)[1:5])
# 
# vec_sp<-envfit(pcoa, otu_rel)
# vec_sp_df<-as.data.frame(cbind(vec_sp$vectors$arrows*sqrt(vec_sp$vectors$r),vec_sp$vectors$pvals))
# vec_sp_df$species<-name_split(rownames(vec_sp_df))
# colnames(vec_sp_df)<-c("Dim1","Dim2","Pvals","species")
# vec_sp_df<-vec_sp_df[vec_sp_df$species%in%as.character(topRA),]

if(paired==T){
ggplot(ds,aes(x=ds$m1,y=ds$m2,color=factor(color),group=group))+
  geom_point(size=3)+
  geom_line(col="black",linetype="dotted")+
  guides(  size=F, color = guide_legend(group_name))+
  labs(x="Component 1", y="Component 2")+
  theme(text = element_text(size=20))+ 
  # geom_segment(data=vec_sp_df,aes(x=0,xend=Dim1,y=0,yend=Dim2),
  #              arrow = arrow(length = unit(0.5, "cm")),
  #              colour="black",
  #              stat="identity",inherit.aes = FALSE) +
  # geom_text(data=vec_sp_df,aes(x=Dim1,y=Dim2),label=unlist(vec_sp_df$species),inherit.aes = F,
  #           size=5)+
  coord_fixed()+
  xlim(-1,1)+
  ylim(-1,1)}

else{
  ggplot(ds,aes(x=ds$m1,y=ds$m2,color=factor(color)))+
    geom_point(size=3)+
    guides(  size=F, color = guide_legend(group_name))+
    labs(x="Component 1", y="Component 2")+
    theme(text = element_text(size=20))+ 
    geom_segment(data=vec_sp_df,aes(x=0,xend=Dim1,y=0,yend=Dim2),
                 arrow = arrow(length = unit(0.5, "cm")),
                 colour="black",
                 stat="identity",inherit.aes = FALSE) +
    geom_text(data=vec_sp_df,aes(x=Dim1,y=Dim2),label=unlist(vec_sp_df$species),inherit.aes = F,
              size=5)+
    coord_fixed()+
    xlim(-1,1)+
    ylim(-1,1)}}



####Paired boxplot
pair_boxplot<-function(data,x,y,group,xlab,ylab,main){
ggplot(data,aes(x=x,y=y))+
  geom_boxplot()+
  geom_point()+
  geom_line(aes(group=group))+
  guides(col=guide_legend("",override.aes = list(size=3)),size=F)+
  xlab(xlab)+
    ggtitle(main)+
  ylab(ylab)+theme(text = element_text(size=30))}


####Prevalence
prevalence<-function(data,group=NULL){
  if(is.null(group)==F){apply(data,2,function(x){tapply(x,group,function(y){sum(ifelse(y==0,0,1))/(dim(data)[1]/2)})})}
  else{apply(data,2,function(x){sum(ifelse(x==0,0,1))/dim(data)[1]})}}


####Wilcoxon Exact tests table

wilcox_table<-function(data,root,RA_cutoff,Prevalence_cutoff,group,paired){
rel_otu=data/root

Taxa_keep=names(na.omit(apply(rel_otu,2,function(x){ifelse(median(x)>RA_cutoff,x,NA)})))

prevalent=prevalence(data)
Prevalent_keep=names(na.omit(sapply(prevalent,function(x){ifelse(max(x)>Prevalence_cutoff,x,NA)})))

names_keep=Taxa_keep[Taxa_keep%in%Prevalent_keep]


g1_rows=grep(levels(group)[1],group)
g2_rows=grep(levels(group)[2],group)

wilx_rel_otu=rel_otu[,colnames(rel_otu)%in%Taxa_keep]
#Calculate Wilcoxon Exact Signed Rank Tests for Each taxa 
wilx=apply(wilx_rel_otu,2,function(x) wilcox.exact(x[g1_rows],x[g2_rows],
                                                   paired=paired,conf.int = T))
significant_wilx=unlist(lapply(wilx, function(x) {if(x$p.value<.05){return(x$p.value)}}))

keep_wilx=wilx[names(wilx)%in%names(significant_wilx)]
#round p-values to 2 digits
wilx.estimate=round(unlist(lapply(keep_wilx, function(x) {return(x$estimate)})),4)*100
table.p=round(p.adjust(unlist(lapply(keep_wilx, function(x){if(x$p.value){return(x$p.value)}})),method="BH"),2)
table.p=ifelse(table.p<=0.00,0.01,table.p)
table.p=ifelse(table.p==0.01,paste("<",table.p,sep=""),table.p)

#Table One
table=cbind.data.frame(names(keep_wilx),
                       wilx.estimate,table.p)

# table.ra.op=table.ra.op[-13,]
#Name columns of Table 1
colnames(table)=c("Genera",
                  "Estimated Change in RA %",
                  "P-Value")

row.names(table)<-seq(nrow(table))
return(table)}


#####Stacked Bar Chart
stacked_barcahrt<-function(data,y,x,fill){
  ggplot(data=data, 
       aes(y=y, 
           x=x, 
           fill=factor(fill))) + 
  geom_bar(stat="identity",width=1)+
  scale_x_discrete()+
  coord_flip()+
  facet_grid(factor(new_collection)~Subject)+
  theme(legend.position = "bottom", 
        legend.box = "horizontal",
        axis.text.x = element_text(angle = 90, hjust = 1),
        legend.text = element_text(size=10))+
  labs(x="Collection", 
       y = "Percentage of Total")+
  labs(fill="")}


