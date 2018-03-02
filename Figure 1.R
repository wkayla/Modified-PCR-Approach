############################ Bacterial Load
library("sas7bdat")
cv_data=read.sas7bdat("C:/Users/wkayla/OneDrive for Business/VAP/Methods/Data/cv_data.sas7bdat")
load_data=read.sas7bdat("C:/Users/wkayla/OneDrive for Business/VAP/Methods/Data/qpcr_2017_02_15.sas7bdat")
load_data=load_data[grepl("T",load_data$Name),]
IDs_load=read.csv("C:/Users/wkayla/OneDrive for Business/VAP/Methods/Data/IDS_load.csv")

######################################### September 28, 2017 Updates

Metadata$Name=strtrim(Metadata$Lib, 5)
libraries=split(Metadata,Metadata$Name)
sequenced=unlist(lapply(libraries,function(x){dim(x)[1]}))

#If only 1 sequence assign to 1, if 2 or more sequences assign to 2
sequenced=ifelse(sequenced==0,0,ifelse(sequenced==1,1,2))    
sequenced_data=cbind.data.frame(names(libraries),sequenced)
colnames(sequenced_data)=c("libs","sequenced")

#Reattach libraries
Meta_seq=merge(Metadata,sequenced_data,by.x="Name",by.y="libs")

#Parse down to only one sequence per sample given multiple amplifications
Meta_seq=Meta_seq[match(unique(Meta_seq$Name),Meta_seq$Name),]

#Merge with load data
seq_load=merge(Meta_seq,load_data,by.x="Name",by.y="Name")


#Assign labels to number of sequences
seq_load$sequenced=factor(seq_load$sequenced,levels=c(1,2),labels=c(" Modified","Both"))

#Find samples with no sequence data
no_sequence=load_data[!(load_data$Name%in%Metadata$Lib),]
no_sequence=merge(no_sequence,IDs_load,by.x="Name",by.y = "Molecular.Id")
no_sequence$sequenced = rep("Neither",dim(no_sequence)[1])

new_seq=seq_load[,c("Subject","Name","sequenced","lq_all")]
no_sequence=no_sequence[,c("Subject","Name","sequenced","lq_all")]

new_seq=data.frame(apply(new_seq,2,as.character))
no_sequence=data.frame(apply(no_sequence,2,as.character))
new_seq_data=rbind.data.frame(new_seq,no_sequence)
new_seq_data$Subject=gsub("s","",gsub("S","",new_seq_data$Subject))

# new_seq_data$Subject=mapvalues(new_seq_data$Subject,from=levels(factor(new_seq_data$Subject)),to=c("01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16"))
# new_seq_data$Subject=factor(new_seq_data$Subject)
ggplot(new_seq_data)+
  geom_point(aes(x=factor(new_seq_data$Subject),y=as.numeric(as.character(new_seq_data$lq_all)),col=relevel(new_seq_data$sequenced,"Neither")))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  guides(col=guide_legend(title="Sequenced"))+
  ylab("Bacterial Load (log10 16S rRNA copy number/reaction)”")+
  xlab("Subject")

tapply(as.numeric(as.character(new_seq_data$lq_all)),new_seq_data$sequenced,summary)

cor.test(as.numeric(new_seq_data$lq_all),as.numeric(new_seq_data$sequenced), method = "spearman")

rel_otu=otu[,grep("Bacteria",colnames(otu))]/otu$root

prev=prevalence(rel_otu)
prev[names(prev)=="Bacteria/Proteobacteria/Gammaproteobacteria/Enterobacteriales/Enterobacteriaceae"  ]
tapply(rel_otu$`Bacteria/Proteobacteria/Gammaproteobacteria/Enterobacteriales/Enterobacteriaceae`,Metadata$modified,summary)

rel_otu$sub=Metadata$Subject
rel_otu$collection=Metadata$Collection
rel_otu$modified=Metadata$modified
table(rel_otu[rel_otu$sub=="S1027",c("Bacteria/Proteobacteria/Gammaproteobacteria/Enterobacteriales/Enterobacteriaceae","collection","modified")])


hist(rank(new_seq_data$lq_all))

new_seq_data$ranked_load=sort(rank(new_seq_data$lq_all),decreasing = T)
mod=lm(rank(new_seq_data$lq_all)~relevel(new_seq_data$sequenced,"Neither"),data = new_seq_data)
sum=summary(mod)
exp(coef(sum))
sum
kruskal.test(new_seq_data$lq_all~new_seq_data$sequenced)

no_neither=new_seq_data[new_seq_data$sequenced!="Neither",]
wilcox.exact(as.numeric(as.character(no_neither$lq_all))~no_neither$sequenced)


no_both=new_seq_data[new_seq_data$sequenced!="Both",]
wilcox.exact(as.numeric(as.character(no_both$lq_all))~no_both$sequenced,paired=F)

tapply(as.numeric(as.character(new_seq_data$lq_all)),new_seq_data$sequenced,median)
