

Paired=Metadata[Metadata$Paired==1,]

paired_otu=otu[otu$Id%in%Paired$Lib,]

div=diversity(paired_otu[,4:438],index="shannon")
pair=rep(1:52,each=2)
div=cbind.data.frame(div,Paired[,6:8],pair)
colnames(div)=c("diversity",colnames(Paired)[6:8],"pairs")

ggplot(div,aes(x=modified,y=diversity)) + 
  geom_boxplot(aes(modified, diversity))+
  geom_point() +
  geom_line(aes(group=pairs,col=Subject))+
  xlab("")+
  ylab("Shannon-H Diversity Index")+
  theme(text = element_text(size=20))


#wilcoxon exact by subject between pairs
subject_pair_data=split(div, 
                        factor(div$Subject), 
                        drop = FALSE)

wilcox_paired_subject=sort(unlist(lapply(subject_pair_data,function(x){
  wilcox.exact(x$diversity~x$modified,paired=T)$p.value})),decreasing=T)
wilcox_paired=wilcox.exact(div$diversity~div$modified,paired=T)

tapply(div$diversity, div$modified,summary)

summary(unlist(lapply(subject_pair_data,dim)))
