
otu_plot=paired_otu[,grep("Bacteria",colnames(otu))]/paired_otu[,"root"]
# BA_otu_plot=otu_plot[as.character(colnames(abundant))]
otu_meta=cbind(paired_meta$Subject,paired_meta$Collection, paired_meta$modified,otu_plot)
colnames(otu_meta)=c("Subject","Collection","Modified",colnames(otu_plot))

xlong<-melt(otu_meta,id.vars=c("Subject","Collection","Modified"))

Modified_long=xlong[xlong$Modified=="Modified",]
Standard_long=xlong[xlong$Modified=="Standard",]

# bland.altman.plot(Modified_long$value,Standard_long$value,  two=.01)
save=bland.altman.stats(as.numeric(Modified_long$value),as.numeric(Standard_long$value),  two=.01)
apply(save$groups,2,summary)
range=max(save$diffs)-min(save$diffs)
differences=sort(save$diffs,decreasing=T)

BAplot=cbind.data.frame(save$diffs,save$means)
ggplot(BAplot,aes(BAplot[,2],BAplot[,1]))+geom_point()+geom_hline(yintercept=.01,lty=2)

BAplot%>%
  ggplot( aes(`save$means`,`save$diffs`))+
  geom_point()+
  # geom_text(label=name.split(as.character(Modified_long$variable)))+
  ylab("Differences (Modified-Standard)")+xlab("Average Relative Abundance")+
  geom_hline(yintercept=c(0.02,-0.02),lty=2)+
  geom_hline(yintercept=0,lty=1)+
  ylim(-.1,.1)


