rel_otu=otu[,grep("Bacteria",colnames(otu))]/otu[,"root"]
rel_otu=cbind.data.frame(Metadata$Subject,otu$Id,rel_otu)



triplicates=Metadata[Metadata$test_spike==1,]
triplicates$Triplicate=rep(1:3,20)


Modified_Triplicates=triplicates[triplicates$modified=="Modified",]
Standard_Triplicates=triplicates[triplicates$modified=="Standard",]

otu_modified_triplicates=rel_otu[rel_otu$'otu$Id'%in%Modified_Triplicates$Lib,]
otu_standard_triplicates=rel_otu[rel_otu$'otu$Id'%in%Standard_Triplicates$Lib,]

#Modified

modified_subject_data=split(otu_modified_triplicates, factor(otu_modified_triplicates$'Metadata$Subject'), drop = FALSE)

sdv_modified=lapply(modified_subject_data,function(x)
  apply(x,2,sd))
mean_modified=lapply(modified_subject_data,function(x)
  apply(apply(apply(x,2,as.character),2,as.numeric),2,mean))
cv_modified=na.omit(unlist(sdv_modified))/na.omit(unlist(mean_modified))


#Standard

standard_subject_data=split(otu_standard_triplicates, factor(otu_standard_triplicates$'Metadata$Subject'), drop = FALSE)

sdv_standard=lapply(standard_subject_data,function(x)
  apply(x,2,sd))
mean_standard=lapply(standard_subject_data,function(x)
  apply(apply(apply(x,2,as.character),2,as.numeric),2,mean))
cv_standard=na.omit(unlist(sdv_standard))/na.omit(unlist(mean_standard))


t.test(cv_modified,cv_standard,paired = T)


#Standard
max_standard=cbind.data.frame(cv_standard,na.omit(unlist(mean_standard)))
max20_standard=subset(max_standard,(!(max_standard$cv_standard)=="NaN") & (max_standard$cv_standard<0.2)& (max_standard$cv_standard>0.195))


#Modified
max_modified=cbind.data.frame(cv_modified,na.omit(unlist(mean_modified)))
max20_modified=subset(max_modified,
                      (!(max_modified$cv_modified)=="NaN") & 
                        (max_modified$cv_modified<0.2)& 
                        (max_modified$cv_modified>0.195))

##################### Coefficient of Variation Plot ###########################
cv_standard_data=cbind.data.frame(na.omit(unlist(mean_standard)),cv_standard)
cv_standard_data=na.omit(cv_standard_data)

cv_modified_data=cbind.data.frame(na.omit(unlist(mean_modified)),cv_modified)
cv_modified_data=na.omit(cv_modified_data)

plot(jitter(cv_standard_data$`na.omit(unlist(mean_standard))`, amount =.005),jitter(cv_standard_data$cv_standard, amount =.005),
     pch=16,
     ylab = "Coefficient of Variation", 
     xlab="Average Relative Abundance",col="Blue",cex.lab=1.65,las=2, xlim=c(0,0.015))
points(jitter(cv_modified_data$`na.omit(unlist(mean_modified))`,amount=.005),jitter(cv_modified_data$cv_modified,amount=.005),
       col="Black")
abline(h=.2,col="red",lty=2)
abline(v=.01,col="red",lty=2)
# axis(side = 1, at = .01, las=2,cex=.5, padj=1)
axis(side = 2, at = .2,las=2)
