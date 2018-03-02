############################### AIM 1 ################################################



#Triplicates
triplicates=Metadata[Metadata$test_spike==1,]

otu_partial=otu[otu$Id%in%Metadata$Lib,]
rel_otu_partial=otu_partial[,4:438]/otu_partial$root
############################# find Most Abundant #############################


#Alternative
sums=colSums(rel_otu_partial)
sums=sort(sums,decreasing = T)
topRA=names(sums)[1:10]

abundant=cbind.data.frame(rel_otu_partial[,as.character(topRA)],otu_partial$Id)
abundant_trip=abundant[abundant$`otu_partial$Id`%in%triplicates$Lib,]
abundant_Paired=abundant[abundant$`otu_partial$Id`%in%paired_meta$Lib,]

#################Sort Metadata by Triplicate
Tripnum=rep(1:3,20)
Tripnum_data=c(Tripnum,rep(0,108))
Meta_trip_sort=Metadata[order(Metadata$test_spike,decreasing=T),]
Meta_trip_sort=cbind.data.frame(Meta_trip_sort,Tripnum_data)
###################################


otu_trip=cbind(triplicates$Subject,Tripnum,triplicates$Collection ,triplicates$modified,abundant_trip[,1:10])

meta=Meta_trip_sort[,c(1,2,6,8,9)]
colnames(abundant)=c(colnames(abundant)[1:10],"Lib")
otu_meta=merge(meta,abundant,"Lib")
otu_meta=otu_meta[,-1]

names=name.split(colnames(abundant))

colnames(otu_trip)=c("Subject","Triplicate","Collection","Modified",names[1:10])

colnames(otu_meta)=c("Collection","Subject","Modified","Triplicate",names[1:10])

xlong_trip<-melt(otu_trip,id.vars=c("Subject","Triplicate","Modified","Collection"))
xlong_trip_sub=xlong_trip[xlong_trip$Subject=="s1003",]


######## relabel collection
otu_meta_graph=otu_meta[otu_meta$Triplicate==0,]

list_subject=lapply(split(otu_meta_graph,factor(otu_meta_graph$Subject)),function(x){x$Collection})

new_list_subject =unlist(lapply(list_subject,function(x){new=match(x,unique(x))}))

otu_meta_graph$new_collection=new_list_subject



xlong<-melt(otu_meta_graph[,2:15],id.vars=c("Subject","Triplicate","new_collection","Modified"))

subject_names=seq(1,14,1)
subject_labeller <- function(variable,value){
  return(subject_names[value])
}

#FIGURE IS CORRECT!
ggplot(data=xlong_trip,
       aes(y=value,
           x=as.factor(Modified),
           fill=factor(variable))) +
  geom_bar(stat="identity",width=1)+
  scale_x_discrete()+
  coord_flip()+
  facet_grid(factor(Triplicate)~Subject,labeller=subject_labeller)+
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        axis.text.x = element_text(angle = 90, hjust = 1),
        legend.text = element_text(size=10))+
  labs(x="Triplicate",
       y = "Percentage of Total")+
  labs(fill="")


xlong$Subject=gsub("s","",gsub("S","",xlong$Subject))
#This figure is correct!
ggplot(data=xlong, 
       aes(y=value, 
           x=as.factor(Modified), 
           fill=factor(variable))) + 
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
  labs(fill="")


