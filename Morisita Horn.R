Mh_values=1-as.matrix(vegdist(paired_otu[,4:438],method="morisita"))


#Creates a block diagonal matrix
rblock<- function(nb) {
  .bdiag(replicate(nb, {
    Matrix(c(0,1,1,0), 2,2) }))
}

a<-rblock(52)
a=as.matrix(a)

#overlay diagonal matrix onto distance matrix
k=cbind.data.frame(Mh_values*a)
f=k[k!=0]

#Delete duplicate MH values; only need one per pair
mh=as.data.frame(unique(f))

boxplot(mh,ylab="Morisita Horn Pairwise Similarity", ylim=c(.994,1),outline=F)
stripchart(mh[,1], vertical = TRUE, data = mh, 
           method = "jitter", add = TRUE, pch = 20, col = 'blue')


summary(mh)