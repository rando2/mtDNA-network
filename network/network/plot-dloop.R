library("pegas")
library("ape")
library("seqinr")
library("ggplot2")

library("adegenet")

#help from https://stackoverflow.com/questions/31220586/r-how-to-plot-correct-pie-charts-in-haplonet-haplotyp-networks-pegas-ape-a
#and https://johnbhorne.wordpress.com/2016/09/15/still-making-haplotype-networks-the-old-way-how-to-do-it-in-r/

#THIS IS SET UP FOR DLOOP

#import data from csv and mega fasta alignment
lithaplo<-read.csv("./dloop-results.csv", stringsAsFactors=FALSE)
lithaplo[is.na(lithaplo)] <- 0
names(lithaplo)[1:3] <- c("Haplo","Acc","Total")
aln <- read.dna("./dloop-results.fas", format="fasta", as.matrix=TRUE)

#Change syntax so they all match 
#for future: gsub with RE to detect ws or _
#newlabs <- gsub(pattern ="isolate_",replacement ="haplotype_",x=  labels(aln))
#newlabs <- gsub(pattern ="isolate ",replacement ="haplotype_",x=  newlabs)
newlabs <- gsub(pattern ="haplotype ",replacement ="haplotype_",x=  labels(aln))
newlabs <- gsub(pattern =" D-loop",replacement ="_D-loop",x=  newlabs)

#Split before and extract only second element
newlabs<- unlist(strsplit(x = newlabs, split = "haplotype_"))[c(TRUE, FALSE)]

#Split after and extract only first element
newlabs<- unlist(strsplit(x = newlabs, split = "_D-loop"))[c(FALSE, TRUE)]
newlabs <- c(279, newlabs)
labelMatch <- data.frame(full=trimws(labels(aln)), haplo=newlabs)

#Organize the frequency table to match the order of sequences
popFreq <- lithaplo[,c(1,4:ncol(lithaplo))]
totFreq <- lithaplo[,c(1,3)]
popFreq <- popFreq[match(labelMatch$haplo, popFreq$Haplo),]
totFreq <- totFreq[match(labelMatch$haplo, totFreq$Haplo),]
rownames(popFreq)<-seq(1,nrow(popFreq)) #this is probably just superstition, but I don't like the labels being out of order

#Set up Network
alnHap <- haplotype(aln)
alnNet <- haploNet(alnHap)

#Set up frequency table using code from 
ind.hap<-with( #This doesn't work for my data, but somehow it is valid for pie, so let's doctor it and see what happens
  stack(setNames(attr(alnHap, "index"), rownames(alnHap))),
  table(hap=ind, individuals=rownames(aln)[values])
)

#Copy frequencies from popFreq
if(ncol(ind.hap) < ncol(popFreq)-1){
  while(ncol(ind.hap) < ncol(popFreq)-1){
    ind.hap <- cbind(ind.hap, rep(0,6))
  }
} else if(ncol(ind.hap) > ncol(popFreq)-1){
  ind.hap <- ind.hap[,1:(ncol(popFreq)-1)]
}
colnames(ind.hap) <- colnames(popFreq)[2:ncol(popFreq)]
for(i in seq(1,nrow(ind.hap))){
  for(j in seq(1,ncol(ind.hap))){
    ind.hap[i,j]=popFreq[i,j+1]
  }
}

#Modify attributes of network to fit our actual counts
attr(alnNet,"freq") <- totFreq$Total
attr(alnNet,"labels") <- popFreq$Haplo
#attr(alnNet,"labels")[1] <-279

#Plot the network
pdf("dloop-farm.results.labeled.fixed.pdf")
plot(alnNet, size=attr(alnNet, "freq"), scale.ratio = 20, cex = 0.8, pie=ind.hap, 
     threshold=0)
legend(-250,250, colnames(ind.hap), col=rainbow(ncol(ind.hap)), pch=20)
dev.off()
jpeg("dloop-farm.results.unlabeled.jpg",width = 8, height = 5, units = 'in', 
     res = 300)
plot(alnNet, size=attr(alnNet, "freq"), labels=FALSE,scale.ratio = 20, 
     cex = 0.8, pie=ind.hap, threshold=0)
legend(-250,250, colnames(ind.hap), col=rainbow(ncol(ind.hap)), pch=20)
dev.off()
