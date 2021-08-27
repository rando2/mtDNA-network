#install.packages("C:/Users/rando2/Dropbox/Work Back-Up/Projects/yms/mtDNA/Paper Write-Up/Haplotypes/network/pegas_0.12-1.tar.gz", repos = NULL, type = "source")
#library("pegas")
library("ape")
library("seqinr")
library("ggplot2")

library("adegenet")
setwd("C:/Users/rando2/Dropbox/Work Back-Up/Projects/yms/mtDNA/Paper Write-Up/Haplotypes/network")
install.packages("./pegas_0.12-1.tar.gz", repos = NULL, type="source")
library("pegas")

#THIS ONE PLOTS CYTOCHROME B
#help from https://stackoverflow.com/questions/31220586/r-how-to-plot-correct-pie-charts-in-haplonet-haplotyp-networks-pegas-ape-a
#and https://johnbhorne.wordpress.com/2016/09/15/still-making-haplotype-networks-the-old-way-how-to-do-it-in-r/

#import data from csv and mega fasta alignment
lithaplo<-read.csv("./cytb-results.csv", stringsAsFactors=FALSE)
lithaplo[is.na(lithaplo)] <- 0
names(lithaplo)[1:3] <- c("Haplo","Acc","Total")
aln <- read.dna("./cytb-farm-results.fas", format="fasta", as.matrix=TRUE)

#Change syntax so they all match 
#for future: gsub with RE to detect ws or _
newlabs <- gsub(pattern ="isolate_",replacement ="haplotype_",x=  labels(aln))
newlabs <- gsub(pattern ="isolate ",replacement ="haplotype_",x=  newlabs)
newlabs <- gsub(pattern ="haplotype ",replacement ="haplotype_",x=  newlabs)
newlabs <- gsub(pattern =" cytochrome",replacement ="_cytochrome",x=  newlabs)

#Split before and extract only second element
newlabs<- unlist(strsplit(x = newlabs, split = "haplotype_"))[c(FALSE, TRUE)]

#Split after and extract only first element
newlabs<- unlist(strsplit(x = newlabs, split = "_cyt"))[c(TRUE, FALSE)]
labelMatch <- data.frame(full=trimws(labels(aln)), haplo=newlabs)

#Organize the frequency table to match the order of sequences
popFreq <- lithaplo[,c(1,4:ncol(lithaplo))]
totFreq <- lithaplo[,c(1,3)]
popFreq <- popFreq[match(labelMatch$haplo, popFreq$Haplo),]
totFreq <- totFreq[match(labelMatch$haplo, totFreq$Haplo),]
rownames(popFreq)<-seq(1,nrow(popFreq)) #this may just be superstition - I don't like the labels being out of order

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

#Plot the network
pdf("cytb-farm.results.labeled.pdf")
#plot(alnNet, pie= ind.hap, scale.ratio=50,  labels=TRUE, legend=c(-100,200),
#     size =attr(alnNet,"freq")) #
plot(alnNet, size=attr(alnNet, "freq"), scale.ratio = 20, cex = 0.8, pie=ind.hap, threshold=0)
legend(-200,200, colnames(ind.hap), col=rainbow(ncol(ind.hap)), pch=20)
dev.off()
pdf("cytb-farm.results.unlabeled.pdf")
#plot(alnNet, pie=ind.hap, scale.ratio=25,  labels=FALSE, legend=c(-100,200),
#     size=attr(alnNet,"freq"), cex=200) #scale.ratio = .05,
plot(alnNet, size=attr(alnNet, "freq"), labels=FALSE,scale.ratio = 20, cex = 0.8,threshold=0,pie=ind.hap)
legend(-200,200, colnames(ind.hap), col=rainbow(ncol(ind.hap)), pch=20)
dev.off()
