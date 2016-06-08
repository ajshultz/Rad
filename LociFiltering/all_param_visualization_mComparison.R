setwd("~/Dropbox/HFRad-Tags/HF_rxstacks_results/PairTestingCombinedResults")

num_loci <- read.csv("Num_Loci.csv")
rownames(num_loci) <- num_loci$sample
num_loci <- num_loci[,2:ncol(num_loci)]
boxplot(num_loci)

num_excess <- read.csv("Num_Excess_Allele_loci.csv")
rownames(num_excess) <- num_excess$sample
num_excess <- num_excess[,2:ncol(num_excess)]
mean_num_excess <- apply(num_excess,2,"mean")
sd_num_excess <- apply(num_excess,2,"sd")
boxplot(num_excess)
boxplot(num_excess,ylim=c(0,100))

plot(mean_num_excess,ylim=c(-20,100))
arrows(x0=c(1:8),y0=mean_num_excess,y1=mean_num_excess-sd_num_excess,code=0)
arrows(x0=c(1:8),y0=mean_num_excess,y1=mean_num_excess+sd_num_excess,code=0)

per_mult <- read.csv("Percentage_Multiple_Pair_loci.csv")
rownames(per_mult) <- per_mult$sample
per_mult <- per_mult[,2:ncol(per_mult)]
mean_per_mult <- apply(per_mult,2,"mean")
sd_per_mult <- apply(per_mult,2,"sd")
boxplot(per_mult,main="Percentage of Loci with Multiple Pairs")
plot(mean_per_mult,ylim=c(0,0.11))
arrows(x0=c(1:8),y0=mean_per_mult,y1=mean_per_mult-sd_per_mult,code=0)
arrows(x0=c(1:8),y0=mean_per_mult,y1=mean_per_mult+sd_per_mult,code=0)
abline(h=0)


per_excess <- read.csv("Percentage_Excess_Allele_loci.csv")
rownames(per_excess) <- per_excess$sample
per_excess <- per_excess[,2:ncol(per_excess)]
mean_per_excess <- apply(per_excess,2,"mean")
sd_per_excess <- apply(per_excess,2,"sd")
boxplot(per_excess,main="Percentage of Loci with More than Two Alleles")
boxplot(per_excess,ylim=c(0,0.008),main="Percentage of Loci with More than Two Alleles")
plot(mean_per_excess,ylim=c(-.001,0.006))
arrows(x0=c(1:8),y0=mean_per_excess,y1=mean_per_excess-sd_per_excess,code=0)
arrows(x0=c(1:8),y0=mean_per_excess,y1=mean_per_excess+sd_per_excess,code=0)
abline(h=0)


layout(matrix(1:2,nrow=2))
boxplot(per_mult,main="Percentage of Loci with Multiple Pairs")
boxplot(per_excess,ylim=c(0,0.008),main="Percentage of Loci with More than Two Alleles")

AllInds <- read.csv("~/Dropbox/HFRad-Tags/HF_rxstacks_results/_mParameter_rxstacks_AllLociComparions.csv",row.names=1)
AllInds <- read.csv("~/Dropbox/HFRad-Tags/HF_rxstacks_results/_mParameter_AllLociComparions.csv",row.names=1)

AllInds <- as.data.frame(t(AllInds))

plot(AllInds$VCFTotalLoci,type="l",ylab="Number of Loci",main="Loci After Filtering",xlab="M Parameter")
points(AllInds$VCFTotalLoci)

plot(AllInds$VCFAvgPi,type="l",ylab="Pi",main="Average Subpopulation Pi value",xlab="M Parameter")
points(AllInds$VCFAvgPi)

plot(AllInds$Blastmorethan1hite-40,col="red",ylim=c(0,900),type="l",xlab="M Parameter",ylab="Number of Loci")
lines(AllInds$VCFTooManyAlleles_m4,add=T)
legend(x="topright",legend=c("More than 1 Blast e-40 hit","Number of loci with more than 2 alleles"),fill=c("red","black"))

plot(AllInds$BlastMorethan1hit,col="red",ylim=c(0,900),type="l",xlab="M Parameter",ylab="Number of Loci")
lines(AllInds$VCFTooManyAlleles_m4,add=T)
legend(x="topright",legend=c("More than 1 Blast e-20 hit","Number of loci with more than 2 alleles"),fill=c("red","black"))

plot(AllInds$Blastmorethan1hite-40,col="red",ylim=c(0,900),type="l",xlab="M Parameter",ylab="Number of Loci")
lines(AllInds$BlastMorethan1hit,add=T)



layout(matrix(1:3,nrow=3))
par(cex.main=1,cex.axis=.75,mar=c(0.5,3.5,1.5,1)+0.1)

labels=c(rownames(AllInds))

plot(AllInds$VCFTotalLoci,type="l",ylab="",main="Number of Loci After Filtering",xlab="",xaxt="n",mgp=c(2.1,1,0))
points(AllInds$VCFTotalLoci)
mtext("Number of Loci",side=2,line=2,cex=.75)
boxplot(per_mult[,c(1,3,5)],main="Percentage of Loci with Multiple Pairs",xlab="",xaxt="n",ylab="")
mtext("Percentage of Loci",side=2,line=2,cex=.75)
par(cex.main=1,cex.axis=.75,mar=c(2.5,3.5,1.5,1)+0.1)
boxplot(per_excess[,c(1,3,5)],ylim=c(0,0.008),main="Percentage of Loci with More than Two Alleles",xlab="",ylab="",las=0)
mtext("Percentage of Loci",side=2,line=2,cex=.75)




layout(matrix(1:4,nrow=4))
par(cex.main=1,cex.axis=.75,mar=c(0.5,3.5,1.5,1)+0.1)

labels=c(rownames(AllInds))

plot(AllInds$VCFTotalLoci,type="l",ylab="",main="Number of Loci After Filtering",xlab="",xaxt="n",mgp=c(2.1,1,0))
points(AllInds$VCFTotalLoci)
mtext("Number of Loci",side=2,line=2,cex=.75)
plot(AllInds$VCFTooManyAlleles_m4,type="l",ylab="",main="Number of Loci with More than Two Alleles",xlab="",xaxt="n",mgp=c(2.1,1,0))
boxplot(per_mult[,c(1,3,5)],main="Percentage of Loci with Multiple Pairs",xlab="",xaxt="n",ylab="")
mtext("Percentage of Loci",side=2,line=2,cex=.75)
par(cex.main=1,cex.axis=.75,mar=c(2.5,3.5,1.5,1)+0.1)
boxplot(per_excess[,c(1,3,5)],ylim=c(0,0.008),main="Percentage of Loci with More than Two Alleles",xlab="",ylab="",las=0)
mtext("Percentage of Loci",side=2,line=2,cex=.75)



