setwd("~/Dropbox/HFRad-Tags/HFPaired_FinalParameterTesting/PairTestingCombinedResults")

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
boxplot(per_mult,main="Percentage of Loci with Multiple Pairs")boxplot(per_excess,ylim=c(0,0.008),main="Percentage of Loci with More than Two Alleles")

AllInds <- read.csv("~/Dropbox/HFRad-Tags/HFPaired_FinalParameterTesting/_MParameterAllLociComparions.csv",row.names=1)

AllInds <- as.data.frame(t(AllInds))

plot(AllInds$VCFTotalLoci,type="l",ylab="Number of Loci",main="Loci After Filtering",xlab="M Parameter")
points(AllInds$VCFTotalLoci)

plot(AllInds$VCFAvgPi,type="l",ylab="Pi",main="Average Subpopulation Pi value",xlab="M Parameter")
points(AllInds$VCFAvgPi)

plot(AllInds$BlastMorethan1hit,col="red",ylim=c(0,900),type="l")
lines(AllInds$VCFTooManyAlleles_m4,add=T)


