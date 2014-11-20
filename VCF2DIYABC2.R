############Convert vcf into a dataframe for DIYABC.
require(VariantAnnotation)
require(snpStats)

vcffile <- "~/Dropbox/HFRad-Tags/HF_m4M3n5_Results/populations_sepphylo_m10r.75p9_only1snp/PhyloPop_m10r.75p9_only1snp_translated.vcf"
indmapfile <- "~/Dropbox/HFRad-Tags/Indiv_m_ConversionTable.csv"


vcf <- readVcf(vcffile,genome="taeGut1")
indmap <- read.csv(indmapfile)


chrom <- as.character(seqnames(rowData(vcf)))


calls <- geno(vcf)$GT
rownames(calls) <- chrom
a0 <- ref(vcf)
a1 <- alt(vcf)

#Take out all Z alleles
names(rownames(calls))<-NULL
Zs <- rownames(calls)!="Z"
callsAuto<- calls[Zs,]
a0Auto <- a0[Zs]
a1Auto <- a1[Zs]

snpMAllAuto <- genotypeToSnpMatrix(x=callsAuto,ref=a0Auto,alt=a1Auto)
snpMAuto <- snpMAllAuto$genotypes

simpleSnp <- matrix(as.numeric(snpMAuto),nrow=nrow(snpMAuto),ncol=ncol(snpMAuto))
rownames(simpleSnp) <- rownames(snpMAuto)
colnames(simpleSnp) <- colnames(snpMAuto)
simpleSnp <- as.data.frame(simpleSnp)
simpleSnp0s <- apply(simpleSnp,2,function(x) sapply(x, function(y) y == 0))
simpleSnp <- apply(simpleSnp,2,function(x) sapply(x,function(y) y-1))
simpleSnp[simpleSnp0s] <- 9

autochr <- rep("A",ncol(simpleSnp))
simpleSnpAlabel <- rbind(autochr,simpleSnp)

indphylo <- c("POP",paste("P",indmap[indmap$Ind %in% rownames(simpleSnp),"PMBasic"],sep=""))
indsex <- c("SEX",as.character(indmap[indmap$Ind %in% rownames(simpleSnp),"Sex"]))

total <- cbind(indsex,indphylo,simpleSnpAlabel)
rownames(total)[1] <- "IND"
total <-

totalsorted <- total[order(total[,"indphylo"]),]
totalsorted <- rbind(totalsorted[nrow(totalsorted),],totalsorted[1:nrow(totalsorted)-1,])
totalsorted <- cbind(rownames(totalsorted),totalsorted)
totalsorted[1,1] <- "IND"

write.table(totalsorted,"~/Dropbox/HFRad-Tags/HF_m4M3n5_Results/populations_sepphylo_m10r.75p9_only1snp/_AutosomalDIYABCInput.txt",col.names=F,quote=F,row.names=F)







########The code below can be used to exclude all individuals not from Pops 1-3, and delete all loci that have all missing data in one of those populations.

totalsorted <- as.data.frame(totalsorted)
pops123 <- totalsorted[totalsorted$indphylo %in% c("P1","P2","P3"),]

keep <- rep(T,ncol(pops123))
for (i in 4:ncol(pops123)){
	testing <- as.numeric(as.character(pops123[,i]))
	grpmeans <- aggregate(testing,by=list(as.character(pops123$indphylo)),FUN="mean")
	if (9 %in% grpmeans$x){
		keep[i] <- F
	}
}

goodloci <- length(keep[keep==T])
length(keep[keep==F])

pops123 <- pops123[,keep]

firstline <- rep("A",ncol(pops123))
firstline[1] <- "IND"
firstline[2] <- "SEX"
firstline[3] <- "POP"

snpfile <- rbind(firstline,as.matrix(pops123))

write.table(snpfile,"~/Dropbox/HFRad-Tags/HF_m4M3n5_Results/populations_sepphylo_m10r.75p9_only1snp/_AutosomalDIYABCInput_nooutgroupsmissing.txt",col.names=F,quote=F,row.names=F)


#Code to pull out a random number of snps
#randomsnps <- c(1,2,3,runif(2000,min=4,max=ncol(totalsorted)))
#totalsortedsmall <- totalsorted[,randomsnps]
#write.table(totalsortedsmall,"_AutosomalDIYABCInput_2000Loci.snp",col.names=F,quote=F)
