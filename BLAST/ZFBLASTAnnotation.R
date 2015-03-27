require(biomaRt)
require(GenomicRanges)
require(rtracklayer)
require(GenomicFeatures)
require(VariantAnnotation)
require(biovizBase)
require(ggbio)
require(snpStats)
require(BSgenome.Tguttata.UCSC.taeGut1)

vcfFile = "~/Dropbox/HFRad-Tags/HF_m4M3n5_Results/populations_sepphylo_m10r.75p9/PhyloPop_m10r.75p9_translated.vcf"

vcfFile = "~/Dropbox/HFRad-Tags/HF_m4M3n5_Results/populations_sepphylo_m10r.75p9_only1snp/PhyloPop_m10r.75p9_only1snp_translated.vcf"

vcfFile = "~/Dropbox/ElisaAlbatross/populations_m=30/ZebraFinch_e-10BlastTranslated.vcf"

indnames <- "~/Dropbox/HFRad-Tags/Indiv_m_ConversionTable.csv"

geocolors <- c("#7570B3","#1B9E77","#D95F02")
geolabs <- c("West","East","Hawaii")
names(geocolors) <- geolabs

#Open blast results from blastn, see how many hits each fragment had
# blast <- read.csv("~/Dropbox/HFRad-Tags/PhyloPopr.5p14BLAST/H9BR64VU015-Alignment-HitTable.csv",header=F)

# fields <- c("queryid", "subjectids", "%identity", "alignmentlength", "mismatches", "gapopens", "q.start", "q.end", "s.start", "s.end", "evalue", "bitscore")
# colnames(blast) <- fields

# numhits <- as.matrix(table(blast$queryid))



#This is a summary of how many hits per framgment, any that have more than 1 hit will be discarded for future analysis
#table(numhits)


#Use biomaRt to extract the zebra finch gene database.
ensembl <- useMart("ensembl",dataset="tguttata_gene_ensembl")

attributes <- listAttributes(ensembl)

ZFgenes <- getBM(attributes= c("ensembl_gene_id","external_gene_id","description","chromosome_name","start_position","end_position"),mart=ensembl)

ZFchr <- getChromInfoFromBiomart(biomart="ensembl",dataset="tguttata_gene_ensembl")

ZFchrlen <- ZFchr[,2]
names(ZFchrlen) <- ZFchr[,1]

#In the case that chrominfo doesn't work above:
zeb <- getIdeogram(genome="taeGut1",cytoband=F)
zeb <- keepSeqlevels(zeb,paste0("chr",c("1","1A","1B","2","3","4","4A",5:28,"LG2","LG5","LGE22","Z")))
ZFchrlen <- seqlengths(zeb)
names(ZFchrlen) <- c("1","1A","1B","2","3","4","4A",5:28,"LG2","LG5","LGE22","Z")



#ZFanno <- GRanges(seqnames=Rle(ZFgenes$chromosome_name),ranges=IRanges(start=ZFgenes$start_position,end=ZFgenes$end_position),gene=ZFgenes$ensembl_gene_id,gene_id=ZFgenes$external_gene_id,seqlengths=ZFchrlen)

#It may actually be more useful to have a TranscriptDB oject to annotate variants. Note that this contains much of the similar info that I extracted previously, but with additional info.

ZFtdb <- makeTranscriptDbFromBiomart(biomart="ensembl",dataset="tguttata_gene_ensembl")

#############################################################################
#Open VCF file with ZF genome positions

HFvcf <- readVcf(vcfFile,genome="taeGut1")
seqlengths(HFvcf) <- ZFchrlen[names(seqlengths(HFvcf))]
#autoplot(HFvcf,layout="karyogram")

#Average sequencing depth:
hist(geno(HFvcf)$DP,breaks=300)
depth <- geno(HFvcf)$DP
loci <- unique(depth)
zeros <- as.numeric(loci)
nozeros <- zeros[zeros!=0]
mean(nozeros)
sd(nozeros)
hist(nozeros,breaks=100,main="Sequencing Depth",xlab="Depth")
boxplot(nozeros)
median(nozeros)
sterr <- sqrt(var(nozeros)/length(nozeros))

#Produce an ideogram
singlesnpvcfFile = "~/Dropbox/HFRad-Tags/HF_m4M3n5_Results/populations_sepphylo_m10r.75p9_only1snp/PhyloPop_m10r.75p9_only1snp_translated.vcf"
singlesnpvcfFile = "~/Dropbox/ElisaAlbatross/populations_m=30/ZebraFinch_e-10BlastTranslated_SingleSnp.vcf"
singlesnpvcfFile = "~/Dropbox/ElisaAlbatross/populations_m=30/Chicken_e-20BlastTranslated_SingleSnp.vcf"
snps <- readVcfLongForm(singlesnpvcfFile,genome="taeGut1")
seqlengths(snps) <- ZFchrlen[names(seqlengths(snps))]
#The above no longer works in the latest version of bioconductor, so have to use the "expand" function.  But, it doesn't work if snpStats is also loaded.
snps <- expand(HFvcf)
seqlengths(snps) <- ZFchrlen[names(seqlengths(snps))]
autoplot(snps,layout="karyogram")
snps <- keepSeqlevels(snps,c("1","1A","2","3","4","4A",5:15,17:24,26:28,"Z"))
p <- autoplot(snps,layout="karyogram")

#How many from Z chromosome
z <- as.character(runValue(seqnames(snps)))
length(z[z=="Z"])
x <- rep(z,runLength(seqnames(snps)))
tx <- table(x)[names(seqlengths(snps))]
lmxz <- lm(tx~seqlengths(snps))
y <- seqlengths(snps)
plot(y,tx,ylab="Number of loci",xlab="Chromosome length",pch=20,axes=F)
points(seqlengths(snps)["Z"],tx["Z"],bg="red",pch=21)
abline(lmxz)
box()
axis(side=2)
axis(side=1)

#Annotate vcf variants, and predict if syn or nonsyn for coding variants.
#Note that this was neccessary before the BSGenome package was available for ZF, but no longer is necessary. So this is the old way.
ZFfa <- open(FaFile(file="~/Stacks/ZebraFinchGenome/SamTools/AllChr.fasta"))
coding <- predictCoding(HFvcf,ZFtdb,ZFfa)


#Need to remove the "chr" from the seqenames of Tguttata so they match.
seqnames(Tguttata) <- sapply(seqnames(Tguttata),function(x) substring(x,4))

rd <- rowData(HFvcf)

#Can either locate all variants at one time, or break them out into different types.
allvar <- locateVariants(rd,ZFtdb,AllVariants())


codingvar <- locateVariants(rd,ZFtdb,CodingVariants())
introns <- locateVariants(rd,ZFtdb,IntronVariants())
intergenic <- locateVariants(rd,ZFtdb,IntergenicVariants())
fiveutr <- locateVariants(rd,ZFtdb,FiveUTRVariants())
threeutr <- locateVariants(rd,ZFtdb,ThreeUTRVariants())
splice <- locateVariants(rd,ZFtdb,SpliceSiteVariants())
promoter <- locateVariants(rd,ZFtdb,PromoterVariants())

coding <- predictCoding(HFvcf,ZFtdb,seqSource=Tguttata)


locs <- unique(mcols(allvar)$LOCATION)

snplocs <- mcols(allvar)$LOCATION

#How you can subset a GRanges object by location.
intronsallvar <- allvar[as.character(snplocs) == "intron",]

#Subset the VCF object by a granges object
intronhits <- findOverlaps(rowData(HFvcf),introns)
HFintrons <- HFvcf[unique(queryHits(intronhits)),]

intergenichits <- findOverlaps(rowData(HFvcf),intergenic)
HFintergenic <- HFvcf[unique(queryHits(intergenichits)),]

codingvarhits <- findOverlaps(rowData(HFvcf),codingvar)
HFcodingvar <- HFvcf[unique(queryHits(codingvarhits)),]

#Note that there are some overlaps with promoters and introns.
promoterhits <- findOverlaps(rowData(HFvcf),promoter)
HFpromoter <- HFvcf[unique(queryHits(promoterhits)),]





#Output coding variants.
codingdf <- as.data.frame(coding,row.names=1:length(coding))
codingdf$loci_id <- names(coding)
write.csv(codingdf,"~/Dropbox/HFRad-Tags/HF_m4M3n5_Results/populations_sepphylo_m10r.75p9/CodingRegions.csv")





#Create a SnpMatrix datatype from the vcf file.
#Playing around with the SFS for my subsampled individuals:

# HFvcf <- readVcf("~/Dropbox/HFRad-Tags/HF_m4M3n5_Results/populations_8indivperpop_10reps/poprep1/batch_1.vcf",genome="taeGut1")
# seqlengths(HFvcf) <- ZFchrlen[names(seqlengths(HFvcf))]

snpMAll <- genotypeToSnpMatrix(HFvcf)
snpM <- snpMAll$genotypes
snpsum <- col.summary(snpM)

indmap <- read.csv(indnames)
#indmap$Ind <- sapply(indmap$Ind,function(x) paste(x,"_m",sep=""))

westerninds <- as.character(indmap[indmap$PMBasic==1,"Ind"])
easterninds <- as.character(indmap[indmap$PMBasic==2,"Ind"])
hawaiiinds <- as.character(indmap[indmap$PMBasic==3,"Ind"])

outgroups <- as.character(indmap[indmap$PMBasic==6,"Ind"])

#Added option to subsample only 8 individuals
# westerninds <- sample(westerninds,8)
# easterninds <- sample(easterninds,8)
# hawaiiinds <- sample(hawaiiinds,8)

#Create a dataset with only individuals from three populations, so that we can remove all snps that are invariant within these three populations.
phyloinds <- c(westerninds,easterninds,hawaiiinds)
snpMPhylo <- snpM[intersect(rownames(snpM),phyloinds)]
snpsumPhylo <- col.summary(snpMPhylo)
snpsumPhylo.var <- rownames(snpsumPhylo)[snpsumPhylo$MAF != 0]

snpM <- snpM[,snpsumPhylo.var]

snpMWest <- snpM[intersect(rownames(snpM),westerninds)]
snpMEast <- snpM[intersect(rownames(snpM),easterninds)]
snpMHawaii <- snpM[intersect(rownames(snpM),hawaiiinds)]

snpMOutgroups <- snpM[intersect(rownames(snpM),outgroups)]

#Added option to subsample only 8 individuals.
# snpMWest <- snpMWest[sample(rownames(snpMWest),8),]
# snpMEast <- snpMEast[sample(rownames(snpMEast),8),]
# snpMHawaii <- snpMHawaii[sample(rownames(snpMHawaii),8),]




snpsumWest <- col.summary(snpMWest)
snpsumEast <- col.summary(snpMEast)
snpsumHawaii <- col.summary(snpMHawaii)
snpsumOutgroups <- col.summary(snpMOutgroups)

####It is also possible that we want to remove all invariant sites from each specific population.  So, I will create datasets that do not contain them.
#Note that the below three lines should only be used if this is desired.
snpsumWest <- snpsumWest[snpsumWest$MAF !=0,]
snpsumEast <- snpsumEast[snpsumEast$MAF !=0,]
snpsumHawaii <- snpsumHawaii[snpsumHawaii$MAF !=0,]


#Plot SFS for all SNPs
geocolorsalpha <- c(rgb(117, 112, 179,alpha=255/2,maxColorValue=255),rgb(27, 158, 119,alpha=255/2,maxColorValue=255),rgb(217, 95, 2,alpha=255/2,maxColorValue=255))
names(geocolorsalpha) <- geolabs

hW <- hist(snpsumWest$MAF,breaks=30,col=geocolorsalpha["West"])
hE <- hist(snpsumEast$MAF,breaks=30,col=geocolorsalpha["East"])
hH <- hist(snpsumHawaii$MAF,breaks=30,col=geocolorsalpha["Hawaii"])

#Smaller breaks for subsampled data
hW <- hist(snpsumWest$MAF,breaks=15,col=geocolorsalpha["West"])
hE <- hist(snpsumEast$MAF,breaks=15,col=geocolorsalpha["East"])
hH <- hist(snpsumHawaii$MAF,breaks=15,col=geocolorsalpha["Hawaii"])


plot(hW,col=geocolorsalpha["West"], ylim=c(0,1000),xlab="MAF",main="Western vs. Eastern SFS" )
plot(hE,col=geocolorsalpha["East"], add=T)

#Plotting as barplots makes it much easier to compare across groups.
barplot(hW$density/sum(hW$density),col=geocolorsalpha["West"],ylim=c(0,1),ylab="Percentage of sites",xlab="MAF",main="Western vs. East SFS",names.arg=hW$mids,cex.names=.75)
barplot(hE$density/sum(hE$density),col=geocolorsalpha["East"],ylim=c(0,1),add=T)


we.test <- ks.test(snpsumWest$MAF,snpsumEast$MAF)

plot(hW,col=geocolorsalpha["West"], ylim=c(0,1000),xlab="MAF",main="Western vs. Hawaii SFS" )
plot(hH,col=geocolorsalpha["Hawaii"], add=T)


#Plotting as barplots makes it much easier to compare across groups.
barplot(hW$density/sum(hW$density),col=geocolorsalpha["West"],ylim=c(0,1),ylab="Percentage of sites",xlab="MAF",main="Western vs. Hawaii SFS",names.arg=hW$mids,cex.names=.75)
barplot(hH$density/sum(hH$density),col=geocolorsalpha["Hawaii"],ylim=c(0,1),add=T)


wh.test <- ks.test(snpsumWest$MAF,snpsumHawaii$MAF)

##########Paper Figure
layout(matrix(1:2,nrow=1))
par(mar=c(4,3,1.5,1)+0.1)
barplot(hW$density/sum(hW$density),col=geocolorsalpha["West"],ylim=c(0,0.6),ylab="Percentage of sites",xlab="Minor Allele Frequency",names.arg=hW$mids,cex.names=.75,cex=.75)
barplot(hH$density/sum(hH$density),col=geocolorsalpha["Hawaii"],ylim=c(0,0.6),add=T,axes=F)

barplot(hW$density/sum(hW$density),col=geocolorsalpha["West"],ylim=c(0,0.6),ylab="Percentage of sites",xlab="Minor Allele Frequency",names.arg=hW$mids,cex.names=.75,cex=.75)
barplot(hE$density/sum(hE$density),col=geocolorsalpha["East"],ylim=c(0,0.6),add=T,axes=F)
legend("topright",legend=c("Hawaii","Western","Eastern"),fill=c(geocolorsalpha["Hawaii"],geocolorsalpha["West"],geocolorsalpha["East"]),title="Population",cex=.6)

#Number of sites for each pop:
length(snpsumWest$MAF)
length(snpsumEast$MAF)
length(snpsumHawaii$MAF)



#Generate SFS for different types of markers



#Intergenic
snpMAllintergenic <- genotypeToSnpMatrix(HFintergenic)
snpMintergenic <- snpMAllintergenic$genotypes
snpsumintergenic <- col.summary(snpMintergenic)

indmap <- read.csv(indnames)
#indmap$Ind <- sapply(indmap$Ind,function(x) paste(x,"_m",sep=""))

westerninds <- as.character(indmap[indmap$PMBasic==1,"Ind"])
easterninds <- as.character(indmap[indmap$PMBasic==2,"Ind"])
hawaiiinds <- as.character(indmap[indmap$PMBasic==3,"Ind"])

snpMWestintergenic <- snpMintergenic[intersect(rownames(snpMintergenic),westerninds)]
snpMEastintergenic <- snpMintergenic[intersect(rownames(snpMintergenic),easterninds)]
snpMHawaiiintergenic <- snpMintergenic[intersect(rownames(snpMintergenic),hawaiiinds)]

#Subsample 8 individuals
# snpMWestintergenic <- snpMWestintergenic[sample(rownames(snpMWestintergenic),8),]
# snpMEastintergenic <- snpMEastintergenic[sample(rownames(snpMEastintergenic),8),]
# snpMHawaiiintergenic <- snpMHawaiiintergenic[sample(rownames(snpMHawaiiintergenic),8),]

snpsumWestintergenic <- col.summary(snpMWestintergenic)
snpsumEastintergenic <- col.summary(snpMEastintergenic)
snpsumHawaiiintergenic <- col.summary(snpMHawaiiintergenic)


####It is also possible that we want to remove all invariant sites from each specific population.  So, I will create datasets that do not contain them.
#Note that the below three lines should only be used if this is desired.
snpsumWestintergenic <- snpsumWestintergenic[snpsumWestintergenic$MAF !=0,]
snpsumEastintergenic <- snpsumEastintergenic[snpsumEastintergenic$MAF !=0,]
snpsumHawaiiintergenic <- snpsumHawaiiintergenic[snpsumHawaiiintergenic$MAF !=0,]



#Plot SFS for all SNPs
geocolorsalpha <- c(rgb(117, 112, 179,alpha=255/2,maxColorValue=255),rgb(27, 158, 119,alpha=255/2,maxColorValue=255),rgb(217, 95, 2,alpha=255/2,maxColorValue=255))
names(geocolorsalpha) <- geolabs

hWintergenic <- hist(snpsumWestintergenic$MAF,breaks=30,col=geocolorsalpha["West"])
hEintergenic <- hist(snpsumEastintergenic$MAF,breaks=30,col=geocolorsalpha["East"])
hHintergenic <- hist(snpsumHawaiiintergenic$MAF,breaks=30,col=geocolorsalpha["Hawaii"])

#Smaller breaks for subsampled data
hWintergenic <- hist(snpsumWestintergenic$MAF,breaks=c(0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50),col=geocolorsalpha["West"])
hEintergenic <- hist(snpsumEastintergenic$MAF,breaks=c(0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50),col=geocolorsalpha["East"])
hHintergenic <- hist(snpsumHawaiiintergenic$MAF,breaks=c(0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50),col=geocolorsalpha["Hawaii"])

plot(hWintergenic,col=geocolorsalpha["West"], ylim=c(0,8000),xlab="MAF",main="Western vs. Eastern SFS" )
plot(hEintergenic,col=geocolorsalpha["East"], add=T)

barplot(hWintergenic$density/sum(hWintergenic$density),col=geocolorsalpha["West"],ylim=c(0,1),ylab="Percentage of sites",xlab="MAF",main="Western vs. East SFS")
barplot(hEintergenic$density/sum(hEintergenic$density),col=geocolorsalpha["East"],ylim=c(0,1),add=T)

we.test <- ks.test(snpsumWestintergenic$MAF,snpsumEastintergenic$MAF)

plot(hWintergenic,col=geocolorsalpha["West"], ylim=c(0,8000),xlab="MAF",main="Western vs. Hawaii SFS" )
plot(hHintergenic,col=geocolorsalpha["Hawaii"], add=T)

barplot(hWintergenic$density/sum(hWintergenic$density),col=geocolorsalpha["West"],ylim=c(0,1),ylab="Percentage of sites",xlab="MAF",main="Western vs. Hawaii SFS")
barplot(hHintergenic$density/sum(hHintergenic$density),col=geocolorsalpha["Hawaii"],ylim=c(0,1),add=T)


wh.test <- ks.test(snpsumWestintergenic$MAF,snpsumHawaiiintergenic$MAF)


##########Paper Figure
layout(matrix(1:2,nrow=1))
barplot(hWintergenic$density/sum(hWintergenic$density),col=geocolorsalpha["West"],ylim=c(0,1),ylab="Percentage of sites",xlab="Minor Allele Frequency",names.arg=hWintergenic$mids,cex.names=.75,main="Intergenic Sites")
barplot(hHintergenic$density/sum(hHintergenic$density),col=geocolorsalpha["Hawaii"],ylim=c(0,1),add=T)

barplot(hWintergenic$density/sum(hWintergenic$density),col=geocolorsalpha["West"],ylim=c(0,1),ylab="Percentage of sites",xlab="Minor Allele Frequency",names.arg=hWintergenic$mids,cex.names=.75,main="Intergenic Sites")
barplot(hEintergenic$density/sum(hEintergenic$density),col=geocolorsalpha["East"],ylim=c(0,1),add=T)
legend("topright",legend=c("Western","Eastern","Hawaii"),fill=geocolorsalpha,title="Population")

#Number of sites for each pop:
length(snpsumWestintergenic$MAF)
length(snpsumEastintergenic$MAF)
length(snpsumHawaiiintergenic$MAF)




#Coding Variants

snpMAllcodingvar <- genotypeToSnpMatrix(HFcodingvar)
snpMcodingvar <- snpMAllcodingvar$genotypes
snpsumcodingvar <- col.summary(snpMcodingvar)

snpMWestcodingvar <- snpMcodingvar[intersect(rownames(snpMcodingvar),westerninds)]
snpMEastcodingvar <- snpMcodingvar[intersect(rownames(snpMcodingvar),easterninds)]
snpMHawaiicodingvar <- snpMcodingvar[intersect(rownames(snpMcodingvar),hawaiiinds)]


#####To subsample 8 individuals per population
# snpMWestcodingvar <- snpMWestcodingvar[sample(rownames(snpMWestcodingvar),8),]
# snpMEastcodingvar <- snpMEastcodingvar[sample(rownames(snpMEastcodingvar),8),]
# snpMHawaiicodingvar <- snpMHawaiicodingvar[sample(rownames(snpMHawaiicodingvar),8),]

snpsumWestcodingvar <- col.summary(snpMWestcodingvar)
snpsumEastcodingvar <- col.summary(snpMEastcodingvar)
snpsumHawaiicodingvar <- col.summary(snpMHawaiicodingvar)


####It is also possible that we want to remove all invariant sites from each specific population.  So, I will create datasets that do not contain them.
#Note that the below three lines should only be used if this is desired.
snpsumWestcodingvar <- snpsumWestcodingvar[snpsumWestcodingvar$MAF !=0,]
snpsumEastcodingvar <- snpsumEastcodingvar[snpsumEastcodingvar$MAF !=0,]
snpsumHawaiicodingvar <- snpsumHawaiicodingvar[snpsumHawaiicodingvar$MAF !=0,]



#Plot SFS for all SNPs
geocolorsalpha <- c(rgb(117, 112, 179,alpha=255/2,maxColorValue=255),rgb(27, 158, 119,alpha=255/2,maxColorValue=255),rgb(217, 95, 2,alpha=255/2,maxColorValue=255))
names(geocolorsalpha) <- geolabs

hWcodingvar <- hist(snpsumWestcodingvar$MAF,breaks=30,col=geocolorsalpha["West"])
hEcodingvar <- hist(snpsumEastcodingvar$MAF,breaks=30,col=geocolorsalpha["East"])
hHcodingvar <- hist(snpsumHawaiicodingvar$MAF,breaks=30,col=geocolorsalpha["Hawaii"])

#Smaller breaks for subsampled data
hWcodingvar <- hist(snpsumWestcodingvar$MAF,breaks=c(0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50),col=geocolorsalpha["West"])
hEcodingvar <- hist(snpsumEastcodingvar$MAF,breaks=c(0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50),col=geocolorsalpha["East"])
hHcodingvar <- hist(snpsumHawaiicodingvar$MAF,breaks=c(0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50),col=geocolorsalpha["Hawaii"])

plot(hWcodingvar,col=geocolorsalpha["West"], ylim=c(0,300),xlab="MAF",main="Western vs. Eastern SFS" )
plot(hEcodingvar,col=geocolorsalpha["East"], add=T)

barplot(hWcodingvar$density/sum(hWcodingvar$density),col=geocolorsalpha["West"],ylim=c(0,1),ylab="Percentage of sites",xlab="MAF",main="Western vs. East SFS")
barplot(hEcodingvar$density/sum(hEcodingvar$density),col=geocolorsalpha["East"],ylim=c(0,1),add=T)

we.test <- ks.test(snpsumWestcodingvar$MAF,snpsumEastcodingvar$MAF)

plot(hWcodingvar,col=geocolorsalpha["West"], ylim=c(0,300),xlab="MAF",main="Western vs. Hawaii SFS" )
plot(hHcodingvar,col=geocolorsalpha["Hawaii"], add=T)

barplot(hWcodingvar$density/sum(hWcodingvar$density),col=geocolorsalpha["West"],ylim=c(0,1),ylab="Percentage of sites",xlab="MAF",main="Western vs. Hawaii SFS")
barplot(hHcodingvar$density/sum(hHcodingvar$density),col=geocolorsalpha["Hawaii"],ylim=c(0,1),add=T)

wh.test <- ks.test(snpsumWestcodingvar$MAF,snpsumHawaiicodingvar$MAF)

##########Paper Figure
layout(matrix(1:2,nrow=1))
barplot(hWcodingvar$density/sum(hWcodingvar$density),col=geocolorsalpha["West"],ylim=c(0,1),ylab="Percentage of sites",xlab="Minor Allele Frequency",names.arg=hWcodingvar$mids,cex.names=.75,main="Coding Sites")
barplot(hHcodingvar$density/sum(hHcodingvar$density),col=geocolorsalpha["Hawaii"],ylim=c(0,1),add=T)

barplot(hWcodingvar$density/sum(hWcodingvar$density),col=geocolorsalpha["West"],ylim=c(0,1),ylab="Percentage of sites",xlab="Minor Allele Frequency",names.arg=hWcodingvar$mids,cex.names=.75,main="Coding Sites")
barplot(hEcodingvar$density/sum(hEcodingvar$density),col=geocolorsalpha["East"],ylim=c(0,1),add=T)
legend("topright",legend=c("Western","Eastern","Hawaii"),fill=geocolorsalpha,title="Population")

#Number of sites for each pop:
length(snpsumWestcodingvar$MAF)
length(snpsumEastcodingvar$MAF)
length(snpsumHawaiicodingvar$MAF)



# Introns

snpMAllintrons <- genotypeToSnpMatrix(HFintrons)
snpMintrons <- snpMAllintrons$genotypes
snpsumintrons <- col.summary(snpMintrons)

snpMWestintrons <- snpMintrons[intersect(rownames(snpMintrons),westerninds)]
snpMEastintrons <- snpMintrons[intersect(rownames(snpMintrons),easterninds)]
snpMHawaiiintrons <- snpMintrons[intersect(rownames(snpMintrons),hawaiiinds)]

#Subsample 8 individuals
# snpMWestintrons <- snpMWestintrons[sample(rownames(snpMWestintrons),8),]
# snpMEastintrons <- snpMEastintrons[sample(rownames(snpMEastintrons),8),]
# snpMHawaiiintrons <- snpMHawaiiintrons[sample(rownames(snpMHawaiiintrons),8),]

snpsumWestintrons <- col.summary(snpMWestintrons)
snpsumEastintrons <- col.summary(snpMEastintrons)
snpsumHawaiiintrons <- col.summary(snpMHawaiiintrons)


####It is also possible that we want to remove all invariant sites from each specific population.  So, I will create datasets that do not contain them.
#Note that the below three lines should only be used if this is desired.
snpsumWestintrons <- snpsumWestintrons[snpsumWestintrons$MAF !=0,]
snpsumEastintrons <- snpsumEastintrons[snpsumEastintrons$MAF !=0,]
snpsumHawaiiintrons <- snpsumHawaiiintrons[snpsumHawaiiintrons$MAF !=0,]



#Plot SFS for all SNPs
geocolorsalpha <- c(rgb(117, 112, 179,alpha=255/2,maxColorValue=255),rgb(27, 158, 119,alpha=255/2,maxColorValue=255),rgb(217, 95, 2,alpha=255/2,maxColorValue=255))
names(geocolorsalpha) <- geolabs

hWintrons <- hist(snpsumWestintrons$MAF,breaks=30,col=geocolorsalpha["West"])
hEintrons <- hist(snpsumEastintrons$MAF,breaks=30,col=geocolorsalpha["East"])
hHintrons <- hist(snpsumHawaiiintrons$MAF,breaks=30,col=geocolorsalpha["Hawaii"])

#Smaller breaks for subsampled data
hWintrons <- hist(snpsumWestintrons$MAF,breaks=c(0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50),col=geocolorsalpha["West"])
hEintrons <- hist(snpsumEastintrons$MAF,breaks=c(0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50),col=geocolorsalpha["East"])
hHintrons <- hist(snpsumHawaiiintrons$MAF,breaks=c(0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50),col=geocolorsalpha["Hawaii"])

plot(hWintrons,col=geocolorsalpha["West"], ylim=c(0,700),xlab="MAF",main="Western vs. Eastern SFS" )
plot(hEintrons,col=geocolorsalpha["East"], add=T)

barplot(hWintrons$density/sum(hWintrons$density),col=geocolorsalpha["West"],ylim=c(0,1),ylab="Percentage of sites",xlab="MAF",main="Western vs. East SFS")
barplot(hEintrons$density/sum(hEintrons$density),col=geocolorsalpha["East"],ylim=c(0,1),add=T)

we.test <- ks.test(snpsumWestintrons$MAF,snpsumEastintrons$MAF)

plot(hWintrons,col=geocolorsalpha["West"], ylim=c(0,700),xlab="MAF",main="Western vs. Hawaii SFS" )
plot(hHintrons,col=geocolorsalpha["Hawaii"], add=T)

barplot(hWintrons$density/sum(hWintrons$density),col=geocolorsalpha["West"],ylim=c(0,1),ylab="Percentage of sites",xlab="MAF",main="Western vs. Hawaii SFS")
barplot(hHintrons$density/sum(hHintrons$density),col=geocolorsalpha["Hawaii"],ylim=c(0,1),add=T)

wh.test <- ks.test(snpsumWestintrons$MAF,snpsumHawaiiintrons$MAF)

##########Paper Figure
layout(matrix(1:2,nrow=1))
barplot(hWintrons$density/sum(hWintrons$density),col=geocolorsalpha["West"],ylim=c(0,1),ylab="Percentage of sites",xlab="Minor Allele Frequency",names.arg=hWintrons$mids,cex.names=.75,main="Intron Sites")
barplot(hHcodingvar$density/sum(hHcodingvar$density),col=geocolorsalpha["Hawaii"],ylim=c(0,1),add=T)

barplot(hWintrons$density/sum(hWintrons$density),col=geocolorsalpha["West"],ylim=c(0,1),ylab="Percentage of sites",xlab="Minor Allele Frequency",names.arg=hWintrons$mids,cex.names=.75,main="Intron Sites")
barplot(hEintrons$density/sum(hEintrons$density),col=geocolorsalpha["East"],ylim=c(0,1),add=T)
legend("topright",legend=c("Western","Eastern","Hawaii"),fill=geocolorsalpha,title="Population")

#Number of sites for each pop:
length(snpsumWestintrons$MAF)
length(snpsumEastintrons$MAF)
length(snpsumHawaiiintrons$MAF)




#Promoter Regions

snpMAllpromoter <- genotypeToSnpMatrix(HFpromoter)
snpMpromoter <- snpMAllpromoter$genotypes
snpsumpromoter <- col.summary(snpMpromoter)

snpMWestpromoter <- snpMpromoter[intersect(rownames(snpMpromoter),westerninds)]
snpMEastpromoter <- snpMpromoter[intersect(rownames(snpMpromoter),easterninds)]
snpMHawaiipromoter <- snpMpromoter[intersect(rownames(snpMpromoter),hawaiiinds)]

#Subsample 8 individuals
# snpMWestpromoter <- snpMWestpromoter[sample(rownames(snpMWestpromoter),8),]
# snpMEastpromoter <- snpMEastpromoter[sample(rownames(snpMEastpromoter),8),]
# snpMHawaiipromoter <- snpMHawaiipromoter[sample(rownames(snpMHawaiipromoter),8),]

snpsumWestpromoter <- col.summary(snpMWestpromoter)
snpsumEastpromoter <- col.summary(snpMEastpromoter)
snpsumHawaiipromoter <- col.summary(snpMHawaiipromoter)


####It is also possible that we want to remove all invariant sites from each specific population.  So, I will create datasets that do not contain them.
#Note that the below three lines should only be used if this is desired.
snpsumWestpromoter <- snpsumWestpromoter[snpsumWestpromoter$MAF !=0,]
snpsumEastpromoter <- snpsumEastpromoter[snpsumEastpromoter$MAF !=0,]
snpsumHawaiipromoter <- snpsumHawaiipromoter[snpsumHawaiipromoter$MAF !=0,]



#Plot SFS for all SNPs
geocolorsalpha <- c(rgb(117, 112, 179,alpha=255/2,maxColorValue=255),rgb(27, 158, 119,alpha=255/2,maxColorValue=255),rgb(217, 95, 2,alpha=255/2,maxColorValue=255))
names(geocolorsalpha) <- geolabs

hWpromoter <- hist(snpsumWestpromoter$MAF,breaks=30,col=geocolorsalpha["West"])
hEpromoter <- hist(snpsumEastpromoter$MAF,breaks=30,col=geocolorsalpha["East"])
hHpromoter <- hist(snpsumHawaiipromoter$MAF,breaks=30,col=geocolorsalpha["Hawaii"])

#Smaller breaks for subsampled data
hWpromoter <- hist(snpsumWestpromoter$MAF,breaks=c(0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50),col=geocolorsalpha["West"])
hEpromoter <- hist(snpsumEastpromoter$MAF,breaks=c(0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50),col=geocolorsalpha["East"])
hHpromoter <- hist(snpsumHawaiipromoter$MAF,breaks=c(0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50),col=geocolorsalpha["Hawaii"])

plot(hWpromoter,col=geocolorsalpha["West"], ylim=c(0,70),xlab="MAF",main="Western vs. Eastern SFS" )
plot(hEpromoter,col=geocolorsalpha["East"], add=T)

barplot(hWpromoter$density/sum(hWpromoter$density),col=geocolorsalpha["West"],ylim=c(0,1),ylab="Percentage of sites",xlab="MAF",main="Western vs. East SFS")
barplot(hEpromoter$density/sum(hEpromoter$density),col=geocolorsalpha["East"],ylim=c(0,1),add=T)

we.test <- ks.test(snpsumWestpromoter$MAF,snpsumEastpromoter$MAF)

plot(hWpromoter,col=geocolorsalpha["West"], ylim=c(0,70),xlab="MAF",main="Western vs. Hawaii SFS" )
plot(hHpromoter,col=geocolorsalpha["Hawaii"], add=T)

barplot(hWpromoter$density/sum(hWpromoter$density),col=geocolorsalpha["West"],ylim=c(0,1),ylab="Percentage of sites",xlab="MAF",main="Western vs. Hawaii SFS")
barplot(hHpromoter$density/sum(hHpromoter$density),col=geocolorsalpha["Hawaii"],ylim=c(0,1),add=T)

wh.test <- ks.test(snpsumWestpromoter$MAF,snpsumHawaiipromoter$MAF)


##########Paper Figure
layout(matrix(1:2,nrow=1))
barplot(hWpromoter$density/sum(hWpromoter$density),col=geocolorsalpha["West"],ylim=c(0,1),ylab="Percentage of sites",xlab="Minor Allele Frequency",names.arg=hWpromoter$mids,cex.names=.75,main="Promoter Sites")
barplot(hHcodingvar$density/sum(hHcodingvar$density),col=geocolorsalpha["Hawaii"],ylim=c(0,1),add=T)

barplot(hWpromoter$density/sum(hWpromoter$density),col=geocolorsalpha["West"],ylim=c(0,1),ylab="Percentage of sites",xlab="Minor Allele Frequency",names.arg=hWpromoter$mids,cex.names=.75,main="Promoter Sites")
barplot(hEpromoter$density/sum(hEpromoter$density),col=geocolorsalpha["East"],ylim=c(0,1),add=T)
legend("topright",legend=c("Western","Eastern","Hawaii"),fill=geocolorsalpha,title="Population")

#Number of sites for each pop:
length(snpsumWestpromoter$MAF)
length(snpsumEastpromoter$MAF)
length(snpsumHawaiipromoter$MAF)



#Plotting to compare different types of sites in the three different populations.

layout(matrix(1:12,nrow=4))
barplot(hHintrons$density/sum(hHintrons$density),col=rgb(1,0,0,.4),add=F)
barplot(hHintergenic$density/sum(hHintergenic$density),col=rgb(0,1,0,.4),ylim=c(0,1),add=F)
barplot(hHpromoter$density/sum(hHpromoter$density),col= rgb(0,0,1,.4),ylim=c(0,1),ylab="Percentage of sites",xlab="MAF",main="Hawaii SFS")
barplot(hHcodingvar$density/sum(hHcodingvar$density),col="black",ylim=c(0,1),add=T,density=20,border=T)
legend("topright",c("Promoters","Introns","Intergenic","Coding"),fill=c(rgb(0,0,1,.4),rgb(1,0,0,.4),rgb(0,1,0,.4),"white"))
legend("topright",c("Promoters","Introns","Intergenic","Coding"),density=c(0,0,0,20))

barplot(hWpromoter$density/sum(hWpromoter$density),col= rgb(0,0,1,.4),ylim=c(0,1),ylab="Percentage of sites",xlab="MAF",main="Western SFS")
barplot(hWintrons$density/sum(hWintrons$density),col=rgb(1,0,0,.4),add=T)
barplot(hWintergenic$density/sum(hWintergenic$density),col=rgb(0,1,0,.4),ylim=c(0,1),add=T)
barplot(hWcodingvar$density/sum(hWcodingvar$density),col="black",ylim=c(0,1),add=T,density=20,border=T)
legend("topright",c("Promoters","Introns","Intergenic","Coding"),fill=c(rgb(0,0,1,.4),rgb(1,0,0,.4),rgb(0,1,0,.4),"white"))
legend("topright",c("Promoters","Introns","Intergenic","Coding"),density=c(0,0,0,20))


barplot(hEpromoter$density/sum(hEpromoter$density),col= rgb(0,0,1,.4),ylim=c(0,1),ylab="Percentage of sites",xlab="MAF",main="Eastern SFS")
barplot(hEintrons$density/sum(hEintrons$density),col=rgb(1,0,0,.4),add=T)
barplot(hEintergenic$density/sum(hEintergenic$density),col=rgb(0,1,0,.4),ylim=c(0,1),add=T)
barplot(hEcodingvar$density/sum(hEcodingvar$density),col="black",ylim=c(0,1),add=T,density=20,border=T)
legend("topright",c("Promoters","Introns","Intergenic","Coding"),fill=c(rgb(0,0,1,.4),rgb(1,0,0,.4),rgb(0,1,0,.4),"white"))
legend("topright",c("Promoters","Introns","Intergenic","Coding"),density=c(0,0,0,20))




#Plotting to compare different types of sites in the three different populations.

layout(matrix(1:12,nrow=4))
par(mar=c(1,4,3,1)+0.1)
barplot(hHintergenic$density/sum(hHintergenic$density),col=geocolorsalpha["Hawaii"],ylim=c(0,.6),add=F,main="Intergenic",ylab="Percentage of sites")
mtext(paste("N=",length(snpsumHawaiiintergenic$MAF),sep=""),side=3,cex=.8)
barplot(hHintrons$density/sum(hHintrons$density),col=geocolorsalpha["Hawaii"],add=F,main="Introns",ylim=c(0,.6),ylab="Percentage of sites")
mtext(paste("N=",length(snpsumHawaiiintrons$MAF),sep=""),side=3,cex=.8)
barplot(hHpromoter$density/sum(hHpromoter$density),col= geocolorsalpha["Hawaii"],ylim=c(0,.6),main="Promoters",ylab="Percentage of sites")
mtext(paste("N=",length(snpsumHawaiipromoter$MAF),sep=""),side=3,cex=.8)
#par(mar=c(5,4,3,1)+0.1)
barplot(hHcodingvar$density/sum(hHcodingvar$density),col=geocolorsalpha["Hawaii"],ylim=c(0,.6),add=F,xlab="Minor Allele Frequency",main="Coding",ylab="Percentage of sites",names.arg=hHcodingvar$mids)
mtext(paste("N=",length(snpsumHawaiicodingvar$MAF),sep=""),side=3,cex=.8)


barplot(hWintergenic$density/sum(hWintergenic$density),col=geocolorsalpha["West"],ylim=c(0,.6),add=F,main="Intergenic")
mtext(paste("N=",length(snpsumWestintergenic$MAF),sep=""),side=3,cex=.8)
barplot(hWintrons$density/sum(hWintrons$density),col=geocolorsalpha["West"],add=F,main="Introns",ylim=c(0,.6))
mtext(paste("N=",length(snpsumWestintrons$MAF),sep=""),side=3,cex=.8)
barplot(hWpromoter$density/sum(hWpromoter$density),col= geocolorsalpha["West"],ylim=c(0,.6),main="Promoters")
mtext(paste("N=",length(snpsumWestpromoter$MAF),sep=""),side=3,cex=.8)
#par(mar=c(5,4,4,2)+0.1)
barplot(hWcodingvar$density/sum(hWcodingvar$density),col=geocolorsalpha["West"],ylim=c(0,.6),add=F,xlab="Minor Allele Frequency",main="Coding",names.arg=hWcodingvar$mids)
mtext(paste("N=",length(snpsumWestcodingvar$MAF),sep=""),side=3,cex=.8)

barplot(hEintergenic$density/sum(hEintergenic$density),col=geocolorsalpha["East"],ylim=c(0,.6),add=F,main="Intergenic")
mtext(paste("N=",length(snpsumEastintergenic$MAF),sep=""),side=3,cex=.8)
legend("topright",legend=c("Hawaii","Western","Eastern"),fill=c(geocolorsalpha["Hawaii"],geocolorsalpha["West"],geocolorsalpha["East"]),title="Population")
barplot(hEintrons$density/sum(hEintrons$density),col=geocolorsalpha["East"],add=F,main="Introns",ylim=c(0,.6))
mtext(paste("N=",length(snpsumEastintrons$MAF),sep=""),side=3,cex=.8)
barplot(hEpromoter$density/sum(hEpromoter$density),col= geocolorsalpha["East"],ylim=c(0,.6),main="Promoters")
mtext(paste("N=",length(snpsumEastpromoter$MAF),sep=""),side=3,cex=.8)
#par(mar=c(5,4,4,2)+0.1)
barplot(hEcodingvar$density/sum(hEcodingvar$density),col=geocolorsalpha["East"],ylim=c(0,.6),add=F,xlab="Minor Allele Frequency",main="Coding",names.arg=hEcodingvar$mids)
mtext(paste("N=",length(snpsumEastcodingvar$MAF),sep=""),side=3,cex=.8)



#West Type Comparisons
ks.test(snpsumWestpromoter$MAF,snpsumWestcodingvar$MAF)
ks.test(snpsumWestintergenic$MAF,snpsumWestcodingvar$MAF)
ks.test(snpsumWestintrons$MAF,snpsumWestcodingvar$MAF)
ks.test(snpsumWestintrons$MAF,snpsumWestintergenic$MAF)
ks.test(snpsumWestpromoter$MAF,snpsumWestintergenic$MAF)
ks.test(snpsumWestintrons$MAF,snpsumWestpromoter$MAF)

#East Type Comparisons
ks.test(snpsumEastpromoter$MAF,snpsumEastcodingvar$MAF)
ks.test(snpsumEastintergenic$MAF,snpsumEastcodingvar$MAF)
ks.test(snpsumEastintrons$MAF,snpsumEastcodingvar$MAF)
ks.test(snpsumEastintrons$MAF,snpsumEastintergenic$MAF)
ks.test(snpsumEastpromoter$MAF,snpsumEastintergenic$MAF)
ks.test(snpsumEastintrons$MAF,snpsumEastpromoter$MAF)

#Hawaii Type Comparisons
ks.test(snpsumHawaiipromoter$MAF,snpsumHawaiicodingvar$MAF)
ks.test(snpsumHawaiiintergenic$MAF,snpsumHawaiicodingvar$MAF)
ks.test(snpsumHawaiiintrons$MAF,snpsumHawaiicodingvar$MAF)
ks.test(snpsumHawaiiintrons$MAF,snpsumHawaiiintergenic$MAF)
ks.test(snpsumHawaiipromoter$MAF,snpsumHawaiiintergenic$MAF)
ks.test(snpsumHawaiiintrons$MAF,snpsumHawaiipromoter$MAF)





#Only intergenic and coding.
layout(matrix(1:3,nrow=1))

barplot(hHcodingvar$density/sum(hHcodingvar$density),col= geocolorsalpha["Hawaii"],ylim=c(0,1),ylab="Percentage of sites",xlab="MAF",main="Hawaii SFS")
legend("topright",legend=c("Coding","Intergenic"),fill=c(geocolorsalpha["Hawaii"],geocolorsalpha["Hawaii"]))
legend("topright",legend=c("Coding","Intergenic"),density=c(20,0))
barplot(hHcodingvar$density/sum(hHcodingvar$density),col="black",ylim=c(0,1),add=T,density=20,border=T)
barplot(hHintergenic$density/sum(hHintergenic$density),col=geocolorsalpha["Hawaii"],ylim=c(0,1),add=T)

barplot(hWcodingvar$density/sum(hWcodingvar$density),col= geocolorsalpha["West"],ylim=c(0,1),ylab="Percentage of sites",xlab="MAF",main="Western SFS")
barplot(hWcodingvar$density/sum(hWcodingvar$density),col="black",ylim=c(0,1),add=T,density=20,border=T)
barplot(hWintrons$density/sum(hWintrons$density),col= geocolorsalpha["West"],add=T)
barplot(hWintergenic$density/sum(hWintergenic$density),col=geocolorsalpha["West"],ylim=c(0,1),add=T)
legend("topright",legend=c("Coding","Intergenic"),fill=c(geocolorsalpha["West"],geocolorsalpha["West"]))
legend("topright",legend=c("Coding","Intergenic"),density=c(20,0))

barplot(hEcodingvar$density/sum(hEcodingvar$density),col= geocolorsalpha["East"],ylim=c(0,1),ylab="Percentage of sites",xlab="MAF",main="Eastern SFS")
legend("topright",legend=c("Coding","Intergenic"),fill=c(geocolorsalpha["East"],geocolorsalpha["East"]))
legend("topright",legend=c("Coding","Intergenic"),density=c(20,0))
barplot(hEcodingvar$density/sum(hEcodingvar$density),col="black",ylim=c(0,1),add=T,density=20,border=T)
barplot(hEintergenic$density/sum(hEintergenic$density),col=geocolorsalpha["East"],ylim=c(0,1),add=T)







