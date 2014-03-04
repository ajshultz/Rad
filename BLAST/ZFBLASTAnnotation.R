require(biomaRt)
require(GenomicRanges)
require(rtracklayer)
require(GenomicFeatures)
require(VariantAnnotation)

#Open blast results from blastn, see how many hits each fragment had
blast <- read.csv("~/Dropbox/HFRad-Tags/PhyloPopr.5p14BLAST/H9BR64VU015-Alignment-HitTable.csv",header=F)

fields <- c("queryid", "subjectids", "%identity", "alignmentlength", "mismatches", "gapopens", "q.start", "q.end", "s.start", "s.end", "evalue", "bitscore")
colnames(blast) <- fields

numhits <- as.matrix(table(blast$queryid))

#This is a summary of how many hits per framgment, any that have more than 1 hit will be discarded for future analysis
table(numhits)


#Use biomaRt to extract the zebra finch gene database.
ensembl <- useMart("ensembl",dataset="tguttata_gene_ensembl")

attributes <- listAttributes(ensembl)

ZFgenes <- getBM(attributes= c("ensembl_gene_id","external_gene_id","description","chromosome_name","start_position","end_position"),mart=ensembl)

ZFchr <- getChromInfoFromBiomart(biomart="ensembl",dataset="tguttata_gene_ensembl")

ZFchrlen <- ZFchr[,2]
names(ZFchrlen) <- ZFchr[,1]

ZFanno <- GRanges(seqnames=Rle(ZFgenes$chromosome_name),ranges=IRanges(start=ZFgenes$start_position,end=ZFgenes$end_position),gene=ZFgenes$ensembl_gene_id,gene_id=ZFgenes$external_gene_id,seqlengths=ZFchrlen)

#It may actually be more useful to have a TranscriptDB oject to annotate variants. Note that this contains much of the similar info that I extracted previously, but with additional info.

ZFtdb <- makeTranscriptDbFromBiomart(biomart="ensembl",dataset="tguttata_gene_ensembl")

