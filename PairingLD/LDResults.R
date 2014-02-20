setwd("~/Dropbox/PythonScripts/RAD")
require(ggplot2)
require(hexbin)

results <- "RLDResults/"

map <- read.table("/Users/allisonshultz/Dropbox/HFRad-Tags/HFPaired_Reduced/_PhyloPop_r.5p14_allsitesoutput/batch_1.plink.map")

pairs <- read.csv("TestResults/MatchedReadPairs.csv",header=F)

pop1 <- read.csv("TestResults/pop_1.csv",header=F)
pop2 <- read.csv("TestResults/pop_2.csv",header=F)
pop3 <- read.csv("TestResults/pop_3.csv",header=F)
pop7 <- read.csv("TestResults/pop_7.csv",header=F)

loci_pos <- map[,2]
locus <- vector()
position <- vector()
for (i in 1:length(loci_pos)){
	sep <- strsplit(as.character(loci_pos[i]),"_")
	locus[i] <- sep[[1]][1]
	position[i]<- sep[[1]][2] 
}

uniqueloci <- unique(locus)

#Once the loci are read in, create a new dataframe with only the values for SNPs on the same fragment.

samelocpop1 <- as.data.frame(matrix(nrow=nrow(pop1),ncol=ncol(pop1)))
samelocpop2 <- as.data.frame(matrix(nrow=nrow(pop2),ncol=ncol(pop2)))
samelocpop3 <- as.data.frame(matrix(nrow=nrow(pop3),ncol=ncol(pop3)))
samelocpop7 <- as.data.frame(matrix(nrow=nrow(pop7),ncol=ncol(pop7)))

for (i in 1:length(uniqueloci)){
	locrows <- locus == uniqueloci[i]
	locvals <- pop1[locrows,locrows]
	samelocpop1[locrows,locrows] <- locvals
}

for (i in 1:length(uniqueloci)){
	locrows <- locus == uniqueloci[i]
	locvals <- pop2[locrows,locrows]
	samelocpop2[locrows,locrows] <- locvals
}

for (i in 1:length(uniqueloci)){
	locrows <- locus == uniqueloci[i]
	locvals <- pop3[locrows,locrows]
	samelocpop3[locrows,locrows] <- locvals
}

for (i in 1:length(uniqueloci)){
	locrows <- locus == uniqueloci[i]
	locvals <- pop7[locrows,locrows]
	samelocpop7[locrows,locrows] <- locvals
}

#Before taking any averages, we want to remove all 0 and 9 values.  0 values occur when actually there is no snp in that population (but it still got calculated because there was a snp between populations or in a different population).  A 9 indicated that r2 could not be calcuated at that position because of mathematical issues, and was put in as a placeholder. Also, I only want to take the lower triangle of the matrix.

no0or9 <- function(x){
	if (!is.na(x)){
		if (x == 0 | x > 1){
			x <- NA
		}
	}
	return(x)
}


samelocpop1 <- as.matrix(samelocpop1)
samelocpop1half <- samelocpop1
samelocpop1half[upper.tri(samelocpop1,diag=T)] <- NA
samelocpop1halfvals <- as.data.frame(apply(samelocpop1half,2,FUN=function(x) {sapply(x,FUN=no0or9)}))

samelocpop2 <- as.matrix(samelocpop2)
samelocpop2half <- samelocpop2
samelocpop2half[upper.tri(samelocpop2,diag=T)] <- NA
samelocpop2halfvals <- as.data.frame(apply(samelocpop2half,2,FUN=function(x) {sapply(x,FUN=no0or9)}))

samelocpop3 <- as.matrix(samelocpop3)
samelocpop3half <- samelocpop3
samelocpop3half[upper.tri(samelocpop3,diag=T)] <- NA
samelocpop3halfvals <- as.data.frame(apply(samelocpop3half,2,FUN=function(x) {sapply(x,FUN=no0or9)}))

samelocpop7 <- as.matrix(samelocpop7)
samelocpop7half <- samelocpop7
samelocpop7half[upper.tri(samelocpop7,diag=T)] <- NA
samelocpop7halfvals <- as.data.frame(apply(samelocpop7half,2,FUN=function(x) {sapply(x,FUN=no0or9)}))

#Need to make a matrix of distances between SNPs located on the same tag.  The values present in the matrix should match those present in smaelocpop7halfvals, as those are the only comparisions we are interested in for this particular analysis.  

#Make numeric vector of all variant positions.   
position <- as.numeric(position)

####pop1
snpdistpop1 <- matrix(nrow=nrow(samelocpop1halfvals),ncol=ncol(samelocpop1halfvals))

for (i in 1:nrow(samelocpop1halfvals)){
	for (j in 1:ncol(samelocpop1halfvals)){
		if (!is.na(samelocpop1halfvals[i,j])){
			snpdistpop1[i,j] <- position[i] - position[j]
		}
	}
}

samelocpop1halfvalsvec <- unlist(samelocpop1halfvals)
samelocpop1vecnona <- samelocpop1halfvalsvec[!is.na(samelocpop1halfvalsvec)]
snpdistpop1vec <- unlist(snpdistpop1)
snpdistpop1vecnona <- snpdistpop1[!is.na(snpdistpop1)]

samelocpop1df <- data.frame(snpdistpop1vecnona,samelocpop1vecnona)
colnames(samelocpop1df) <- c("SNPdist","r2val")

write.csv(samelocpop1df,file=paste(results,"sameloc_pop1.csv",sep=""),row.names=F)

save(samelocpop1,samelocpop1halfvals,snpdistpop1,file=paste(results,"pop1samelocr2data.rdat",sep=""))


#####pop2

snpdistpop2 <- matrix(nrow=nrow(samelocpop2halfvals),ncol=ncol(samelocpop2halfvals))

for (i in 1:nrow(samelocpop2halfvals)){
	for (j in 1:ncol(samelocpop2halfvals)){
		if (!is.na(samelocpop2halfvals[i,j])){
			snpdistpop2[i,j] <- position[i] - position[j]
		}
	}
}

samelocpop2halfvalsvec <- unlist(samelocpop2halfvals)
samelocpop2vecnona <- samelocpop2halfvalsvec[!is.na(samelocpop2halfvalsvec)]
snpdistpop2vec <- unlist(snpdistpop2)
snpdistpop2vecnona <- snpdistpop2[!is.na(snpdistpop2)]

samelocpop2df <- data.frame(snpdistpop2vecnona,samelocpop2vecnona)
colnames(samelocpop2df) <- c("SNPdist","r2val")

write.csv(samelocpop2df,file=paste(results,"sameloc_pop2.csv",sep=""),row.names=F)

save(samelocpop2,samelocpop2halfvals,snpdistpop2,file=paste(results,"pop2samelocr2data.rdat",sep=""))

########pop3


snpdistpop3 <- matrix(nrow=nrow(samelocpop3halfvals),ncol=ncol(samelocpop3halfvals))

for (i in 1:nrow(samelocpop3halfvals)){
	for (j in 1:ncol(samelocpop3halfvals)){
		if (!is.na(samelocpop3halfvals[i,j])){
			snpdistpop3[i,j] <- position[i] - position[j]
		}
	}
}

samelocpop3halfvalsvec <- unlist(samelocpop3halfvals)
samelocpop3vecnona <- samelocpop3halfvalsvec[!is.na(samelocpop3halfvalsvec)]
snpdistpop3vec <- unlist(snpdistpop3)
snpdistpop3vecnona <- snpdistpop3[!is.na(snpdistpop3)]

samelocpop3df <- data.frame(snpdistpop3vecnona,samelocpop3vecnona)
colnames(samelocpop3df) <- c("SNPdist","r2val")

write.csv(samelocpop3df,file=paste(results,"sameloc_pop3.csv",sep=""),row.names=F)

save(samelocpop3,samelocpop3halfvals,snpdistpop3,file=paste(results,"pop3samelocr2data.rdat",sep=""))


########pop7

snpdistpop7 <- matrix(nrow=nrow(samelocpop7halfvals),ncol=ncol(samelocpop7halfvals))

for (i in 1:nrow(samelocpop7halfvals)){
	for (j in 1:ncol(samelocpop7halfvals)){
		if (!is.na(samelocpop7halfvals[i,j])){
			snpdistpop7[i,j] <- position[i] - position[j]
		}
	}
}

samelocpop7halfvalsvec <- unlist(samelocpop7halfvals)
samelocpop7vecnona <- samelocpop7halfvalsvec[!is.na(samelocpop7halfvalsvec)]
snpdistpop7vec <- unlist(snpdistpop7)
snpdistpop7vecnona <- snpdistpop7[!is.na(snpdistpop7)]

samelocpop7df <- data.frame(snpdistpop7vecnona,samelocpop7vecnona)
colnames(samelocpop7df) <- c("SNPdist","r2val")

write.csv(samelocpop7df,file=paste(results,"sameloc_pop7.csv",sep=""),row.names=F)

save(samelocpop7,samelocpop7halfvals,snpdistpop7,file=paste(results,"pop7samelocr2data.rdat",sep=""))

############plotting
#Note, for ggplot2 funcitons, labs() can be used to change main title and axis titles.
plot(snpdistpop1vecnona,samelocpop1vecnona,col=rgb( 27, 158, 119,50,maxColorValue=255), pch=16,ylim=c(0,1),xlab="Distance between SNPs",ylab="r^2 LD value", main="Western Population")
plot(snpdistpop1vecnona,samelocpop1vecnona,col=rgb( 117, 112, 179,50,maxColorValue=255), pch=16,ylim=c(0,1),xlab="Distance between SNPs",ylab="r^2 LD value", main="Eastern Population")
plot(snpdistpop1vecnona,samelocpop1vecnona,col=rgb( 217, 95, 2,50,maxColorValue=255), pch=16,ylim=c(0,1),xlab="Distance between SNPs",ylab="r^2 LD value", main="Hawaii Population")


ggplot(samelocpop1df,aes(x=SNPdist,y=r2val)) + stat_bin2d(colour="white") + ylim(c(0,1)) + scale_fill_gradientn(colours=c("white","#7570B3")) + labs(title="Western Population",x="SNP Distance",y="r^2 Value")
ggplot(samelocpop2df,aes(x=SNPdist,y=r2val)) + stat_bin2d(colour="white") + ylim(c(0,1)) + scale_fill_gradientn(colours=c("white","#1B9E77")) + labs(title="Eastern Population",x="SNP Distance",y="r^2 Value")
ggplot(samelocpop3df,aes(x=SNPdist,y=r2val)) + stat_bin2d(colour="white") + ylim(c(0,1)) + scale_fill_gradientn(colours=c("white","#D95F02")) + labs(title="Hawaii Population",x="SNP Distance",y="r^2 Value")


samelocpop1df$pop = rep(1,nrow(samelocpop1df))
samelocpop2df$pop = rep(2,nrow(samelocpop2df))
samelocpop3df$pop = rep(3,nrow(samelocpop3df))

samelocallpopdf <- rbind(samelocpop1df,samelocpop2df,samelocpop3df)

a <- ggplot(data=samelocpop1df,aes(x=SNPdist,y=r2val))

a + geom_point(colour='#7570B3',alpha=1/3) + stat_smooth(method= "loess",colour='#7570B3',fill='#7570B3',alpha=.3) + theme_bw()

d <- qplot(samelocpop1df$SNPdist,samelocpop1df$r2val,data=samelocpop1df)
d + stat_summary(fun.data = "mean_cl_boot",colour="#7570B3",geom="smooth") + theme_bw()


geocolors <- c("#7570B3","#1B9E77","#D95F02")
geolabs <- c("West","East","Hawaii")

all <- ggplot(data=samelocallpopdf,aes(y=r2val,x=SNPdist,colour=factor(pop)))

all + geom_point(alpha=.6) + stat_smooth(method="loess",aes(fill=factor(pop)), size=2) + scale_color_manual(name="Pop",labels=geolabs,values=geocolors) + scale_fill_manual(values=geocolors,guide="none") + theme_bw() + guides(color=guide_legend(override.aes=list(fill=geocolors))) + labs(title="All Populations loess smoothed r^2 values",x="SNP Distance",y="r^2 Value")

all + stat_smooth(method="loess",aes(fill=factor(pop)), size=2) + scale_color_manual(name="Pop",labels=geolabs,values=geocolors) + scale_fill_manual(values=geocolors,guide="none") + theme_bw() + guides(color=guide_legend(override.aes=list(fill=geocolors)))



all + geom_point(alpha=.6) + stat_summary(fun.data="mean_cl_boot",aes(fill=factor(pop)), size=1,geom="smooth") + scale_color_manual(name="Pop",labels=geolabs,values=geocolors) + scale_fill_manual(values=geocolors,guide="none") + theme_bw() + guides(color=guide_legend(override.aes=list(fill=geocolors))) + labs(title="All Populations Mean CI r^2 values",x="SNP Distance",y="r^2 Value")


all + stat_summary(fun.data="mean_cl_boot",aes(fill=factor(pop)), size=1,geom="smooth") + scale_color_manual(name="Pop",labels=geolabs,values=geocolors) + scale_fill_manual(values=geocolors,guide="none") + theme_bw() + guides(color=guide_legend(override.aes=list(fill=geocolors)))


pop1samelocmeans <- aggregate(samelocpop1df,by=list(samelocpop1df$SNPdist),mean)
pop1samelocsds <- aggregate(samelocpop1df,by=list(samelocpop1df$SNPdist),sd)
pop2samelocmeans <- aggregate(samelocpop2df,by=list(samelocpop2df$SNPdist),mean)
pop2samelocsds <- aggregate(samelocpop2df,by=list(samelocpop2df$SNPdist),sd)
pop3samelocmeans <- aggregate(samelocpop3df,by=list(samelocpop3df$SNPdist),mean)
pop3samelocsds <- aggregate(samelocpop3df,by=list(samelocpop3df$SNPdist),sd)





#############################################################################

#The next part of this script will calculated LD for all paired loci (from same fragment).  For each locus pair, an average r^2 will be calculated based on the matrix of SNP values.  Again, this will be done separately for each population.

########
pop1pairsr2 <- vector()

for (i in 1:nrow(pairs)){
	r2vec <- unlist(pop1[grep(pairs[i,1],locus),grep(pairs[i,2],locus)])
	r2vec2 <- sapply(r2vec,FUN=no0or9)
	pop1pairsr2[i] <- mean(r2vec2[!is.na(r2vec2)])
}

pop2pairsr2 <- vector()

for (i in 1:nrow(pairs)){
	r2vec <- unlist(pop2[grep(pairs[i,1],locus),grep(pairs[i,2],locus)])
	r2vec2 <- sapply(r2vec,FUN=no0or9)
	pop2pairsr2[i] <- mean(r2vec2[!is.na(r2vec2)])
}

pop3pairsr2 <- vector()

for (i in 1:nrow(pairs)){
	r2vec <- unlist(pop3[grep(pairs[i,1],locus),grep(pairs[i,2],locus)])
	r2vec2 <- sapply(r2vec,FUN=no0or9)
	pop3pairsr2[i] <- mean(r2vec2[!is.na(r2vec2)])
}

pop7pairsr2 <- vector()

for (i in 1:nrow(pairs)){
	r2vec <- unlist(pop7[grep(pairs[i,1],locus),grep(pairs[i,2],locus)])
	r2vec2 <- sapply(r2vec,FUN=no0or9)
	pop7pairsr2[i] <- mean(r2vec2[!is.na(r2vec2)])
}

allpoppairsr2 <- data.frame(pop1pairsr2,pop2pairsr2,pop3pairsr2)

boxplot(allpoppairsr2)

rand214 <- sample(1:length(uniqueloci),428)
fakepair1 <- uniqueloci[rand214[1:214]]
fakepair2 <- uniqueloci[rand214[215:428]]
fakepairs <- data.frame(fakepair1,fakepair2)

pop1fakepairsr2 <- vector()

for (i in 1:nrow(fakepairs)){
	r2vec <- unlist(pop1[grep(fakepairs[i,1],locus),grep(fakepairs[i,2],locus)])
	r2vec2 <- sapply(r2vec,FUN=no0or9)
	pop1fakepairsr2[i] <- mean(r2vec2[!is.na(r2vec2)])
}

pop2fakepairsr2 <- vector()

for (i in 1:nrow(fakepairs)){
	r2vec <- unlist(pop2[grep(fakepairs[i,1],locus),grep(fakepairs[i,2],locus)])
	r2vec2 <- sapply(r2vec,FUN=no0or9)
	pop2fakepairsr2[i] <- mean(r2vec2[!is.na(r2vec2)])
}

pop3fakepairsr2 <- vector()

for (i in 1:nrow(fakepairs)){
	r2vec <- unlist(pop3[grep(fakepairs[i,1],locus),grep(fakepairs[i,2],locus)])
	r2vec2 <- sapply(r2vec,FUN=no0or9)
	pop3fakepairsr2[i] <- mean(r2vec2[!is.na(r2vec2)])
}

pop7fakepairsr2 <- vector()

for (i in 1:nrow(fakepairs)){
	r2vec <- unlist(pop7[grep(fakepairs[i,1],locus),grep(fakepairs[i,2],locus)])
	r2vec2 <- sapply(r2vec,FUN=no0or9)
	pop7fakepairsr2[i] <- mean(r2vec2[!is.na(r2vec2)])
}

allpoppairsfakepairsr2 <- data.frame(pop1pairsr2,pop1fakepairsr2,pop2pairsr2,pop2fakepairsr2,pop3pairsr2,pop3fakepairsr2)

boxplot(allpoppairsfakepairsr2,col = rep(geocolors,each=2),ylab="r^2 (LD)")
legend("topleft",legend=c("West","East","Hawaii"),fill=geocolors)
