basicphylo <- read.table("~/Dropbox/HFRad-Tags/HF_m4M3n5_Results/PopMapBasicPhyloNoOutgroups.txt")

pops <- c(1,2,3)

indpops <- list()
numinds <- vector()
for (i in 1:length(pops)){
	indpops[[i]] <- as.character(basicphylo[basicphylo[,"V2"]==i,1])
	numinds[i] <- (length(indpops[[i]]))
	}
	
subsamplesize <- 8


#Use script below to do once
newindpops <- list()
newpopdes <- list()
for (i in 1:length(indpops)){
	newindpops[[i]] <- sample(indpops[[i]],subsamplesize)
	newpopdes[[i]] <- rep(pops[i],subsamplesize)
}

indvec <- unlist(newindpops)
popvec <- unlist(newpopdes)

newpopmap <- data.frame(indvec,popvec)
write.table(newpopmap,file="~/Dropbox/HFRad-Tags/HFPaired_Reduced/PopMapBasicPhyloSameNumInds.txt",quote=F,sep="\t",row.names=F,col.names=F)




#To do the same as above but yet 10 different map replicates:
for (j in 1:10){
	newindpops <- list()
	newpopdes <- list()
	for (i in 1:length(indpops)){
		newindpops[[i]] <- sample(indpops[[i]],subsamplesize)
		newpopdes[[i]] <- rep(pops[i],subsamplesize)
	}
	
	indvec <- unlist(newindpops)
	popvec <- unlist(newpopdes)
	
	newpopmap <- data.frame(indvec,popvec)
	write.table(newpopmap,file=paste("~/Dropbox/HFRad-Tags/HF_m4M3n5_Results/ReplicatePopmaps_8indivs/PopMapBasicPhyloSameNumInds_",j,".txt",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
}

