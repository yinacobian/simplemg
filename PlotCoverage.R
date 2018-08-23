#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


genomename <- args[1]
samplename <- args[2]
path <- args[3] 

covfile <- paste (path, "/","cov_mpileup_map_", genomename, "_", samplename, ".tab", sep="")
cov <- read.table (covfile)
avgcov <- sum(cov$V2)/(length(cov$V2))
okavgcov <-(format (avgcov, digits=2))
xlabname <- paste (genomename, "     Average coverage ", okavgcov, "X", sep="")
imagename <- paste (path, "/", "covplot_", samplename, "_", genomename,".png", sep="")
png(imagename,width=25,height=10,units="cm",res=1200)
plot (cov$V2, type= "h", col="darkblue", xlab=xlabname, ylab="Coverage/nt", main=samplename)
dev.off()


