#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


genomeid <- args[1]
samplename <- args[2]
path <- args[3] 

covfile <- paste (path, "/","cov_mpileup_map_", genomeid, "_", samplename, ".tab", sep="")
coverage <- read.table (covfile)
namesfile <- paste (path, "/","ok_id_names_topspecies_", samplename, ".txt", sep="")
namesfile <- read.table (namesfile, sep='\t', row.names=1)
genomename <- namesfile [genomeid, ]
avgcov <- (sum(coverage$V2))/(length(coverage$V2))
okavgcov <-(format (avgcov, digits=2))

tb <- table((coverage$V2))
pnocoverage <- (((tb[1])*100)/(length(coverage$V2)))
okpnocoverage <-(format (pnocoverage, digits=2))

xlabname <- paste (genomename, " ", "(",genomeid,")", sep="")
labcoverage <- paste ("Average coverage ", okavgcov, "X", ",", " regions with no coverage ", okpnocoverage, "%", sep="")

#imagename <- paste (path, "/", "log","covplot_", samplename, "_", genomeid,".png", sep="")
imagename <- paste (path, "/","covplot_", samplename, "_", genomeid,".png", sep="")

png(imagename,width=25,height=10,units="cm",res=1200)
#plot (log(coverage$V2), type= "h", col="darkblue", xlab=xlabname, ylab="Coverage/nt", main=samplename, sub=labcoverage)
plot (coverage$V2, type= "h", col="darkblue", xlab=xlabname, ylab="Coverage/nt", main=samplename, sub=labcoverage)


dev.off()



