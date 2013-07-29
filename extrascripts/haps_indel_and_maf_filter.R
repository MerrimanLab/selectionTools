#
# Murray Cadzow
# July 2013
# University of Otago

args<-commandArgs(TRUE)
#read in haps file from shapeit
hapsPop=read.table(file=args[1])
maf=as.numeric(args[2])
af = apply(hapsPop[,6:length(hapsPop[1,])], 1, FUN=function(x){sum(x)/length(x)}  )
hapsPop=hapsPop[(af > maf) & (af < (1-maf)),]
print(paste("af >",maf))
print(table((af > maf ) & (af < (1-maf))))

write.table(hapsPop, file=paste(args[1],".mod",sep="", quote=FALSE, row.names=FALSE, col.names=FALSE))
