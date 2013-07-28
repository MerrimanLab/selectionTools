args<-commandArgs(TRUE)
#read in haps file from shapeit
hapsPop=read.table(file=args[1])
write.table(hapsPop, file=paste(args1,".mod",sep="", quote=FALSE, row.names=FALSE, col.names=FALSE))