#
# Murray Cadzow
# July 2013
# University of Otago

args<-commandArgs(TRUE)
#read in haps file from shapeit
hapsPop=read.table(file=args[1])
maf=as.numeric(args[2])
output=as.character(args[3])

hapsPop=read.table(file=args[1])
hapsPop=hapsPop[nchar(as.character(hapsPop[,4]))==1 &
                     as.character(hapsPop[,4]) != "-" &
                     nchar(as.character(hapsPop[,5]))==1 &
                     as.character(hapsPop[,5]) != "-", ] #remove indels
hapsPop[,1]= hapsPop[,2]


af = apply(hapsPop[,6:length(hapsPop[1,])], 1, FUN=function(x){sum(x)/length(x)}  )
hapsPop=hapsPop[(af > maf) & (af < (1-maf)),]
print(paste("af >",maf))
print(table((af > maf ) & (af < (1-maf))))

write.table(hapsPop, file=output, quote=FALSE, row.names=FALSE, col.names=FALSE)
