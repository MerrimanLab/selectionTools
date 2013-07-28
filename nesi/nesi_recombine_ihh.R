args<-commandArgs(TRUE)
pop=as.character(args[1])
chr=as.numeric(args[2])
window=as.numeric(args[3]) # size of the small windows
smallWindowOverlap=as.numeric(args[4])) # overlap for the small windows
overlap=as.numeric(args[5]) # overlap between the large windows
cores=as.numeric(args[6])



##use for testing
#window=5000000
#overlap = 1000000
#cores = 10
bigWindow= (window-smallWindowOverlap) * (cores-1) + window

setwd(workingdir)
fileList=dir(pattern="*wd.*ihh", recursive=TRUE)
ihh_list=list()
for(i in 1:length(fileList)){ihh_list[[i]]<-read.table(fileList[i])}
names(ihh_list)<-fileList
fileNumber=1:length(ihh_list)

for (n in fileNumber){
  print((n -1)* (bigWindow-overlap))
  print(paste("window",n,": is from:",( (bigWindow-overlap) * (n-1)), "to:", ((bigWindow - overlap)* n + overlap) , sep=" "))
  print(paste("merge window",n,": is from:",((bigWindow-overlap)* (n-1) + 1/2*overlap), "to:", ((bigWindow -overlap)* n + (1/2 * overlap)), sep=" "))
  
  if(n == 1){ # from start to first half of overlaped region (first chunk)
    print("n=1")
    results = ihh_list[[n]][ihh_list[[n]][,2] <= (n * bigWindow - 1/2 *overlap) ,] #correct window
    #print(max(results[,2]))
  } else {
    if(n == max(fileNumber)){ #take second half of overlap at start and go until the end (final chunk)
      #print("max")
      a= results
      b = ihh_list[[n]][ ((bigWindow-overlap)* (n-1) + 1/2*overlap) < ihh_list[[n]][,2]  ,]
      print(max(results[,2]))
      print(min(b[,2]))
      results = rbind(a,b)
      print(max(results[,2]))
    } else { #start =take second half of overlap, end = take first half (middle regions)
      #print("middle")
      
      a = results
      b = ihh_list[[n]][ ((bigWindow-overlap)* (n-1) + 1/2*overlap) < ihh_list[[n]][,2]  & ihh_list[[n]][,2] <=  ((bigWindow -overlap)* n + (1/2 * overlap)), ]
      #print(max(a[,2]))
      #print(min(b[,2]))
      results = rbind(a,b )
      #print(max(results[,2]))
    }
  } 
}

write.table(results, file=paste(pop,"_chr",chr,".iHH",sep=""),quote=FALSE, sep="\t")