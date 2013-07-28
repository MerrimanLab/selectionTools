#script to split chromosome into x sized segments to compute iHH on
args<-commandArgs(TRUE)
#read in haps file from shapeit
pop1=as.character(args[1])
hapsPop=read.table(file=args[2])
chr=as.numeric(args[3])
window=as.numeric(args[4])
overlap=as.numeric(args[5])
cores=as.numeric(args[6])
working_dir=as.character(args[7])
offset=as.numeric(args[8])
#size of each region
#window=500000
#overlap = 100000
#haps file
#hapsPop=read.table("CEU.haps")

#want to create overlapping bins
#column 3 is base position
setwd(working_dir)
#pseudo code
i=0
while(i * (window - overlap) <= hapsPop[length(hapsPop[,3]),3]){
  i = i + 1
  if(i == 1){
    bin = hapsPop[hapsPop[,3] < (window * i),]
  } else {
    #next window is "shifted" back by the overlap distance
    # starts at the old max (i-1 windows) - overlap, ends at new max (i windows) - overlap
    bin = hapsPop[ hapsPop[,3] > (window - overlap) * (i-1) & hapsPop[,3] < (window - overlap)* i + overlap , ]
  }
  if(length(bin[,3]) > 0){
  hapsPop.only=bin[,6:ncol(bin)]
  allelesPop=bin[,4:5]
  hapsPop.only[hapsPop.only == 1] = 2
  hapsPop.only[hapsPop.only == 0] = 1
  t.hapsPop=t(hapsPop.only)
  
  ##Construct the ind file
  ind=matrix(ncol=5,nrow=nrow(bin))
  ind[,1] = as.character(bin[,1])
  ind[,2] = chr
  ind[,3] = bin[,3]
  ind[,4] = 1
  ind[,5] = 2
  ind=subset(ind,!duplicated(bin[,1])&&!duplicated(bin[,3]))
  ind=subset(ind,!duplicated(ind[,1]))
  ind=subset(ind,!duplicated(ind[,3]))
  indPop1=ind
  #write out entire chromosome
  write.table(indPop1,file=paste("ind_",pop1,".test",i,sep=""),col.names=F,row.names=F)
  write.table(t.hapsPop,file=paste("t_",pop1,".haps",i, sep=""),col.names=F)
  }
  
}
#column 3 of haps file is position



#d_pop1 = data2haplohh(map_file=paste("ind_",pop1,".test", sep=""),hap_file=paste("t_",pop1,".haps",sep=""))
#result_ihh_pop1 = scan_hh(d_pop1)
#result_pop1 = ihh2ihs(result_ihh_pop1)



## works below here ##

# Using rehh package to compute iHS  
library(rehh) 
library(multicore) #package for run in parallel in R.


#hap_file="neutral_data_rehh/hap_neutral_"; 
#map_file="neutral_data_rehh/map_neutral_"; 


fileNumber = offset:i 
map_file=paste("ind_",pop1,".test",sep="")
hap_file=paste("t_",pop1,".haps", sep="")

flag = 0; 
para = list(); 
for( i in fileNumber){     
  p = c(paste(hap_file,i,sep=""), paste(map_file,i,sep=""))   
  
  if(flag==0){        
    para = list(p)      
  }else{        
    para = c(para,list(p))     
  }     
  flag = 1;  
}  

my_scan_hh = function(x){     
  d = data2haplohh(hap_file=x[1],map_file=x[2])     
  res = scan_hh(d)
  write.table(res,paste(x[1],".iHH",sep=""))
}  

# run in parallel, using the number of cores specified by the arguments. 
neutral_res = mclapply(para,my_scan_hh,mc.cores=cores)  

for ( j in fileNumber){
	neutral_res[[j]] = read.table(paste(hap_file,i,'.iHH',sep=''))
}
save(neutral_res,file="neutral_res.RData")


#bin regions
#i=0
#w = window
#while(i * (window-overlap) <= hapsPop[length(hapsPop[,3]),3]){
#  i = i + 1
#  print(paste("window:",i, sep=" "))
#  if(i == 1){
#    print(window * i)
#  } else {
#    #next window is "shifted" back by the overlap distance
#    print((window - overlap) * (i-1))
#    print((window -overlap) * i + overlap)
#    print("")
#  }
#  
#}



#combine iHH results from window

for (n in seq(fileNumber)){
  print((n -1)* (window-overlap))
  print(paste("window",n,": is from:",( (window-overlap) * (n-1)), "to:", ((window - overlap)* n + overlap) , sep=" "))
  print(paste("merge window",n,": is from:",((window-overlap)* (n-1) + 1/2*overlap), "to:", ((window -overlap)* n + (1/2 * overlap)), sep=" "))
  
  if(n == 1){ # from start to first half of overlaped region (first chunk)
  print("n=1")
    results = neutral_res[[n]][neutral_res[[n]][,2] <= (n * window - 1/2 *overlap) ,] #correct window
    print(max(results[,2]))
   } else {
      if(n == max(fileNumber)){ #take second half of overlap at start and go until the end (final chunk)
        print("max")
        a= results
        b = neutral_res[[n]][ ((window-overlap)* (n-1) + 1/2*overlap) < neutral_res[[n]][,2]  ,]
        print(max(results[,2]))
        print(min(b[,2]))
        results = rbind(a,b)
        print(max(results[,2]))
      } else { #start =take second half of overlap, end = take first half (middle regions)
        print("middle")
        
        a = results
        b = neutral_res[[n]][ ((window-overlap)* (n-1) + 1/2*overlap) < neutral_res[[n]][,2]  & neutral_res[[n]][,2] <=  ((window -overlap)* n + (1/2 * overlap)), ]
        print(max(a[,2]))
        print(min(b[,2]))
        results = rbind(a,b )
        print(max(results[,2]))
     }
   } 
}
save.image("multi_core_rehh.RData") 
