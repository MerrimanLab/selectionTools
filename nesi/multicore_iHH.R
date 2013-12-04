#
# Murray Cadzow and James Boocock
# July 2013
# University of Otago
#

require(getopt)

require(rehh) 
require(multicore) #package for run in parallel in R.
#script to split chromosome into x sized segments to compute iHH on
args<-commandArgs(TRUE)
spec = matrix(c(
	'help', 	'h', 	0, 	"logical",
	'input',   'i', 1,  "character",
	'chr',     'c', 1,  "integer",
	'window',  'w', 1,  "integer",
	'overlap', 'o', 1,  "integer",
	'cores',   'r', 1,  "integer",
	'working_dir', 'd', 1, "character",
	'offset',   's' , 1, "integer",
	"maf"  ,    'm' , 1, "integer",
	"pop"  ,    'p', 1, "character",
	"ihs"  ,    'I', 0, "logical",
	"big_gap",  'b', 1, "integer",
	"small_gap", 'S', 1, "integer",
	"small_gap_penalty", 'P', 1, "integer",
    "haplo-hh",    "H",   1,   "character"
), byrow=T, ncol=4) 
opt = getopt(spec)


if (!is.null(opt$help)){
	cat(getopt(spec,usage=TRUE));
	q(status=1);
}
#read in haps file from shapeit
pop1=as.character(opt$pop)
#print("*")
hapsPop=read.table(opt$input)
hapsPop=hapsPop[nchar(as.character(hapsPop[,4]))==1 & nchar(as.character(hapsPop[,5]))==1, ] #remove indels
chr=as.numeric(opt$chr)
window=as.numeric(opt$window)
overlap=as.numeric(opt$overlap)
cores=as.numeric(opt$cores)
working_dir=as.character(opt$working_dir)
offset=as.numeric(opt$offset)
maf=as.numeric(opt$maf)
#size of each region
#window=500000
#overlap = 100000
#haps file
#hapsPop=read.table("CEU.haps")


#print("Why are you not working")
#want to create overlapping bins
#column 3 is base position
setwd(working_dir)

#calculate offset from file - 

# first position in file
offset=ceiling(hapsPop[1,3]/(window-overlap))
#pseudo code
i=ceiling(hapsPop[1,3]/(window-overlap))
while((i-1) * (window - overlap) <= hapsPop[length(hapsPop[,3]),3]){
  if(i == 1){
    bin = hapsPop[hapsPop[,3] < (window * i),]
  } else {
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
  i = i + 1
  
}
fileNumber = offset:i 
map_file=paste("ind_",pop1,".test",sep="")
hap_file=paste("t_",pop1,".haps", sep="")


#print(fileNumber)
flag = 0; 
para = list(); 
new_file_number = 0
for( i in fileNumber){  
  if(file.exists(paste(hap_file,i,sep=""))){ 
    p = c(paste(hap_file,i,sep=""), paste(map_file,i,sep=""))   
    new_file_number = new_file_number + 1 
    if(flag==0){        
        para = list(p)      
    }else{        
        para = c(para,list(p))     
    }     
    flag = 1;  
    }
}  
fileNumber = offset:(offset+new_file_number-1)

my_scan_hh = function(x){     
  d = data2haplohh(hap_file=x[1],map_file=x[2],min_maf=maf)    
  res = scan_hh(d,big_gap=opt$big_gap,small_gap=opt$small_gap,small_gap_penalty=opt$small_gap_penalty)
  write.table(res,paste(x[1],".iHH",sep=""))
  return(res)
}  

neutral_res = mclapply(para,my_scan_hh,mc.cores=cores)  

index = 1
for ( j in fileNumber){
	neutral_res[[index]] = read.table(paste(hap_file,j,'.iHH',sep=''))
    index = index + 1
}
#save(neutral_res,file="neutral_res.RData")
#save.image(file="working_data.RData")

results=data.frame()
for (n in fileNumber){
  i=n-(offset-1)
  if(n == 1){ # from start to first half of overlaped region (first chunk)
    results = neutral_res[[i]][neutral_res[[i]][,2] <= ((n+offset-1) * window - 1/2 *overlap) ,] #correct window
   } else {
      if(n == max(fileNumber)){ #take second half of overlap at start and go until the end (final chunk)
        a= results
        b = neutral_res[[i]][ ((window-overlap)* (n-1) + 1/2*overlap) <= neutral_res[[i]][,2]  ,]
        results = rbind(a,b)
      } else { #start =take second half of overlap, end = take first half (middle regions)
        a = results
        b = neutral_res[[i]][ ((window-overlap)* (n-1) + 1/2*overlap) <= neutral_res[[i]][,2]  & neutral_res[[i]][,2] <  ((window -overlap)* (n) + (1/2 * overlap)), ]
        results = rbind(a,b )
     }
   } 
}
write.table(results,paste(pop1,"chr", chr,"wd",working_dir,".ihh",sep="_"))
if (!is.null(opt$ihs)){
		ihs =ihh2ihs(results)
		write.table(ihs$res.ihs,paste(pop1,"chr", chr,"wd",working_dir,".ihs",sep="_"))
}
