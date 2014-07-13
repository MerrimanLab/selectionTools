#!/bin/env Rscript
# Murray Cadzow and James Boocock
# July 2013
# University of Otago
#

require(getopt)

require(rehh) 
require(parallel) #package for run in parallel in R.
#script to split chromosome into x sized segments to compute iHH on
args<-commandArgs(TRUE)
missing_code='.'
spec = matrix(c(
    'help',     'h',    0,  "logical",
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
    "haplo_hh",    "H",   0,   "logical",
    "missing_code", "M",  1,   "character",
    "physical_map_haps",  "g" , 1, "character"
), byrow=T, ncol=4) 
opt = getopt(spec)
if (!is.null(opt$help)){
    cat(getopt(spec,usage=TRUE));
    q(status=1);
}
#read in haps file from shapeit
pop1=as.character(opt$pop)
#print("*")
hapsPop=read.table(opt$input,stringsAsFactors=F,header=F)
hapsPop=hapsPop[nchar(as.character(hapsPop[,4]))==1 & nchar(as.character(hapsPop[,5]))==1, ] #remove indels
chr=as.numeric(opt$chr)
window=as.numeric(opt$window)
overlap=as.numeric(opt$overlap)
cores=as.numeric(opt$cores)
working_dir=as.character(opt$working_dir)
offset=as.numeric(opt$offset)
maf=as.numeric(opt$maf)
if(!is.null(opt$missing_code)){
    missing_code = opt$missing_code
}

if (!is.null(opt$physical_map_haps)){
    map_positions=read.table(opt$physical_map_haps, header=F)
    map_positions=as.numeric(map_positions[,1])
}else{
    map_positions=hapsPop[,3]
}
#size of each region
#window=500000
#overlap = 100000
#haps file
#hapsPop=read.table("CEU.haps")
if(!is.null(opt$haplo_hh)){
  bin = hapsPop
  hapsPop.only=bin[,6:ncol(bin)]
  allelesPop=bin[,4:5]
  hapsPop.only[hapsPop.only == "1"] = 2
  hapsPop.only[hapsPop.only == "0"] = 1
  hapsPop.only[hapsPop.only == missing_code] = 0
  t.hapsPop=t(hapsPop.only)
  ##Construct the ind file
  ind=matrix(ncol=5,nrow=nrow(bin))
  ind[,1] = as.character(bin[,2])
  ind[,2] = chr
  ind[,3] = bin[,3]
  ind[,4] = 1
  ind[,5] = 2
  indPop1=ind
  write.table(indPop1,file=paste(pop1,"chr", chr,"wd",working_dir,".map",sep="_"),col.names=F,row.names=F)
  write.table(t.hapsPop,file=paste(pop1,"chr", chr,"wd",working_dir,".haps",sep="_"),col.names=F)
  d = data2haplohh(hap_file=paste(pop1,"chr", chr,"wd",working_dir,".haps",sep="_"),map_file=paste(pop1,"chr", chr,"wd",working_dir,".map",sep="_"),min_maf=maf)   
  save(d, file=paste(pop1,"chr", chr,"wd",working_dir,".RData",sep="_"))
}
#print("Why are you not working")
#want to create overlapping bins
#column 3 is base position
setwd(working_dir)

#calculate offset from file - 
# first position in file
offset=ceiling(map_positions[1]/(window-overlap))
#pseudo code
i=ceiling(map_positions[1]/(window-overlap))

while((i-1) * (window - overlap) <= map_positions[length(map_positions)]){
  if(i == 1){
    if (!is.null(opt$physical_map_haps)){
        genetic_pos = map_positions[map_positions < (window *i)]
    }
    bin = hapsPop[map_positions < (window * i),]
  } else {
    if (!is.null(opt$physical_map_haps)){
        genetic_pos = map_positions[map_positions >  (window - overlap) * (i-1) & map_positions < (window - overlap)* i + overlap]
    }
    bin = hapsPop[ map_positions > (window - overlap) * (i-1) & map_positions < (window - overlap)* i + overlap , ]
  }
  if(length(bin[,3]) > 0){
  hapsPop.only=bin[,6:ncol(bin)]
  allelesPop=bin[,4:5]
  hapsPop.only[hapsPop.only == "1"] = 2
  hapsPop.only[hapsPop.only == "0"] = 1
  hapsPop.only[hapsPop.only == missing_code] = 0
  t.hapsPop=t(hapsPop.only)
  ##Construct the ind file
  ind=matrix(ncol=5,nrow=nrow(bin))
  ind[,1] = as.character(bin[,2])
  ind[,2] = chr
  ind[,3] = bin[,3]
  ind[,4] = 1
  ind[,5] = 2
  indPop1=ind
  #write out entire chromosome
  if (!is.null(opt$physical_map_haps)){
        write.table(genetic_pos,file=paste('gene_',pop1,'.map',i,sep=''),col.names=F,row.names=F)
  }
  write.table(indPop1,file=paste("ind_",pop1,".test",i,sep=""),col.names=F,row.names=F)
  write.table(t.hapsPop,file=paste("t_",pop1,".haps",i, sep=""),col.names=F)
  }
  i = i + 1
}
fileNumber = offset:i 
if (!is.null(opt$physical_map_haps)){
    genetic_map_file =paste0('gene_',pop1,'.map')
}
map_file=paste("ind_",pop1,".test",sep="")
hap_file=paste("t_",pop1,".haps", sep="")
#print(fileNumber)
flag = 0; 
para = list(); 
new_file_number = 0
for( i in fileNumber){  
  if(file.exists(paste(hap_file,i,sep=""))){ 
    if (!is.null(opt$physical_map_haps)){
        p = c(paste(hap_file,i,sep=""), paste(map_file,i,sep=""),paste(genetic_map_file,i,sep=''))   
    }else{
        p = c(paste(hap_file,i,sep=""), paste(map_file,i,sep=""),-1)   
    }
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
  if(!is.null(opt$physical_map_haps)){
    physical_positions = read.table(x[3],header=F)
    physical_positions = as.numeric(physical_positions[,1])
    res = scan_hh(d,big_gap=opt$big_gap,small_gap=opt$small_gap,small_gap_penalty=opt$small_gap_penalty,physical_positions=physical_positions)
  }else{
    res = scan_hh(d,big_gap=opt$big_gap,small_gap=opt$small_gap,small_gap_penalty=opt$small_gap_penalty)
  }
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
save.image(file="working_data.RData")

results=data.frame()
if(!is.null(opt$physical_map_haps)){
for (n in fileNumber){
  i=n-(offset-1)
  temp_physical_map= as.numeric(read.table(paste0('gene_',pop1,'.map',n),header=F)[,1])
  if(n == 1){ # from start to first half of overlaped region (first chunk)
    results = neutral_res[[i]][temp_physical_map[i] <= ((n+offset-1) * window - 1/2 *overlap) ,] #correct window
   } else {
      if(n == max(fileNumber)){ #take second half of overlap at start and go until the end (final chunk)
        a= results
        b = neutral_res[[i]][ ((window-overlap)* (n-1) + 1/2*overlap) <= temp_physical_map[i]  ,]
        results = rbind(a,b)
      } else { #start =take second half of overlap, end = take first half (middle regions)
        a = results
        b = neutral_res[[i]][ ((window-overlap)* (n-1) + 1/2*overlap) <= temp_physical_map[i]  & temp_physical_map[i] <  ((window -overlap)* (n) + (1/2 * overlap)), ]
        results = rbind(a,b )
     }
   } 
}
if (!is.null(opt$physical_map_haps)){
    # Need to replace the second column of 
    for( i in 1:nrow(results)){
        results[i,2] = map_positions[i]     
    }
}
write.table(results,paste(pop1,"chr", chr,"wd",working_dir,".ihh",sep="_"))
if (!is.null(opt$ihs)){
        ihs =ihh2ihs(results)
        write.table(ihs$res.ihs,paste(pop1,"chr", chr,"wd",working_dir,".ihs",sep="_"))
    }
}else{
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

}

