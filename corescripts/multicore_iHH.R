#!/bin/env Rscript
# Murray Cadzow and James Boocock
# July 2013 updated March 2015
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
    'input',   'i', 1,  "character", #haps file
    'chr',     'c', 1,  "integer",
    'window',  'w', 1,  "integer", #window size in bp
    'overlap', 'o', 1,  "integer", #overlap between windows in bp
    'cores',   'r', 1,  "integer", # cores to use in parallel
    'working_dir', 'd', 1, "character",
    'offset',   's' , 1, "integer",
    "maf"  ,    'm' , 1, "integer", # allele frequency threshold below which to exclude snp from calculations
    "pop"  ,    'p', 1, "character",
    "ihs"  ,    'I', 0, "logical", #calculate ihs
    "big_gap",  'b', 1, "integer", # size of gap in bp to stop calculating ihh (200,000 suggested, voight et al 2006)
    "small_gap", 'S', 1, "integer", #size of gap to start applying gap penalty (20,000 suggested)
    "small_gap_penalty", 'P', 1, "integer", # gap penalty to apply (20,000 suggested)
    "haplo_hh",    "H",   0,   "logical", #do bifurcation diagram pre-calculation
    "missing_code", "M",  1,   "character",
    "physical_map_haps",  "g" , 1, "character" #file with physical positions if haps file has genetic positions
), byrow=T, ncol=4) 
opt = getopt(spec)
if (!is.null(opt$help)){
    cat(getopt(spec,usage=TRUE));
    q(status=1);
}
#default sizes of each region that should be used
#window=10000000
#overlap = 2000000

#read in haps file from shapeit
pop1=as.character(opt$pop)
hapsPop=read.table(opt$input,stringsAsFactors=F,header=F) #haplotype file, ideally coded ancestral and derived
hapsPop=hapsPop[nchar(as.character(hapsPop[,4]))==1 & nchar(as.character(hapsPop[,5]))==1, ] #remove indels
chr=as.numeric(opt$chr)
window=as.numeric(opt$window)
overlap=as.numeric(opt$overlap)
cores=as.numeric(opt$cores)
working_dir=as.character(opt$working_dir)
offset=as.numeric(opt$offset) #how many windows to offset the start by. default should be 1
maf=as.numeric(opt$maf) #filter of derived allele freq (MAF) for ihs calculation
if(!is.null(opt$missing_code)){
    missing_code = opt$missing_code
}

if (!is.null(opt$physical_map_haps)){ #physical_map_haps option is for use after interpolation script of pipeline
    map_positions=read.table(opt$physical_map_haps, header=F)
    map_positions=as.numeric(map_positions[,2])
}else{
    map_positions=hapsPop[,3]
}

#haplo_hh option was so bification diagrams could be used
# uncomment d= ... to create haplo_hh object
# not practical for genomewide data
# not sure if this chunk should remain but still needed for next chunk to work
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
# little hack to rename duplicated rsIDs that may pop up due to imputation
  indPop1[duplicated(indPop1[,1]),1] = paste(indPop1[duplicated(indPop1[,1]),2],indPop1[duplicated(indPop1[,1]),3],sep=":")

  write.table(indPop1,file=paste(pop1,"chr", chr,"wd",working_dir,".map",sep="_"),col.names=F,row.names=F)
  write.table(t.hapsPop,file=paste(pop1,"chr", chr,"wd",working_dir,".haps",sep="_"),col.names=F)
  #d = data2haplohh(hap_file=paste(pop1,"chr", chr,"wd",working_dir,".haps",sep="_"),map_file=paste(pop1,"chr", chr,"wd",working_dir,".map",sep="_"),min_maf=maf)   
  #save(d, file=paste(pop1,"chr", chr,"wd",working_dir,".RData",sep="_"))
	rm(ind,indPop1,bin,t.hapsPop, hapsPop.only)
}


#column 3 is base position
setwd(working_dir)

###
### ACTUALLY STARTED DOING STUFF FROM HERE:
###
## goal is to break chromosome into overlapping windows
## and process each window in parallel then join them back together
## at the end

# first window number in file
offset=ceiling(map_positions[1]/(window-overlap))

# number of windows needed
i=ceiling(map_positions[1]/(window-overlap))

# create overlapping windows
# for each window need to create both a haps file and a map file
# physical map = pos in bp
# genetic map is in centimorgans
# i used as window index

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
  indPop1[duplicated(indPop1[,1]),1] = paste(indPop1[duplicated(indPop1[,1]),2],indPop1[duplicated(indPop1[,1]),3],sep=":")
  #write out entire chromosome
  if (!is.null(opt$physical_map_haps)){
        write.table(genetic_pos,file=paste('gene_',pop1,'.map',i,sep=''),col.names=F,row.names=F)
  }
  write.table(indPop1,file=paste("ind_",pop1,".test",i,sep=""),col.names=F,row.names=F)
  write.table(t.hapsPop,file=paste("t_",pop1,".haps",i, sep=""),col.names=F)
  }
  i = i + 1
}
rm(ind,indPop1,bin,t.hapsPop, hapsPop.only)

#the indices of the window files created
fileNumber = offset:i 
if (!is.null(opt$physical_map_haps)){
    genetic_map_file =paste0('gene_',pop1,'.map')
}

map_file=paste("ind_",pop1,".test",sep="")
hap_file=paste("t_",pop1,".haps", sep="")
#print(fileNumber)

# set up a list of all the files that need to have ihh calculated
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
actualfileNumber = offset:(offset+new_file_number-1)

# function to perform the ihh calculation
# ideally this will be updated at some stage
my_scan_hh = function(x){     
  if(length(read.table(x[1]) != 0)){
    d = data2haplohh(hap_file=x[1],map_file=x[2],min_maf=maf)   
    if(!is.null(opt$physical_map_haps)){
        physical_positions = read.table(x[3],header=F)
        physical_positions = as.numeric(physical_positions[,1])
        res = scan_hh(d,big_gap=opt$big_gap,small_gap=opt$small_gap,small_gap_penalty=opt$small_gap_penalty,physical_positions=physical_positions)
    }else{
        res = scan_hh(d,big_gap=opt$big_gap,small_gap=opt$small_gap,small_gap_penalty=opt$small_gap_penalty)
    }
    write.table(res,paste(x[1],".iHH",sep=""))    
    
    } 
    return(NULL)
}
 
# calculate ihh in parallel 
neutral_res = mclapply(para,my_scan_hh,mc.cores=cores)  

# re-read in all the results files - this is to avoid index problems due to missing windows
index = 0
neutral_res=list()
for ( j in fileNumber){
index = index + 1
    if(file.exists(paste(hap_file,j,'.iHH',sep=''))){
        neutral_res[[j]] = read.table(paste(hap_file,j,'.iHH',sep=''))
    } else{
        print(paste(hap_file,j,'.iHH does not exist! Continuing without',sep=""))
    }    
}
#save(neutral_res,file="neutral_res.RData")
save.image(file="working_data.RData")

#combine all the windows into a single chromosome again
results=data.frame()
if(!is.null(opt$physical_map_haps)){
  print(fileNumber)
  for (n in fileNumber){
    i=n-(offset-1) #keeps track of if window is the first
    if(file.exists(paste0('gene_',pop1,'.map',n))){
      temp_physical_map= as.numeric(read.table(paste0('gene_',pop1,'.map',n),header=F)[,1])
      if(i == 1){ # from start to first half of overlaped region (first chunk)
        results = neutral_res[[n]][temp_physical_map <= ((n+offset-1) * window - 1/2 *overlap) ,] #take correct window when window1 != file1             
      } else {
        if(n == max(fileNumber)){ #take second half of overlap at start and go until the end (final chunk)
          results = rbind(results, neutral_res[[n]][ ((window-overlap)* (n-1) + 1/2*overlap) <= temp_physical_map  ,])                  
        } else { #start =take second half of overlap, end = take first half (middle regions)
               results = rbind(results,neutral_res[[n]][ ((window-overlap)* (n-1) + 1/2*overlap) <= temp_physical_map  & temp_physical_map <  ((window -overlap)* (n) + (1/2 * overlap)), ])      
        }
      } 
    }else{
      print(paste("File: gene_",pop1,'.map',n," DOES NOT EXIST, skipping region.", sep=""))
    }
  }
  #names(map)= c("name", "name2", "gen_pos", "a1", "a2","phys_pos")
  results$name = rownames(results[,])
  ##### replace genetic positions with physical positions
  if (!is.null(opt$physical_map_haps)){
    #results[,2] = map_positions
    m= read.table(opt$physical_map_haps, header=F)
    m[,3]=chr
    m[duplicated(m[,1]),1] = paste(m[duplicated(m[,1]),3],m[duplicated(m[,1]),2],sep=":")
    results$name = rownames(results[,])
    z = merge(results, m, by.x="name", by.y="V1")
    results = z[,c("CHR","V2", "FREQ_a","IHHa","IHHd", "IES")]
    names(results) = c("CHR","POSITION", "FREQ_a","IHHa","IHHd", "IES")
    rownames(results) = z$name    
    rm(z)  
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
                results = rbind(results, neutral_res[[i]][ ((window-overlap)* (n-1) + 1/2*overlap) <= neutral_res[[i]][,2]  ,])
            } else { #start =take second half of overlap, end = take first half (middle regions)
                results = rbind(results, neutral_res[[i]][ ((window-overlap)* (n-1) + 1/2*overlap) <= neutral_res[[i]][,2]  & neutral_res[[i]][,2] <  ((window -overlap)* (n) + (1/2 * overlap)),])
            }
        } 
    }
    write.table(results,paste(pop1,"chr", chr,"wd",working_dir,".ihh",sep="_"))
    if (!is.null(opt$ihs)){
        ihs =ihh2ihs(results)
        write.table(ihs$res.ihs,paste(pop1,"chr", chr,"wd",working_dir,".ihs",sep="_"))
    }   
}

