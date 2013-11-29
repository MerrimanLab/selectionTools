#
#  Murray Cadzow and James Boocock
#  November 2013
#  University of Otago
#
#

require(getopt)

require(rehh)

args=commandArgs(TRUE)
spec = matrix(c(
	'pop1', 	'p', 1, 'character',
	'pop2',   'P', 1, 'character',
	'chr',    'c', 1, 'character',
	'pop1file','-i',1 ,'character',
	'pop2file,','-I',1,'character'
),byrow=T,ncol=4)
opt = getopt(spec)	

if(!is.null(opt$help)){
	cat(getopt(spec,usage=TRUE));
	q(status=1);
}
pop1_data = read.table(opt$pop1file,header=T,row.names=1)
pop2_data = read.table(opt$pop2file,header=T,row.names=1)

rsb_out = ies2rsb(pop1_data,pop2_data,popname1=opt$pop1,popname2=opt$pop2)
write.table(rsb_out$res.rsb,paste(opt$chr,opt$pop1,opt$pop2,'.rsb',sep=""))
