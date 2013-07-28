#Set directory to Rehh directory
#setwd('~/R/rehh')
args<-commandArgs(TRUE)

#read in haps file from shapeit
pop1=as.character(args[1])
hapsPop1=read.table(file=args[2])
pop2=as.character(args[3])
hapsPop2=read.table(file=args[4])
chr=as.numeric(args[5])

dir.create(pop1)
dir.create(pop2)



#haps=read.table('../PPARGC1A/PPARGC1A_merriman_sequences_phased.haps')
#samples=read.table('pparg.haps',header=T)
#haps=subset(haps,!duplicated(haps[,1]))
#haps=subset(haps,!duplicated(haps[,3]))
setwd(pop1)
hapsPop1.only=hapsPop1[,6:ncol(hapsPop1)]
allelesPop1=hapsPop1[,4:5]
hapsPop1.only[hapsPop1.only == 1] = 2
hapsPop1.only[hapsPop1.only == 0] = 1
t.hapsPop1=t(hapsPop1.only)
library(rehh)
## change some rehh graphs functions to plot
ihsplotMod = function (data, plot.pval = "TRUE", ylim.scan = 2, pch = 16, 
                       main = "iHS") 
{
  tmp_chr = unique(data[, 1])
  col_chr = 1:length(tmp_chr)
  names(col_chr) = tmp_chr
  pos_chr = rep(0, length(tmp_chr))
  tmp_nmrk = table(data[, 1])
  pos_mrk = cumsum(tmp_nmrk)
  pos_chr[1] <- floor(pos_mrk[1]/2)
  if (length(tmp_chr) > 1) {
    for (i in 2:length(tmp_chr)) {
      pos_chr[i] = pos_mrk[i - 1] + floor((tmp_nmrk[i]/2))
    }
  }
  plot(data[,2],data[, 3], pch = pch, las = 1, col = col_chr[as.character(data[, 
                                                                      1])],# xaxt = "n"
       , xlab = paste("Chromosome", chr,sep=" "), ylab = "iHS", 
       main = main)
  abline(h = ylim.scan, lty = 2)
  abline(h = -1 * ylim.scan, lty = 2)
  #axis(1, at = pos_chr, labels = tmp_chr, las = 1)
  if (plot.pval) {
    plot(data[, 4], pch = pch, las = 1, col = col_chr[as.character(data[, 
                                                                        1])], xaxt = "n", xlab = "Chromosome", main = "Pvalue", 
         ylab = expression("-" * log[10] * "[" ~ "1-2|" * 
                             Phi[scriptstyle(italic(iHS))] * "-0.5|" ~ "]"))
    abline(h = ylim.scan, lty = 2)
    abline(h = -1 * ylim.scan, lty = 2)
    axis(1, at = pos_chr, labels = tmp_chr, las = 1)
    abline(h = ylim.scan, lty = 2)
    abline(h = -1 * ylim.scan, lty = 2)
  }
}




##Construct the ind file
ind=matrix(ncol=5,nrow=nrow(hapsPop1))
ind[,1] = as.character(hapsPop1[,1])
ind[,2] = chr
ind[,3] = hapsPop1[,3]
ind[,4] = 1
ind[,5] = 2
ind=subset(ind,!duplicated(hapsPop1[,1])&&!duplicated(hapsPop1[,3]))
ind=subset(ind,!duplicated(ind[,1]))
ind=subset(ind,!duplicated(ind[,3]))
indPop1=ind
write.table(indPop1,file=paste("ind_",pop1,".test",sep=""),col.names=F,row.names=F)
write.table(t.hapsPop1,file=paste("t_",pop1,".haps", sep=""),col.names=F)
d_pop1 = data2haplohh(map_file=paste("ind_",pop1,".test", sep=""),hap_file=paste("t_",pop1,".haps",sep=""))
result_ihh_pop1 = scan_hh(d_pop1)
result_pop1 = ihh2ihs(result_ihh_pop1)
png(file=paste("ihsplot",pop1,".png",sep=""))
ihsplot(result_pop1$res.ihs,plot.pval=TRUE,ylim.scan=2,main="iHS")
dev.off()
write.table(result_pop1$res.ihs,file=paste("res.ihs.",pop1,".txt",sep=""), sep="\t")
save.image(file="ihs.RData")

sig=rownames(result_pop1$res.ihs[result_pop1$res.ihs[,4] > 2 & result_pop1$res.ihs[result_pop1$res.ihs[,4] < -2],])
sig=subset(sig, sig != "NA")
sigPop1=sig

#make bification plot for specified marker
i= 1
for (marker in sigPop1){
  png(file=paste(marker,pop1,".png", sep=""))
bifurcation.diagram(d_pop1,mrk_foc=which(hapsPop1[,1]==marker))
dev.off()
}


#################
#################
###### pop2 #####

setwd(paste("../",pop2, sep=""))

hapsPop2.only=hapsPop2[,6:ncol(hapsPop2)]
allelesPop2=hapsPop2[,4:5]
hapsPop2.only[hapsPop2.only == 1] = 2
hapsPop2.only[hapsPop2.only == 0] = 1
t.hapsPop2=t(hapsPop2.only)


##Construct the ind file
ind=matrix(ncol=5,nrow=nrow(hapsPop2))
ind[,1] = as.character(hapsPop2[,1])
ind[,2] = 4
ind[,3] = hapsPop2[,3]
ind[,4] = 1
ind[,5] = 2
ind=subset(ind,!duplicated(hapsPop2[,1])&&!duplicated(hapsPop2[,3]))
ind=subset(ind,!duplicated(ind[,1]))
ind=subset(ind,!duplicated(ind[,3]))
indPop2=ind
write.table(indPop2,file=paste("ind_",pop2,".test", sep=""),col.names=F,row.names=F)
write.table(t.hapsPop2,file=paste("t_", pop2,".haps",sep=""),col.names=F)
d_pop2 = data2haplohh(map_file=paste("ind_",pop2,".test",sep=""),hap_file=paste("t_",pop2,".haps", sep=""))
result_ihh_pop2 = scan_hh(d_pop2)
result_pop2 = ihh2ihs(result_ihh_pop2)
png(file=paste("ihsplot_",pop2,".png",sep=""))
ihsplot(result_pop2$res.ihs,plot.pval=TRUE,ylim.scan=2,main="iHS")
dev.off()
write.table(result_pop2$res.ihs,file=paste("res.ihs.",pop2,".txt",sep=""), sep="\t")
save.image(file="ihs.RData")


sig=rownames(result_pop2$res.ihs[result_pop2$res.ihs[,4] > 2 & result_pop2$res.ihs[result_pop2$res.ihs[,4] < -2] ,])
sig=subset(sig, sig != "NA")
sigPop2=sig

#make bification plot for specified marker
i= 1
for (marker in sigPop2){
  png(file=paste(marker,"_",pop2,".png", sep=""))
  bifurcation.diagram(d_pop2,mrk_foc=which(hapsPop2[,1]==marker))
  dev.off()
}
setwd("../")
hhpop1 = result_ihh_pop1
hhpop2 = result_ihh_pop2
rsb_pop1_pop2=ies2rsb(hhpop1,hhpop2,popname1=pop1,popname2=pop2, method="bilateral")
png(file=paste("rsbplot_",pop1,"_",pop2,".png",sep=""))
rsbplot(rsb_pop1_pop2$res.rsb,plot.pval="FALSE",ylim.scan=2,pch=16,main=NA)
dev.off()

#create Fst graphs
weirFst=read.table("./out.weir.fst",sep="\t", header=TRUE)
weirMean=mean(weirFst[,3][weirFst[,3] != "NaN"])
weirSd=sd(weirFst[,3][weirFst[,3] != "NaN"])
weirThresUpper = weirMean + 3 * weirSd
weirThresLower = weirMean - 3 * weirSd
hapFst=read.table("./out.hapmap.fst",sep="\t", header=TRUE)
hapMean=mean(hapFst[,3])
hapSd=sd(hapFst[,3])
hapThresUpper = hapMean + 3 * hapSd
hapThresLower = hapMean - 3 * hapSd
png(file=paste(pop1,"_",pop2,"_weir_fst.png",sep=""), height=450, width=1350)
plot(weirFst[,3]~weirFst[,2], type="p", xlab=chr, ylab="Weir Fst")
abline(h=weirThresUpper, lty=2, col = "red")
abline(h=weirThresLower, lty=2, col="red")
abline(h=weirMean, lty=2,col="green")
dev.off()
png(file=paste(pop1,"_",pop2,"_hapmap_fst.png",sep=""), height=450, width=1350)
plot(hapFst[,3]~hapFst[,2], type="p", xlab=chr, ylab="Hapmap Fst")
abline(h=hapThresUpper, lty=2, col = "red")
abline(h=hapThresLower, lty=2, col="red")
abline(h=hapMean, lty=2,col="green")
dev.off()

png(file=paste(pop1,"_",pop2,"combined_ihs_fst_rsb",".png",sep=""))
par(mfcol=c(5,1), mar=c(4,4,0.5,0.5), oma=c(1.5,2,1,1))
plot(hapFst[,3]~hapFst[,2], type="l", xlab=chr, ylab="Hapmap Fst")
abline(h=hapAverage, col = "red")
plot(weirFst[,3]~weirFst[,2], type="l", xlab=chr, ylab="Weir Fst")
ihsplot(result_pop1$res.ihs,plot.pval=FALSE,ylim.scan=2,main=paste("iHS_",pop1,sep=""))
ihsplot(result_pop2$res.ihs,plot.pval=FALSE,ylim.scan=2,main=paste("iHS_",pop2,sep=""))
rsbplot(rsb_pop1_pop2$res.rsb,plot.pval="FALSE",ylim.scan=2,pch=16,main=NA)
dev.off()

png(file=paste(pop1,"_",pop2,"_iHS_rsb_distributions.png",sep=""))
par(mfrow=c(3,1))
distribplot(result_pop1$res.ihs[,3],col=c("blue","red"), main=paste("iHS distribution",pop1,sep=" "),xlab="iHS")
distribplot(result_pop2$res.ihs[,3],col=c("blue","red"), main=paste("iHS distribution",pop2,sep=" "),xlab="iHS")
distribplot(rsb_pop1_pop2$res.rsb[,3],col=c("blue","red"), main=paste("rsb distribution",pop1,pop2,sep=" "),xlab="rsb")
dev.off()
save.image(file="ihs.RData")


