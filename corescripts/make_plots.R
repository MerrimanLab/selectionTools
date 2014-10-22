#!/bin/env Rscript

require(getopt)

spec = matrix(c(
    'help',     'h',    0,  "logical",
    'plot_type','p',    1,  "character",
    'data_file', 'd',   1,  "character",
    'plot_output_name', 'n', 1, 'character'
 ), byrow=T, ncol=4)

opt = getopt(spec)

if (!is.null(opt$help)){
    cat(getopt(spec,usage=TRUE));
    q(status=1);
}

if (opt$plot_type == "taj"){
    tajimaD=read.table(file=opt$data_file, header=TRUE)
    png(opt$plot_output_name)
    plot(tajimaD[,4] ~ tajimaD[,2],pch='.',cex=2,xlab="Chromosome position (bp)", ylab="D statistic")
    dev.off()
}else if( opt$plot_type == "ihs"){
    ihs = read.table(file=opt$data_file)
    png(opt$plot_output_name)
    plot(ihs[,4] ~ ihs[,2],,pch='.',cex=2,ylab=expression("-" * log[10] * "[" ~ "1-2|" * Phi[scriptstyle(italic(iHS))] * "-0.5|" ~ "]"),
         xlab="Chromosome Position BP")
    dev.off()
}else if(opt$plot_type == "fay"){
    fay = read.table(file=opt$data_file,comment.char="#") 
    png(opt$plot_output_name)
    plot(CEUFay[,15] ~ CEUFay[,1],xlab='Chromosome position (bp)',ylab="H Statistic"))
    dev.off()
}else if(opt$plot_type == "rsb"){

}else if(opt$plot_type == "fst_weir"){

}else if(opt$plot_type == "fst_hap"){

}

