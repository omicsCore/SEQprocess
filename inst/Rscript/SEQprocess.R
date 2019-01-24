#!/usr/bin/Rscript --vanilla
suppressPackageStartupMessages(library(SEQprocess))

argnames<-commandArgs(trailingOnly=TRUE)
if(length(argnames)==0) {
  system(paste0(path.package("SEQprocess"), "/Rscript/SEQprocess.R", " -h"))
  stop("Required arguments must be supplied")
}

# -h, --help option
source(system.file("Rscript", "functions_help.R", package="SEQprocess"))
get.option()
argList=as.list(grep("--",argnames, value=T, invert=T))
names(argList)=sub("--","",grep("--",argnames, value=T, invert=F))

######### Run #########
SEQprocess(argList = argList)
