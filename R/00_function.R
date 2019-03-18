

#' @title get.fns
#' @description Gets path of input files
#' @usage get.fns(input.dir, idx)
#' @param input.dir Path to directory including input files 
#' @param idx Suffix of input files
#' @return Path to the input files
#' @import R.utils
#' @export
get.fns=function(input.dir, idx=".1.fastq$|.2.fastq$"){
  fns=dir(input.dir, full.names=TRUE, recursive=TRUE)
  fns=grep(idx, fns, value=TRUE)
  
  # library(R.utils)
  #if(!isDirectory(input.dir)) fns=input.dir
  fns
}



#' @title print_message 
#' @description Show command line
#' @usage print_message(cmd)
#' @param cmd What users want to show as a message 
#' @return message()
#' @export
print_message=function(cmd){
  lapply(cmd, message)
}



#' @title get.process.names
#' @description Process names to be used in report
#' @usage get.process.names(qc, trim.method, align.method, bwa.method , rm.dup, realign, variant.call.method, annotation.method, rseq.quant.method)
#' @param qc As the quality check progresses, qc is added to the report processes.
#' @param trim.method As the trimming progresses, trimming method is added to the report processes.
#' @param align.method As the alignment progresses, alignment method is added to the report processes.
#' @param bwa.method When alignment is performed with bwa, bwa method is added to the report processes
#' @param rm.dup As the removal of duplicates progresses, removal method is added to the report processes.
#' @param realign As the re-alignment progresses, re-alignment is added to the report processes.
#' @param variant.call.method As the variant calling progresses, variant calling method is added to the report processes.
#' @param annotation.methodd As the variant annotation progresses, annotation method is added to the report processes.
#' @param rseq.quant.method As the RNA quantification progresses, RNA quantification method is added to the report processes.
#' @export 
get.process.names=function(qc, trim.method, align.method, bwa.method , rm.dup, realign, variant.call.method, annotation.method, rseq.abundance.method, make.eSet){
  proc.names=NULL
  if(qc) proc.names=c(proc.names, "qc")
  if(trim.method=="trim.galore") proc.names=c(proc.names,"trim") else if(trim.method=="cutadapt") proc.names=c(proc.names,"cutadapt")
  if(align.method=="bwa"&bwa.method=="mem") proc.names=c(proc.names, "bwa-mem") else if(bwa.method=="aln") proc.names=c(proc.names, "bwa-aln") else if(align.method!="none") proc.names=c(proc.names,align.method) 
  if(rm.dup =="MarkDuplicates") proc.names=c(proc.names, "rmdu")else if(rm.dup =="BARCODE") proc.names=c(proc.names, "rmdu_b") 
  if(realign) proc.names=c(proc.names, "realign")
  if(variant.call.method!="none") proc.names=c(proc.names, variant.call.method)
  if(annotation.method=="annovar") proc.names=c(proc.names, "annot")else if(annotation.method=="vep") proc.names=c(proc.names, "vep")
  if(rseq.abundance.method!="none") proc.names=c(proc.names, rseq.abundance.method)
  if(make.eSet) proc.names=c(proc.names,"make.Set")
  proc.names
}
