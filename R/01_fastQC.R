
#' @title fastQC
#' @description A wrapper function to run fastQC
#' @usage fastqc(fq1, fq2, output.dir, run.cmd=TRUE)
#' @param fq1 Path to read1 fastq files
#' @param fq2 Path to read2 fastq files
#' @param output.dir Output directory
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details FastQC aims to provide a QC report that detects problems originating from either the sequencer or the starting library material.
#' @return Quality check report for sequence data. (e.g., .html)
#' @import parallel
#' @references FastQC: a quality control tool for high throughput sequence data. Andrews S. (2010). 
#' @seealso \url{http://www.bioinformatics.babraham.ac.uk/projects/fastqc}
#' @export
fastQC=function(fq1, fq2, output.dir, run.cmd=TRUE, mc.cores=1){
  
  fq2=ifelse(grepl(fq2.idx, fq2), fq2, "")
  # command
  cmd=paste(fastqc.path, "-o", output.dir, "--extract", fq1, fq2)
  
  # run
  message("[[",Sys.time(),"]] Run fastQC --------")
  
  print_message(cmd)
  
  if(run.cmd) mclapply(cmd, system, mc.cores=mc.cores)
  cat(cmd, file=file.path(output.dir, "run.fastqc.log"), sep="\n", append = FALSE)
}
