
#' @title picard.collectmetrics
#' @description Provide read alignment information.
#' @usage picard.collectmetrics(fns.bam, out.fns, ref.fa, run.cmd=T, mc.cores=1)
#' @param fns.bam Path to BAM files
#' @param output.dir Output directory
#' @param ref.fa Reference fasta file path
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details Provides a summary of the alignment process using the bam files
#' @return txt files
#' @import parallel
#' @author Ji-Hye Lee
#' @references Currently there is no journal reference for picard.
#' @seealso \url{http://broadinstitute.github.io/picard/}
#' @export
picard.collectmetrics=function(fns.bam, 
                               ref.fa,
                               run.cmd=TRUE,
                               mc.cores=1){
  
  out.fns=sub(".bam", "", fns.bam)
  #command
  cmd=paste0("java -jar ", picard.path, " CollectMultipleMetrics", " INPUT=", fns.bam, " OUTPUT=", out.fns, " REFERENCE_SEQUENCE=", ref.fa, " PROGRAM=null PROGRAM=CollectAlignmentSummaryMetrics")
  
  #run
  message("[[",Sys.time(),"]] Run picard.CollectMultipleMetrics ---- ")
  print_message(cmd)
  if(run.cmd) mclapply(cmd, system,mc.cores=mc.cores)
  output.dir=sub(paste0("/",basename(fns.bam[1])),"",file.path(fns.bam[1]))
  cat(cmd, file=file.path(output.dir, "run.CollectMultipleMetrics.log"), sep="\n")
  
  collect.output=paste0(out.fns,".alignment_summary_metrics")
  collect.output
}

