
#' @title picard.rmdu
#' @description A wrapper function to run Picard (MarkDuplicates)
#' @usage picard.rmdu(fns.bam, output.dir, sample.name, out.metrics, CREATE_INDEX=TRUE, REMOVE_DUPLICATES=TRUE, VALIDATION_STRINGENCY="LENIENT", tmp.dir=file.path(output.dir, "tmp"), BARCODE_TAG=FALSE, run.cmd=TRUE, mc.cores=1)
#' @param fns.bam Path to BAM files
#' @param output.dir Output directory
#' @param sample.name A character vector for the sample names
#' @param CREATE_INDEX  A parameter value for CREATE_INDEX in picard. Logical, whether to create bam index files (default=TRUE)
#' @param REMOVE_DUPLICATES  A parameter value for REMOVE_DUPLICATES in picard. Logical, whether to remove duplicates (default=TRUE)
#' @param VALIDATION_STRINGENCY  A parameter value for VALIDATION_STRINGENCY in picard. A character value of validation stringency (default="LENIENT")
#' @param tmp.dir Temporary directory (default= ./tmp)
#' @param BARCODE_TAG  A parameter value for BARCODE_TAG in picard. If barcode sequencing data, set this option TRUE. Duplicated BARCODE is removed. 
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.#' 
#' @details The MarkDuplicates tool works by comparing sequences in the 5 prime positions of both reads and read-pairs in a 
#'          SAM/BAM file. An BARCODE_TAG option is available to facilitate duplicate marking using molecular barcodes. After 
#'          duplicate reads are collected, the tool differentiates the primary and duplicate reads using an algorithm that 
#'          ranks reads by the sums of their base-quality scores.
#' @return Removing duplicate reads in bam files (e.g., .rmdu.bam)
#' @import parallel
#' @seealso \url{http://broadinstitute.github.io/picard/}
#' @export
picard.rmdu=function(fns.bam, 
                     output.dir,
                     type,
                     sample.name, 
                     CREATE_INDEX=TRUE, 
                     REMOVE_DUPLICATES=TRUE,
                     VALIDATION_STRINGENCY="LENIENT",
                     tmp.dir=file.path(output.dir, "tmp"),
                     BARCODE_TAG=FALSE,
                     run.cmd=TRUE,
                     mc.cores=1){
  
  out.dirs=file.path(output.dir, sample.name)
  lapply(1:length(out.dirs), function(a) dir.create(out.dirs[a], recursive = TRUE, showWarnings=FALSE))
  out.fns=file.path(out.dirs, paste0(sample.name,".rmdu.bam"))
  out.metrics.fns=sub(".bam",".metrics", out.fns)
  
  cmd.add=ifelse(BARCODE_TAG, " BARCODE_TAG=RX", "")
  cmd=paste0("java -XX:ParallelGCThreads=4 -XX:ConcGCThreads=1 -Xmx8G -jar ", picard.path, " MarkDuplicates", " INPUT=", fns.bam, " OUTPUT=", out.fns, " METRICS_FILE=", out.metrics.fns, 
             " CREATE_INDEX=", CREATE_INDEX, " VALIDATION_STRINGENCY=", VALIDATION_STRINGENCY, 
             " REMOVE_DUPLICATES=", REMOVE_DUPLICATES, " TMP_DIR=", tmp.dir, cmd.add)
  
  #run
  message("[[",Sys.time(),"]] Run picard.rmdu---- ")
  print_message(cmd)
  
  if(run.cmd) mclapply(cmd, system,mc.cores=mc.cores)
  cat(cmd, file=file.path(output.dir, "run.rmdu.log"), sep="\n", append = FALSE)
  
  out.fns
}
