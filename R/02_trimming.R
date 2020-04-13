


#' @title trim.gal
#' @description A wrapper function to run Trim Galore
#' @usage trim.gal(fq1, fq2, trim.quality=30, trim.clip_R1=13, trim.clip_R2=13, output.dir, run.cmd=TRUE, mc.cores=1)
#' @param fq1 Path to read1 fastq files
#' @param fq2 Path to read2 fastq files
#' @param output.dir Output directory
#' @param trim.quality A parameter value for --quality in trimgalore. A numeric value of phred score cutoff to trim (default=30)
#' @param trim.clip_R1 A parameter value for --clip_R1 in trimgalore. A numeric value of bp to remove adaptor from the 5-prime end of read1 files (default=13)
#' @param trim.clip_R2 A parameter value for --clip_R2 in trimgalore. A numeric value of bp to remove adaptor from the 5-prime end of read2 files (default=13)
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details Trims low quality bases, and cleans up adapter sequences for paired-end files.
#' @return Trimmed fastq files (e.g., .val_1.fastq, .val_1.fastq)
#' @import parallel
#' @references Krueger F. Trim Galore!
#' @seealso \url{https://www.bioinformatics.babraham.ac.uk/projects/trim_galore}
#' @export
trim.gal=function(fq1, fq2, output.dir, trim.quality=30, trim.clip_R1=13, trim.clip_R2=13, run.cmd=TRUE, mc.cores=1){
  
  sample.name=sub(fq1.idx, "", basename(fq1))
  out.dirs=file.path(output.dir, sample.name)
  lapply(1:length(out.dirs), function(a) dir.create(out.dirs[a], recursive = TRUE, showWarnings=FALSE))
  
  # command
  cmd=paste(trim_galore.path, "--quality", trim.quality, "--paired", "--clip_R1", trim.clip_R1, "--clip_R2", trim.clip_R2, "-o", out.dirs, "--dont_gzip", fq1, fq2)
  
  # run
  message("[[",Sys.time(),"]] Run trim_galore --------")
  print_message(cmd)
  
  if(run.cmd) mclapply(cmd, system, mc.cores = mc.cores)
  cat(cmd, file=file.path(output.dir, "run.trim.log"), sep="\n", append = FALSE)
  
  fq.pair=list(fq1=get.fns(output.dir, fq1.idx), fq2=get.fns(output.dir, fq2.idx))
  fq.pair
}


#' @title cutadapt
#' @description A wrapper function to run Cutadapt
#' @usage cutadpat(fq1, output.dir, adpat.seq="TGGAATTCTCGGGTGCCAAGG", m=1, mc.cores=1, run.cmd=TRUE)
#' @param fq1 Path to fastq files
#' @param output.dir Output directory
#' @param adapt.seq A parameter value for -b in cutadapt. Adapter sequences user wants to remove
#' @param m A parameter value for -m in cutadapt. Discards processed reads that are shorter than m option. Reads that are too short before adapter removal are also discarded.
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details Cutadapt finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequence.
#' @import parallel
#' @references Cutadapt Removes Adapter Sequences from High-Throughput Sequencing Reads
#' @seealso \url{https://cutadapt.readthedocs.io/en/stable/}
#' @export
cutadapt=function(fq1,
                  output.dir,
                  sample.name,
                  adapt.seq="TGGAATTCTCGGGTGCCAAGG", 
                  m=1,
                  run.cmd=TRUE,
                  mc.cores=1){
  
  out.dirs=file.path(output.dir, sample.name)
  lapply(1:length(out.dirs), function(a) dir.create(out.dirs[a], recursive = TRUE, showWarnings=FALSE))
  
  fq.fn=file.path(out.dirs, paste0(sample.name, ".1_val_1.fq"))
  cutadapt.txt=file.path(out.dirs, paste0(sample.name, ".cutadapt.txt"))
  
  cmd=paste0(cutadapt.path, " -b ", adapt.seq, " -m ", m, " -f fastq ", fq1 ," > ", fq.fn ," 2> ",cutadapt.txt)
  message("[[",Sys.time(),"]] Cut adapter sequence ---- ")
  print_message(cmd)
  
  if(run.cmd) mclapply(cmd, system,mc.cores=mc.cores)
  cat(cmd, file=file.path(output.dir, "run.cutadapt.log"), sep="\n", append = FALSE)
  
  fq.fn
}
