
#-


#' @title gatk.baserecalibrator
#' @description A wrapper function to run GATK (BaseRecalibrator)
#' @usage gatk.baserecalibrator(fns.bam, output.dir, sample.name, ref.fa, ref.dbSNP, ref.gold_indels, unsafe=FALSE, run.cmd=TRUE, mc.cores=1)
#' @param fns.bam Path to input BAM files
#' @param output.dir Output directory
#' @param sample.name A character vector for the sample names
#' @param ref.fa Reference fasta file
#' @param ref.dbSNP Known SNP sites reference(VCF)
#' @param ref.gold_indels Known sites to indel variant call format(VCF)
#' @param unsafe A parameter value for -U ALLOW_N_CIGAR_READS in GATK. This parameter must be TRUE in RNA-seq data. 
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details First pass of the base quality score recalibration. Generates a recalibration table based on various covariates. 
#'          The default covariates are read group, reported quality score, machine cycle, and nucleotide context. This walker 
#'          generates tables based on specified covariates via by-locus traversal operating only at sites that are in the known sites VCF. 
#' @return GATK report file contained recalibration table by read group, quality scores and all the optional covariates. (e.g., .grp)
#' @import parallel
#' @references The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data.
#' @seealso \url{https://software.broadinstitute.org/gatk/}
#' @export
gatk.baserecalibrator=function(fns.bam, 
                                 output.dir, 
                                 sample.name,
                                 
                                 ref.fa, 
                                 ref.dbSNP, 
                                 ref.gold_indels,
                                 
                                 unsafe=FALSE,
                                 gatk.thread.number=4,
                                 run.cmd=TRUE,
                                 mc.cores=1){
  
  out.dirs=file.path(output.dir, sample.name)
  lapply(1:length(out.dirs), function(a) dir.create(out.dirs[a], recursive = TRUE, showWarnings=FALSE))
  out.fns=file.path(out.dirs, paste0(sample.name, ".table"))
  
  #BaseRecalibrator
  message("[[", Sys.time(),"]] Run BaseRecalibrator -------")
  
  cmd.add=ifelse(unsafe, "-U ALLOW_N_CIGAR_READS", "")
  cmd=paste(GATK4.path, "--java-options -Xmx64G BaseRecalibrator", "-I", fns.bam, "-R", ref.fa, "--known-sites", ref.dbSNP, "--known-sites", ref.gold_indels, "-O", out.fns, cmd.add)
  
  print_message(cmd)
  
  if(run.cmd) mclapply(cmd, system,mc.cores=mc.cores)
  cat(cmd, file=file.path(output.dir, "baserecal.run.log"), sep="\n", append=FALSE)
  
  out.fns
}




#' @title gatk.applyBQSR
#' @description A wrapper function to run GATK (PrintReads)
#' @usage gatk.applyBQSR(fn.realign.bam, output.dir, sample.name, ref.fa, fns.grp, unsafe=FALSE, run.cmd=TRUE, mc.cores=1)
#' @param fns.bam Path to input BAM files
#' @param output.dir Output directory
#' @param sample.name A character vector for the sample names
#' @param ref.fa Reference fasta file
#' @param fns.grp GATK report file created by BaseRecalibrator
#' @param unsafe A parameter value for -U ALLOW_N_CIGAR_READS in GATK. This parameter must be TRUE in RNA-seq data. 
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details Writes a new file using reads from SAM format file (SAM/BAM/CRAM) that pass criteria. Improves the accuracy of variant 
#'          calling based on Base Quality Score Recalibration.
#' @return GATK PrintReads returns a Base quality score recalibrated BAM files (eg. recal.bam)
#' @references The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data.
#' @seealso \url{https://software.broadinstitute.org/gatk/}
#' @export
gatk.applyBQSR=function(fns.bam, 
                           fns.grp,
                           output.dir,
                           sample.name,
                           ref.fa,
                           unsafe=FALSE,
                           gatk.thread.number=4,
                           run.cmd=TRUE,
                           mc.cores=1){
  
  out.dirs=file.path(output.dir, sample.name)
  lapply(1:length(out.dirs), function(a) dir.create(out.dirs[a], recursive = TRUE, showWarnings=FALSE))
  out.fns= file.path(out.dirs, paste0(sample.name, ".recal.bam"))
  
  #PrintReads
  message("[[", Sys.time(),"]] Run PrintReads-------")
  
  cmd.add=ifelse(unsafe, "-U ALLOW_N_CIGAR_READS", "")
  cmd=paste(GATK4.path, "--java-options -Xmx64G ApplyBQSR", "-I", fns.bam, "-R", ref.fa, "--bqsr-recal-file", fns.grp, "-O", out.fns, cmd.add)
  
  print_message(cmd)
  
  if(run.cmd) mclapply(cmd, system, mc.cores=mc.cores)
  cat(cmd, file=file.path(output.dir, "applyBQSR.run.log"), sep="\n", append=FALSE)
  
  out.fns
}

#' @title gatk.baserecalibrator.o
#' @description A wrapper function to run GATK (BaseRecalibrator)
#' @usage gatk.baserecalibrator(fns.bam, output.dir, sample.name, ref.fa, ref.dbSNP, ref.gold_indels, unsafe=FALSE, run.cmd=TRUE, mc.cores=1)
#' @param fns.bam Path to input BAM files
#' @param output.dir Output directory
#' @param sample.name A character vector for the sample names
#' @param ref.fa Reference fasta file
#' @param ref.dbSNP Known SNP sites reference(VCF)
#' @param ref.gold_indels Known sites to indel variant call format(VCF)
#' @param unsafe A parameter value for -U ALLOW_N_CIGAR_READS in GATK. This parameter must be TRUE in RNA-seq data. 
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details First pass of the base quality score recalibration. Generates a recalibration table based on various covariates. 
#'          The default covariates are read group, reported quality score, machine cycle, and nucleotide context. This walker 
#'          generates tables based on specified covariates via by-locus traversal operating only at sites that are in the known sites VCF. 
#' @return GATK report file contained recalibration table by read group, quality scores and all the optional covariates. (e.g., .grp)
#' @import parallel
#' @references The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data.
#' @seealso \url{https://software.broadinstitute.org/gatk/}
#' @export
gatk.baserecalibrator.o=function(fns.bam, 
                                 output.dir, 
                                 sample.name,
                                 
                                 ref.fa, 
                                 ref.dbSNP, 
                                 ref.gold_indels,
                                 
                                 unsafe=FALSE,
                                 gatk.thread.number=4,
                                 run.cmd=TRUE,
                                 mc.cores=1){
  
  out.dirs=file.path(output.dir, sample.name)
  lapply(1:length(out.dirs), function(a) dir.create(out.dirs[a], recursive = TRUE, showWarnings=FALSE))
  out.fns=file.path(out.dirs, paste0(sample.name, ".grp"))
  
  #BaseRecalibrator
  message("[[", Sys.time(),"]] Run BaseRecalibrator -------")
  
  cmd.add=ifelse(unsafe, "-U ALLOW_N_CIGAR_READS", "")
  cmd=paste("java -jar", GATK.path, "-T BaseRecalibrator", "-I", fns.bam, "-R", ref.fa, "-knownSites", ref.dbSNP, "-knownSites", ref.gold_indels, "-o", out.fns, "-nct", gatk.thread.number, cmd.add)
  
  print_message(cmd)
  
  if(run.cmd) mclapply(cmd, system,mc.cores=mc.cores)
  cat(cmd, file=file.path(output.dir, "baserecal.run.log"), sep="\n", append=FALSE)
  
  out.fns
}



#' @title gatk.printreads.o
#' @description A wrapper function to run GATK (PrintReads)
#' @usage gatk.printreads(fn.realign.bam, output.dir, sample.name, ref.fa, fns.grp, unsafe=FALSE, run.cmd=TRUE, mc.cores=1)
#' @param fns.bam Path to input BAM files
#' @param output.dir Output directory
#' @param sample.name A character vector for the sample names
#' @param ref.fa Reference fasta file
#' @param fns.grp GATK report file created by BaseRecalibrator
#' @param unsafe A parameter value for -U ALLOW_N_CIGAR_READS in GATK. This parameter must be TRUE in RNA-seq data. 
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details Writes a new file using reads from SAM format file (SAM/BAM/CRAM) that pass criteria. Improves the accuracy of variant 
#'          calling based on Base Quality Score Recalibration.
#' @return GATK PrintReads returns a Base quality score recalibrated BAM files (eg. recal.bam)
#' @references The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data.
#' @seealso \url{https://software.broadinstitute.org/gatk/}
#' @export
gatk.printreads.o=function(fns.bam, 
                           fns.grp,
                           output.dir,
                           sample.name,
                           ref.fa,
                           unsafe=FALSE,
                           gatk.thread.number=4,
                           run.cmd=TRUE,
                           mc.cores=1){
  
  out.dirs=file.path(output.dir, sample.name)
  lapply(1:length(out.dirs), function(a) dir.create(out.dirs[a], recursive = TRUE, showWarnings=FALSE))
  out.fns= file.path(out.dirs, paste0(sample.name, ".recal.bam"))
  
  #PrintReads
  message("[[", Sys.time(),"]] Run PrintReads-------")
  
  cmd.add=ifelse(unsafe, "-U ALLOW_N_CIGAR_READS", "")
  cmd=paste("java -jar", GATK.path, "-T PrintReads", "-I", fns.bam, "-R", ref.fa, "-BQSR", fns.grp, "-o", out.fns, "-nct", gatk.thread.number, cmd.add)
  
  print_message(cmd)
  
  if(run.cmd) mclapply(cmd, system, mc.cores=mc.cores)
  cat(cmd, file=file.path(output.dir, "printreads.run.log"), sep="\n", append=FALSE)
  
  out.fns
}




#' @title gatk.baserecalibrator.o
#' @description A wrapper function to run GATK (BaseRecalibrator)
#' @usage gatk.baserecalibrator(fns.bam, output.dir, sample.name, ref.fa, ref.dbSNP, ref.gold_indels, unsafe=FALSE, run.cmd=TRUE, mc.cores=1)
#' @param fns.bam Path to input BAM files
#' @param output.dir Output directory
#' @param sample.name A character vector for the sample names
#' @param ref.fa Reference fasta file
#' @param ref.dbSNP Known SNP sites reference(VCF)
#' @param ref.gold_indels Known sites to indel variant call format(VCF)
#' @param unsafe A parameter value for -U ALLOW_N_CIGAR_READS in GATK. This parameter must be TRUE in RNA-seq data. 
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details First pass of the base quality score recalibration. Generates a recalibration table based on various covariates. 
#'          The default covariates are read group, reported quality score, machine cycle, and nucleotide context. This walker 
#'          generates tables based on specified covariates via by-locus traversal operating only at sites that are in the known sites VCF. 
#' @return GATK report file contained recalibration table by read group, quality scores and all the optional covariates. (e.g., .grp)
#' @import parallel
#' @references The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data.
#' @seealso \url{https://software.broadinstitute.org/gatk/}
#' @export
gatk.baserecalibrator.o=function(fns.bam, 
                                 output.dir, 
                                 sample.name,
                                 
                                 ref.fa, 
                                 ref.dbSNP, 
                                 ref.gold_indels,
                                 
                                 unsafe=FALSE,
                                 gatk.thread.number=4,
                                 run.cmd=TRUE,
                                 mc.cores=1){
  
  out.dirs=file.path(output.dir, sample.name)
  lapply(1:length(out.dirs), function(a) dir.create(out.dirs[a], recursive = TRUE, showWarnings=FALSE))
  out.fns=file.path(out.dirs, paste0(sample.name, ".grp"))
  
  #BaseRecalibrator
  message("[[", Sys.time(),"]] Run BaseRecalibrator -------")
  
  cmd.add=ifelse(unsafe, "-U ALLOW_N_CIGAR_READS", "")
  cmd=paste("java -jar", GATK.path, "-T BaseRecalibrator", "-I", fns.bam, "-R", ref.fa, "-knownSites", ref.dbSNP, "-knownSites", ref.gold_indels, "-o", out.fns, "-nct", gatk.thread.number, cmd.add)
  
  print_message(cmd)
  
  if(run.cmd) mclapply(cmd, system,mc.cores=mc.cores)
  cat(cmd, file=file.path(output.dir, "baserecal.run.log"), sep="\n", append=FALSE)
  
  out.fns
}



#' @title gatk.printreads.o
#' @description A wrapper function to run GATK (PrintReads)
#' @usage gatk.printreads(fn.realign.bam, output.dir, sample.name, ref.fa, fns.grp, unsafe=FALSE, run.cmd=TRUE, mc.cores=1)
#' @param fns.bam Path to input BAM files
#' @param output.dir Output directory
#' @param sample.name A character vector for the sample names
#' @param ref.fa Reference fasta file
#' @param fns.grp GATK report file created by BaseRecalibrator
#' @param unsafe A parameter value for -U ALLOW_N_CIGAR_READS in GATK. This parameter must be TRUE in RNA-seq data. 
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details Writes a new file using reads from SAM format file (SAM/BAM/CRAM) that pass criteria. Improves the accuracy of variant 
#'          calling based on Base Quality Score Recalibration.
#' @return GATK PrintReads returns a Base quality score recalibrated BAM files (eg. recal.bam)
#' @references The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data.
#' @seealso \url{https://software.broadinstitute.org/gatk/}
#' @export
gatk.printreads.o=function(fns.bam, 
                           fns.grp,
                           output.dir,
                           sample.name,
                           ref.fa,
                           unsafe=FALSE,
                           gatk.thread.number=4,
                           run.cmd=TRUE,
                           mc.cores=1){
  
  out.dirs=file.path(output.dir, sample.name)
  lapply(1:length(out.dirs), function(a) dir.create(out.dirs[a], recursive = TRUE, showWarnings=FALSE))
  out.fns= file.path(out.dirs, paste0(sample.name, ".recal.bam"))
  
  #PrintReads
  message("[[", Sys.time(),"]] Run PrintReads-------")
  
  cmd.add=ifelse(unsafe, "-U ALLOW_N_CIGAR_READS", "")
  cmd=paste("java -jar", GATK.path, "-T PrintReads", "-I", fns.bam, "-R", ref.fa, "-BQSR", fns.grp, "-o", out.fns, "-nct", gatk.thread.number, cmd.add)
  
  print_message(cmd)
  
  if(run.cmd) mclapply(cmd, system, mc.cores=mc.cores)
  cat(cmd, file=file.path(output.dir, "printreads.run.log"), sep="\n", append=FALSE)
  
  out.fns
}



#' @title gatk.baserecalibrator.o
#' @description A wrapper function to run GATK (BaseRecalibrator)
#' @usage gatk.baserecalibrator(fns.bam, output.dir, sample.name, ref.fa, ref.dbSNP, ref.gold_indels, unsafe=FALSE, run.cmd=TRUE, mc.cores=1)
#' @param fns.bam Path to input BAM files
#' @param output.dir Output directory
#' @param sample.name A character vector for the sample names
#' @param ref.fa Reference fasta file
#' @param ref.dbSNP Known SNP sites reference(VCF)
#' @param ref.gold_indels Known sites to indel variant call format(VCF)
#' @param unsafe A parameter value for -U ALLOW_N_CIGAR_READS in GATK. This parameter must be TRUE in RNA-seq data. 
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details First pass of the base quality score recalibration. Generates a recalibration table based on various covariates. 
#'          The default covariates are read group, reported quality score, machine cycle, and nucleotide context. This walker 
#'          generates tables based on specified covariates via by-locus traversal operating only at sites that are in the known sites VCF. 
#' @return GATK report file contained recalibration table by read group, quality scores and all the optional covariates. (e.g., .grp)
#' @import parallel
#' @references The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data.
#' @seealso \url{https://software.broadinstitute.org/gatk/}
#' @export
gatk.baserecalibrator.o=function(fns.bam, 
                               output.dir, 
                               sample.name,
                               
                               ref.fa, 
                               ref.dbSNP, 
                               ref.gold_indels,
                               
                               unsafe=FALSE,
                               gatk.thread.number=4,
                               run.cmd=TRUE,
                               mc.cores=1){
  
  out.dirs=file.path(output.dir, sample.name)
  lapply(1:length(out.dirs), function(a) dir.create(out.dirs[a], recursive = TRUE, showWarnings=FALSE))
  out.fns=file.path(out.dirs, paste0(sample.name, ".grp"))
  
  #BaseRecalibrator
  message("[[", Sys.time(),"]] Run BaseRecalibrator -------")
  
  cmd.add=ifelse(unsafe, "-U ALLOW_N_CIGAR_READS", "")
  cmd=paste("java -jar", GATK.path, "-T BaseRecalibrator", "-I", fns.bam, "-R", ref.fa, "-knownSites", ref.dbSNP, "-knownSites", ref.gold_indels, "-o", out.fns, "-nct", gatk.thread.number, cmd.add)
  
  print_message(cmd)
  
  if(run.cmd) mclapply(cmd, system,mc.cores=mc.cores)
  cat(cmd, file=file.path(output.dir, "baserecal.run.log"), sep="\n", append=FALSE)
  
  out.fns
}



#' @title gatk.printreads.o
#' @description A wrapper function to run GATK (PrintReads)
#' @usage gatk.printreads(fn.realign.bam, output.dir, sample.name, ref.fa, fns.grp, unsafe=FALSE, run.cmd=TRUE, mc.cores=1)
#' @param fns.bam Path to input BAM files
#' @param output.dir Output directory
#' @param sample.name A character vector for the sample names
#' @param ref.fa Reference fasta file
#' @param fns.grp GATK report file created by BaseRecalibrator
#' @param unsafe A parameter value for -U ALLOW_N_CIGAR_READS in GATK. This parameter must be TRUE in RNA-seq data. 
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details Writes a new file using reads from SAM format file (SAM/BAM/CRAM) that pass criteria. Improves the accuracy of variant 
#'          calling based on Base Quality Score Recalibration.
#' @return GATK PrintReads returns a Base quality score recalibrated BAM files (eg. recal.bam)
#' @references The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data.
#' @seealso \url{https://software.broadinstitute.org/gatk/}
#' @export
gatk.printreads.o=function(fns.bam, 
                         fns.grp,
                         output.dir,
                         sample.name,
                         ref.fa,
                         unsafe=FALSE,
                         gatk.thread.number=4,
                         run.cmd=TRUE,
                         mc.cores=1){
  
  out.dirs=file.path(output.dir, sample.name)
  lapply(1:length(out.dirs), function(a) dir.create(out.dirs[a], recursive = TRUE, showWarnings=FALSE))
  out.fns= file.path(out.dirs, paste0(sample.name, ".recal.bam"))
  
  #PrintReads
  message("[[", Sys.time(),"]] Run PrintReads-------")
  
  cmd.add=ifelse(unsafe, "-U ALLOW_N_CIGAR_READS", "")
  cmd=paste("java -jar", GATK.path, "-T PrintReads", "-I", fns.bam, "-R", ref.fa, "-BQSR", fns.grp, "-o", out.fns, "-nct", gatk.thread.number, cmd.add)
  
  print_message(cmd)
  
  if(run.cmd) mclapply(cmd, system, mc.cores=mc.cores)
  cat(cmd, file=file.path(output.dir, "printreads.run.log"), sep="\n", append=FALSE)
  
  out.fns
}

