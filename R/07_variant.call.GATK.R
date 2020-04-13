
#' @title gatk.haplotypecaller
#' @description A wrapper function to run GATK (HaplotypeCaller)
#' @usage gatk.haplotypecaller(fns.bam, output.dir, sample.name, ref.fa, genotyping_mode="DISCOVERY", output_mode="EMIT_VARIANTS_ONLY", stand_call_conf_number=30, unsafe=FALSE, run.cmd=TRUE, mc.cores=1)
#' @param fns.bam Path to BAM files
#' @param output.dir Output directory
#' @param sample.name A character vector for the sample names
#' @param ref.fa Referance fasta file
#' @param genotyping_mode A parameter value for --genotyping_mode in GATK. A character vector to determine the alternate alleles to use for genotyping (default: DISCOVERY)
#' @param output_mode A parameter value for --output_mode in GATK. A character vector to produces variant calls (default: EMIT_VARIANTS_ONLY) 
#' @param stand_call_conf_number A parameter value for -stand_call_conf in GATK. A numeric value of The minimum phred-scaled confidence threshold at which variants should be called (default: 30)
#' @param unsafe A parameter value for -U ALLOW_N_CIGAR_READS in GATK. This parameter must be TRUE in RNA-seq data. 
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details HaplotypeCaller is capable of calling SNPs and indels simultaneously via local de-novo assembly of haplotypes in an active 
#'          region. 
#' @return Variant calling format files (.vcf)
#' @references The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data.
#' @seealso \url{https://software.broadinstitute.org/gatk/}
#' @export
gatk.haplotypecaller=function(fns.bam, 
                              output.dir, 
                              sample.name, 
                              ref.fa, 
                              genotyping_mode= "DISCOVERY", output_mode= "EMIT_VARIANTS_ONLY", stand_call_conf_number= 30, 
                              unsafe=FALSE, 
                              run.cmd=TRUE,
                              mc.cores=1){
  
  out.dirs=file.path(output.dir, sample.name)
  out.fns=file.path(out.dirs, paste0(sample.name, ".vcf"))
  
  #run
  message("[[",Sys.time(),"]] Run haplotypecaller -------")
  
  cmd.add=ifelse(unsafe, "-U ALLOW_N_CIGAR_READS", "")
  cmd=paste(GATK4.path, "--java-options -Xmx64G HaplotypeCaller", "-R", ref.fa, "-I", fns.bam, "-O", out.fns, " --output-mode", output_mode, "-stand-call-conf", stand_call_conf_number, cmd.add)
  
  print_message(cmd)
  
  if(run.cmd) mclapply(cmd, system, mc.cores=mc.cores)
  cat(cmd, file=file.path(output.dir, "haplotypecall.run.log"), sep="\n", append = FALSE)
  
  out.fns
}


#' @title gatk.variantfilter
#' @description A wrapper function to run GATK (VariantFiltration)
#' @usage gatk.variantfilter(fns.vcf, output.dir, sample.name, ref.fa, FS=30.0, QD=2.0, QUAL=50, DP=5, gatk.window=35, cluster=3, run.cmd=TRUE, mc.cores=1)
#' @param fns.vcf Path to VCF files
#' @param output.dir Output directory
#' @param sample.name A character vector for the sample names
#' @param ref.fa Reference fasta file 
#' @param FS A parameter value for FS in GATK. FisherStrand. (default=30.0)
#' @param QD A parameter value for QD in GATK. Quality by Depth. (default=2.0)
#' @param QUAL A parameter value for QUAL in GATK. Low quality. (default=50)
#' @param DP A parameter value for DP in GATK. Low depth. (default=5)
#' @param gatk.window A parameter value for -window in GATK. The window size (in bases) in which to evaluate clustered SNPs.
#' @param cluster A parameter value for -cluster in GATK. The number of SNPs which make up a cluster. Must be at least 2. 
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details Filter variant calls based on INFO and/or FORMAT annotations. This tool is designed for hard-filtering variant calls based 
#'          on certain criteria. 
#' @return Filtered VCF file (eg. .f.vcf)
#' @references The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data.
#' @seealso \url{https://software.broadinstitute.org/gatk/}
#' @export
gatk.variantfilter=function(fns.vcf,
                            output.dir,
                            ref.fa,
                            FS=30.0,
                            QD=2.0,
                            QUAL=50,
                            DP=5,
                            gatk.window= 35,
                            cluster= 3,
                            run.cmd=TRUE,
                            mc.cores=1){
  
  FS=trimws(format(round(as.numeric(FS), 1), nsmall = 1))
  QD=trimws(format(round(as.numeric(QD), 1), nsmall = 1))
  
  out.fns=sub("vcf$", "f.vcf", fns.vcf)
  
  fil.opt=paste0("--filter-name \"FS\" --filter-expression", " \" FS > ", FS,"\"", " --filter-name \"QD\" --filter-expression", " \"QD < ", QD,"\"", " --filter-name \"LowQual\" --filter-expression ", "\"QUAL < ", QUAL,"\"", " --filter-name \"LowDepth\" --filter-expression ", "\"DP < ", DP,"\"")
  message("[[",Sys.time(),"]] Run variant filtration -------")
  cmd=paste(GATK4.path, "--java-options -Xmx64G VariantFiltration", "-R", ref.fa, "-V", fns.vcf, "-O", out.fns, "-window", gatk.window, "-cluster", cluster, fil.opt) 
  
  print_message(cmd)
  
  if(run.cmd) mclapply(cmd, system,mc.cores=mc.cores) 
  cat(cmd, file=file.path(output.dir, "variantfilt.run.log"), sep="\n", append = FALSE)
  
  out.fns
}
