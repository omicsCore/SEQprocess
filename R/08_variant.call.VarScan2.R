
#' @title multiple.reads.pileup
#' @description A wrapper function to run samtools (mpileup)
#' @usage multiple.reads.pileup(ref.fa, normal.bam, tumor.bam, sample.name, output.dir, mapQ=1, run.cmd=TRUE, mc.cores=1)
#' @param ref.fa Reference fasta file path
#' @param normal.bam Path to normal sample recalibration bam files
#' @param tumor.bam Path to tumor sample recalibration bam files as tumor-normal pair
#' @param sample.name A character vector for the sample names
#' @param output.dir Output directory
#' @param mapQ A parameter value for mapQ in varscan2. Skip alignments with mapQ smaller than mapQ value (default:1)
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @import parallel
#' @details Generate VCF, BCF or pileup for one or multiple BAM files. Alignment records are grouped by sample (SM) identifiers in @RG header 
#'          lines. If sample identifiers are absent, each input file is regarded as one sample.
#' @seealso \url {http://www.htslib.org/doc/samtools.html}
#' @export
multiple.reads.pileup=function(ref.fa,
                               normal.bam,
                               tumor.bam,
                               sample.name,
                               output.dir,
                               run.cmd=TRUE,
                               mc.cores=1,
                               mapQ=1){
  
  out.dirs=file.path(output.dir, sample.name)
  lapply(1:length(out.dirs), function(a) dir.create(out.dirs[a], recursive = TRUE, showWarnings=FALSE))
  out.fns=file.path(out.dirs, paste0(sample.name, ".pileup"))
  
  #command-line
  cmd = paste(samtools.path, "mpileup", "-f", ref.fa, "-q", mapQ, "-B", normal.bam, tumor.bam, ">", out.fns)
  #run
  message("[[",Sys.time(),"]] Run SAMtools pileup ---- ")
  if(run.cmd) mclapply(cmd, system, mc.cores=mc.cores)
  cat(cmd, file=file.path(output.dir, "run.samtools.pileup.log"), sep="\n", append = FALSE)
  #return
  out.fns
}



#' @title varscan
#' @description A wrapper function to run VarScan2
#' @usage varscan(fn.pileup, output.dir, sample.name, min_coverage_normal=8, min_coverage_tumor=6, min_var_freq=0.10, min_freq_for_hom=0.75, somatic_p_value=0.05, strand_filter=0, run.cmd=TRUE, mc.cores=1)
#' @param fn.pileup samtools mpileup output file path
#' @param output.dir Output directory
#' @param sample.name A character vector for the sample names
#' @param min_coverage_normal A parameter value for --min-coverage-normal in VarScan2. Minimum coverage in normal to call somatic (default:8)
#' @param min_coverage_tumor A parameter value for --min-coverage-tumor in VarScan2. Minimum coverage in tumor to call somatic (default:6)
#' @param min_var_freq A parameter value for --min-var-freq in VarScan2. Minimum variant frequency to call a heterozygote (default:0.10)
#' @param min_freq_for_hom A parameter value for --min-freq-for-hom in VarScan2. Minimum frequency to call homozygote (default:0.75)
#' @param somatic_p_value A parameter value for --somatic-p-value in VarScan2. P-value threshold to call a somatic site (default:0.05)
#' @param strand_filter A parameter value for --strand-fiter in VarScan2. If set to 1, removes variants with > 90 percent strand bias(default:0)
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @import parallel
#' @details VarScan is a platform-independent mutation caller for targeted, exome, and whole-genome sequencing data.
#' @references VarScan 2: somatic mutation and copy number alteration discovery in cancer by exome sequencing.
#' @seealso \url {http://varscan.sourceforge.net/}
#' @export
varscan = function(fn.pileup, 
                   output.dir,
                   sample.name,
                   min_coverage_normal=8,
                   min_coverage_tumor=6, 
                   min_var_freq=0.10,
                   min_freq_for_hom=0.75, 
                   somatic_p_value=0.05,
                   strand_filter=0,
                   
                   run.cmd=TRUE,
                   mc.cores=1){
  
  out.dirs=file.path(output.dir, sample.name)
  lapply(1:length(out.dirs), function(a) dir.create(out.dirs[a], recursive = TRUE, showWarnings=FALSE))
  out.fns=file.path(out.dirs, sample.name)
  
  #command-line
  cmd = paste("java -jar", varscan.path,
              "somatic",
              fn.pileup,
              out.fns,
              "--mpileup", 1,
              "--min-coverage-normal", min_coverage_normal,
              "--min-coverage-tumor", min_coverage_tumor,
              "--min-var-freq", min_var_freq,
              "--min-freq-for-hom", min_freq_for_hom,
              "--normal-purity", 1,
              "--tumor-purity", 1,
              "--p-value", 0.99,
              "--somatic-p-value", somatic_p_value,
              "--strand-filter", strand_filter,
              " --output-vcf")
  #run
  message("[[",Sys.time(),"]] Run varscan ---- ")
  print_message(cmd)
  
  if(run.cmd) mclapply(cmd, system, mc.cores=mc.cores)
  cat(cmd, file=file.path(output.dir, "run.varscan.log"), sep="\n", append = FALSE)
  #return
  out.fns
}



#' @title processSomatic
#' @description A wrapper function to run VarScan2 (processSomatic)
#' @usage processSomatic(fns.vcf, output.dir, min_tumor_freq=0.1, max_normal_freq=0.05, p_value=0.07, run.cmd=TRUE, mc.cores=1)
#' @param fns.vcf varscan output file path
#' @param output.dir Output directory
#' @param min_tumor_freq A parameter value for --min-tumor-freq in varscan2. Minimum variant allele frequency in tumor (default:0.10)
#' @param max_normal_freq A parameter value for --max-normal-freq in varscan2. Maximum variant allele frequency in normal (default:0.05)
#' @param p_value A parameter value for --p-value in varscan2. P-value for high-confidence calling (default:0.07)
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @import parallel
#' @details processSomatic will separate a somatic output file by somatic_status (Germline, Somatic, LOH). Somatic mutations will 
#'          further be classified as high-confidence (.hc) or low-confidence (.lc).
#' @references VarScan 2: somatic mutation and copy number alteration discovery in cancer by exome sequencing.
#' @seealso \url {http://varscan.sourceforge.net/}
#' @export
processSomatic = function(fns.vcf,
                          output.dir,
                          min_tumor_freq=0.1,
                          max_normal_freq=0.05,
                          p_value=0.07,
                          
                          run.cmd = TRUE,
                          mc.cores = 1){
  
  #command-line
  cmd = paste("java -jar", varscan.path,
              "processSomatic",
              fns.vcf,
              "--min-tumor-freq", min_tumor_freq,
              "--max-normal-freq", max_normal_freq,
              "--p-value", p_value)
  #run
  message("[[",Sys.time(),"]] Run varscan ProcessSomatic---- ")
  print_message(cmd)
  
  if(run.cmd) mclapply(cmd, system, mc.cores=mc.cores)
  cat(cmd, file=file.path(output.dir, "run.processSomatic.log"), sep="\n", append = FALSE)
  
}
