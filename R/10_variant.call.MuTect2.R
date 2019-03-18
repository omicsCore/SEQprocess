

#' @title gatk.mutect2.normal
#' @description A wrapper function to run GATK (MuTect2)
#' @usage gatk.mutect2(normal.bam, sample.name, ref.dbSNP, cosmic.vcf, output.dir, run.cmd=TRUE, mc.cores=1)
#' @param normal.bam BAM files of normal samples
#' @param sample.name A character vector for the sample names
#' @param ref.dbSNP Known SNP sites reference vcf
#' @param cosmic.vcf Known variant sites of cosmic database vcf file
#' @param output.dir Output directory
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details MuTect2 is a somatic SNP and indel caller that combines the DREAM challenge-winning somatic genotyping engine of the original 
#'          MuTect (Cibulskis et al., 2013) with the assembly-based machinery of HaplotypeCaller. This function takes normal samples as 
#'          input to make the panel of normal (pon).
#' @return Only normal sample vcf files.
#' @references Sensitive detection of somatic point mutations in impure and heterogeneous cancer samples
#' @seealso \url{https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_cancer_m2_MuTect2.php}
#' @export
gatk.mutect2.normal=function(normal.bam,
                      sample.name,
                      ref.dbSNP,
                      cosmic.vcf,
                      output.dir,
                      run.cmd=TRUE,
                      mc.cores=1
                      ){
  #normal vcf list
  out.dirs=file.path(output.dir, sample.name)
  lapply(1:length(out.dirs), function(a) dir.create(out.dirs[a], recursive = TRUE, showWarnings=FALSE))
  out.vcf=file.path(out.dirs, paste0(sample.name, ".normal.vcf"))
  
  cmd=paste("java", "-jar", "/data/program/gatk/3.7/GenomeAnalysisTK.jar",
            "-T", "MuTect2",
            "-R", ref.fa,
            "-I:normal", normal.bam,
            "--dbsnp", ref.dbSNP,
            "--cosmic", cosmic.vcf,
            "-o", out.vcf)
  
  message("[[",Sys.time(),"]] Run GATK MuTect2 Normal only mode --------")
  print_message(cmd)
  
  if(run.cmd) mclapply(cmd, system, mc.cores=mc.cores)
  cat(cmd, file=file.path(output.dir, "run.mutect.log"), sep="\n", append=FALSE)
  out.vcf
}



#' @title gatk.combinevariants
#' @description A wrapper function to run GATK (CombineVariants)
#' @usage gatk.combinevariants(ref.fa, normal.vcf, minN=2, filteredrecordsmergetype="KEEP_IF_ANY_UNFILTERED", output.dir, run.cmd=TRUE, mc.cores=1)
#' @param ref.fa Reference fasta file
#' @param normal.vcf Normal sample vcfs list
#' @param minN Parameter value for -minN in GATK CombineVariants. Minimum number of samples to call the variant (default=2)
#' @param filteredrecordsmergetype  A parameter value for --filteredrecordsmergetype in GATK CombineVariants. Determines how to handle records seen at the same site in the VCF
#' @param output.dir Output directory
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details The MuTect2 pipeline employs a "Panel of Normal" to identify additional germline mutations. This method enables a higher 
#'          level of confidence to be assigned to somatic variants that are called by the MuTect2 pipeline.
#' @return pon(panel of normal) vcf file
#' @references Sensitive detection of somatic point mutations in heterogeneous cancer samples
#' @seealso \url{https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_cancer_m2_MuTect2.php}
#' @export
gatk.combinevariants=function(ref.fa,
                              normal.vcf,
                              minN=2,
                              filteredrecordsmergetype="KEEP_IF_ANY_UNFILTERED",
                              output.dir,
                              run.cmd=TRUE,
                              mc.cores=1){
  
  pon.vcf=file.path(output.dir, paste0("NormalPanel.vcf"))
  
  vcf.list=paste("-V", normal.vcf)
  for(i in 2:length(normal.vcf)) vcf.list=paste(vcf.list, vcf.list[i])
  
  cmd=paste("java", "-jar", GATK.path,
            "-T CombineVariants",
            "-R", ref.fa,
            vcf.list[1],
            "-minN", minN,
            "--setKey null",
            "--filteredAreUncalled",
            "--filteredrecordsmergetype", filteredrecordsmergetype,
            "--genotypemergeoption UNSORTED",
            "-o", pon.vcf)
  
  message("[[",Sys.time(),"]] Run GATK CombineVariants --------")
  print_message(cmd)
  if(run.cmd) mclapply(cmd, system, mc.cores=mc.cores)
  cat(cmd, file=file.path(output.dir, "run.mutect.log"), sep="\n", append=FALSE)
  pon.vcf
}


#' @title mutect2
#' @description A wrapper function to run GATK (MuTect2) Processed through the variant calling as tumor-normal pairs. 
#' @usage run.mutect2(output.dir, ref.fa, tumor.bam, normal.bam, pon.vcf, cosmic.vcf, ref.dbSNP, contamination_fraction_to_filter=0.02, run.cmd=TRUE, mc.cores=1)
#' @param output.dir Output directory
#' @param ref.fa Reference fasta file
#' @param tumor.bam Tumor sample bam files
#' @param normal.bam Bam files form normal samples obtained from a function gatk.mutect.normal
#' @param pon.vcf Panel of normal samples obtained from a function gatk.combinedvariant
#' @param cosmic.vcf Known variant sites of cosmic database vcf file
#' @param ref.dbSNP Known SNP sites reference vcf
#' @param contamination_fraction_to_filter A parameter value for --contamination_fraction_to_filter in GATK MuTect2. Fraction of contamination to aggressively remove (default=0.02)
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details MuTect2 is designed to produce somatic variant calls only, and includes some logic to skip variant sites that are very 
#'          clearly germline based on the evidence present in the Normal sample compared to the Tumor sample. 
#' @return VCF files
#' @references Sensitive detection of somatic point mutations in impure and heterogeneous cancer samples
#' @seealso \url{https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_cancer_m2_MuTect2.php}
#' @export
mutect2=function(output.dir,
                 ref.fa,
                 tumor.bam,
                 normal.bam,
                 pon.vcf,
                 cosmic.vcf,
                 ref.dbSNP,
                 contamination_fraction_to_filter=0.02,
                 run.cmd=TRUE,
                 mc.cores=1){
  
  out.vcf=file.path(output.dir, paste0("MuTect_variants.vcf"))
  
  cmd=paste("java -jar", GATK.path,
            "-T MuTect2",
            "-R", ref.fa,
            "-I:tumor", tumor.bam,
            "-I:normal", normal.bam,
            "--normal_panel", pon.vcf,
            "--cosmic", cosmic.vcf,
            "--dbsnp", ref.dbSNP,
            "--contamination_fraction_to_filter", contamination_fraction_to_filter,
            "-o", out.vcf,
            "--output_mode EMIT_VARIANTS_ONLY",
            "--disable_auto_index_creation_and_locking_when_reading_rods")
  
  message("[[",Sys.time(),"]] Run GATK MuTect2 --------")
  print_message(cmd)
  if(run.cmd) mclapply(cmd, system, mc.cores=mc.cores)
  cat(cmd, file=file.path(output.dir, "run.mutect.log"), sep="\n", append=FALSE)
  out.vcf
}
