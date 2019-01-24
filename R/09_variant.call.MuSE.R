

#' @title MuSE.call
#' @description A wrapper function to run MuSE (call)
#' @usage MuSE.call(tumor.bam, normal.bam, output.dir, sample.name, ref.fa, run.cmd=TRUE, mc.cores=1)
#' @param normal.bam path to normal sample recalibration bam files
#' @param tumor.bam path to tumor sample recalibration bam files as tumor-normal pair
#' @param output.dir Output directory
#' @param sample.name A character vector for the sample names
#' @param ref.fa Reference fasta file
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @import parallel
#' @details The first step of MuSE, MuSE.call takes as input indexed reference fasta file and BAM files. The BAM files require aligning all the 
#'          sequence reads against the reference genome using the Burrows-Wheeler alignment tool BWA-mem algorithm. 
#'          In addition, the BAM files need to be processed by following the GATK-MarkDuplicates, realigning the paired tumor-normal BAMs 
#'          jointly and recalibrating base quality scores.
#' @return MuSE.call output txt file. 
#' @references MuSE: accounting for tumor heterogeneity using a sample-specific error model improves sensitivity and specificity in mutation calling from sequencing data
#' @seealso \url{http://bioinformatics.mdanderson.org/main/MuSE}
#' @export
MuSE.call=function(tumor.bam,
                   normal.bam,
                   output.dir,
                   sample.name,
                   ref.fa,
                   run.cmd=TRUE,
                   mc.cores=1){
  
  out.dirs=file.path(output.dir, sample.name)
  lapply(1:length(out.dirs), function(a) dir.create(out.dirs[a], recursive = TRUE, showWarnings=FALSE))
  out.prefix=file.path(out.dirs,sample.name)
  out.txt=file.path(out.dirs,paste0(sample.name,".MuSE.txt"))

  message("[[", Sys.time(),"]] Run MuSE call -------")
  cmd=paste(MuSE.path, "call -O", out.prefix, "-f", ref.fa, tumor.bam, normal.bam)
  print_message(cmd)
  
  if(run.cmd) mclapply(cmd, system,mc.cores=mc.cores)
  cat(cmd, file=file.path(output.dir, "MuSE.call.run.log"), sep="\n", append = FALSE)
  
  out.txt
}



#' @title MuSE.sump
#' @description A wrapper function to run MuSE (sump)
#' @usage MuSE.sump(MuSE.txt, output.dir, ref.dbSNP, ref.gold_indels, data.type="E", run.cmd=TRUE, mc.cores=1)
#' @param MuSE.txt Path to MuSE.call output text file
#' @param output.dir Output directory
#' @param MuSE.data.type E is used for WXS data and G can be used for WGS data
#' @param ref.dbSNP Known SNP sites reference
#' @param ref.gold_indels Known Indel sites reference
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details The second step, MuSE sump, takes as input the output file from MuSE.call and dbSNP variant call format file. MuSE provide 
#'          two options for building the sample-specific error model. One is applicable to WES data (option ‘-E’), and the other to WGS 
#'          data (option -G).
#' @return vcf files included variant call information
#' @references MuSE: accounting for tumor heterogeneity using a sample-specific error model improves sensitivity and specificity in mutation calling from sequencing data
#' @seealso \url{http://bioinformatics.mdanderson.org/main/MuSE}
#' @export
MuSE.sump=function(MuSE.txt,
                   output.dir,
                   MuSE.data.type=c("E", "G"),
                   ref.dbSNP,
                   ref.gold_indels,
                   
                   run.cmd=TRUE,
                   mc.cores=1){
  sample.name=sub(".MuSE.txt", "", basename(MuSE.txt))
  out.dirs=file.path(output.dir, sample.name)
  lapply(1:length(out.dirs), function(a) dir.create(out.dirs[a], recursive = TRUE, showWarnings=FALSE))
  out.fns=file.path(out.dirs, paste0(sub(".MuSE.txt","", basename(MuSE.txt)), ".muse_variants.vcf"))
  
  MuSE.data.type=paste0("-",MuSE.data.type)
  
  
  message("[[", Sys.time(),"]] Run MuSE sump -------")
  
  # check ref
  ref.dbSNP.gz=paste0(ref.dbSNP, ".gz")
  ref.dbSNP.tbi=paste0(ref.dbSNP.gz, ".tbi")
  if(!file.exists(ref.dbSNP.gz)) vcf2gz(ref.dbSNP)
  if(!file.exists(ref.dbSNP.tbi)) tabix.vcf(ref.dbSNP.gz)
  
  ref.gold_indels.gz=paste0(ref.gold_indels, ".gz")
  ref.gold_indels.tbi=paste0(ref.gold_indels.gz, ".tbi")
  if(!file.exists(ref.gold_indels.gz)) vcf2gz(ref.gold_indels)
  if(!file.exists(ref.gold_indels.tbi)) tabix.vcf(ref.gold_indels.gz)
  
  cmd=paste(MuSE.path, "sump -I", MuSE.txt, "-O", out.fns, MuSE.data.type, "-D", ref.dbSNP.gz, ref.gold_indels.gz)
  print_message(cmd)
  
  if(run.cmd) mclapply(cmd, system,mc.cores=mc.cores)
  
  cat(cmd, file=file.path(output.dir, "MuSE.sump.run.log"), sep="\n", append = FALSE)
  
  out.fns
}



#' @title vcf2gz
#' @description Compress VCF file to gz file
#' @param vcf Path to dbSNP, indel reference vcf file
#' @export
vcf2gz=function(vcf){
  vcf.gz=paste0(vcf,".gz")
  cmd=paste("bgzip -c", vcf, ">", vcf.gz )
  system(cmd)
  vcf.gz
}


#' @name Muse.tabix2.vcf
#' @title VCFtabix
#' @description Tabix indexes a TAB-delimited genome position file in.tab.bgz.
#' @param vcf.gz.file Path to dbSNP, indel reference vcf.gz file
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @export
tabix.vcf=function(vcf.gz.file, run.cmd=TRUE){
  cmd=paste("tabix -p vcf", vcf.gz.file)
  
  if(run.cmd) system(cmd)
  out.idx=paste0(file.path(vcf.gz.file),".tbi")
  out.idx
}
