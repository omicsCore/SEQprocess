
#' @title somaticsniper
#' @description A wrapper function to run SomaticSniper
#' @usage somaticsniper(ref.fa, tumor.bam, normal.bam, output.dir, sample.name, mapQual=1, LOH=TRUE, Genotype=TRUE, somaticQual=15, somaticMutation=0.01, Theta=0.85, Hap.number=2, Hap.diff=0.001, out.format="vcf", run.cmd=TRUE, mc.cores=1)
#' @param ref.fa Reference fasta file
#' @param tumor.bam Tumor sample bam files
#' @param normal.bam Normal sample bam files
#' @param output.dir Output directory
#' @param sample.name A character vector for the sample names
#' @param mapQual A parameter value for -q in SomaticSniper Filtering reads with mapping quality less than INT (default:1)
#' @param LOH A parameter value for -L in SomaticSniper. Do not report LOH variants as determined by genotypes (logical)
#' @param Genotype A parameter value for -G in SomaticSniper. Do not report Gain of Referene variants as determined by genotypes (logical)
#' @param somaticQual A parameter value for -Q in SomaticSniper. Filtering somatic SNV output with somatic quality less than INT (default:15)
#' @param somaticMutation A parameter value for -s in SomaticSniper. Prior probability of a somatic mutation (default:0.01)
#' @param Theta A parameter value for -T in SomaticSniper. Theta in maq consensus calling model (default:0.85)
#' @param Hap.number A parameter value for -N in SomaticSniper. Number of haplotypes in the sample (default:2)
#' @param Hap.diff A parameter value for -r in SomaticSniper. Prior of a difference between two haplotypes (default:0.001)
#' @param out.format A parameter value for -F in SomaticSniper. Select output format (vcf or classic) (default:vcf)
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details The purpose of this program is to identify single nucleotide positions that are different between tumor and normal
#'          (or in theory, any two bam files). It takes a tumor bam and a normal bam and compares the two to determine the differences. 
#' @return VCF files        
#' @references SomaticSniper: identification of somatic point mutations in whole genome sequencing data.
#' @seealso \url{http://gmt.genome.wustl.edu/packages/somatic-sniper/}
#' @export
somaticsniper=function(tumor.bam,
                       normal.bam,
                       output.dir,
                       sample.name=sub(bam.idx, "", basename(tumor.bam)),
                       
                       ref.fa,
                       
                       mapQual= 1,
                       LOH= TRUE, #"-L", 
                       Genotype= TRUE, #"-G",
                       somaticQual= 15,
                       somaticMutation= 0.01,
                       Theta= 0.85,
                       Hap.number= 2,
                       Hap.diff= 0.001,
                       out.format= "vcf",
                       
                       run.cmd=TRUE,
                       mc.cores=1){
  
  out.dirs=file.path(output.dir, sample.name)
  lapply(1:length(out.dirs), function(a) dir.create(out.dirs[a], recursive = TRUE, showWarnings=FALSE))
  out.fns= file.path(out.dirs, paste0(sample.name, "_variants.vcf"))
  
  #SomaticSniper
  message("[[", Sys.time(),"]] Run SomaticSniper -------")
  
  LOH=ifelse(LOH, "-L", "")
  Genotype=ifelse(Genotype, "-G", "")
  
  cmd=paste(somaticsniper.path, "-q", mapQual, LOH, Genotype, "-Q", somaticQual, "-s", somaticMutation, "-T", Theta, "-N", Hap.number, "-r", Hap.diff, "-F", out.format, "-f", ref.fa, tumor.bam, normal.bam, out.fns)
  
  print_message(cmd)
  
  if(run.cmd) mclapply(cmd, system,mc.cores=mc.cores)
  cat(cmd, file=file.path(output.dir, "somaticsniper.run.log"), sep="\n", append = FALSE)
  
  out.fns
}

  
