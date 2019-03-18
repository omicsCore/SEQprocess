
#' @title vcf2annovar
#' @description A wrapper function to run vcf2annovar.pl in ANNOVAR
#' @usage vcf2annovar(fns.vcf, output.dir, sample.name, format="vcf4", coverage=0, run.cmd=TRUE, mc.cores=1)
#' @param fns.vcf Path to VCF files
#' @param output.dir Output directory
#' @param sample.name A character vector for the sample names
#' @param format A parameter value for -format in ANNOVAR Input files format (.vcf)
#' @param coverage A parameter value for -coverage in ANNOVAR Read coverage threshold in pileup file (default:0)
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details The convert2annovar.pl script convert other "genotype calling" format into ANNOVAR format.
#'          Additionally, the program can generate ANNOVAR input files from a list of dbSNP identifiers, or from transcript identifiers, 
#'          or from a genomic region. 
#' @return Converted annovar files from variant calling format (e.g., .annovar)
#' @references ANNOVAR: functional annotation of genetic variants from high-throughput sequencing data
#' @seealso \url{http://annovar.openbioinformatics.org/en/latest/user-guide/}
#' @export
vcf2annovar = function(fns.vcf,
                       output.dir,
                       sample.name,
                       format = "vcf4",
                       coverage = 0,
                       run.cmd = TRUE,
                       mc.cores=1){
  
  out.dirs=file.path(output.dir, sample.name)
  lapply(1:length(out.dirs), function(a) dir.create(out.dirs[a], recursive = TRUE, showWarnings=FALSE))
  out.fns = file.path(out.dirs, paste0(sample.name, ".annovar"))
  
  message("[[",Sys.time(),"]] Run vcf2annovar -------")
  
  cmd = paste(vcf2annovar.pl, fns.vcf, "-format", format, "-coverage", coverage, "-includeinfo >", out.fns)
  
  print_message(cmd)
  
  if(run.cmd) mclapply(cmd, system,mc.cores=mc.cores)
  cat(cmd, file=file.path(output.dir, "vcf2annovar.run.log"), sep="\n")
  out.fns
}




#' @title table.annovar
#' @description A wrapper function to run table_annovar.pl in ANNOVAR
#' @usage table.annovar(fn.annovar, output.dir, sample.name, annovar.db, ref="hg38", protocol="refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_all,exac03,avsnp147,ljb26_all,cosmic70", protocol.type="g,r,r,f,f,f,f,f,f", nastring=".", run.cmd=TRUE, mc.cores=1)
#' @param fns.annovar Path to annovar files 
#' @param output.dir Output directory
#' @param sample.name A character vector for the sample names
#' @param annovar.db.dir Path to directory with ANNOVAR database
#' @param ref A parameter value for -buildver in ANNOVAR. Specify the genome build version (default: hg38)
#' @param protocol A parameter value for -protocol in ANNOVAR Database names in ANNOVAR (default: "refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_all,exac03,avsnp147,ljb26_all,cosmic70")
#' @param protocol.type A parameter value for -operation in ANNOVAR. Strings separated by commas that specify the types of operation for each protocol (g: genome, r: region, f: filter, default="g,r,r,f,f,f,f,f,f")
#' @param nastring A parameter value for -nastring in ANNOVAR. Strings to display when a score is not available (default: ".")
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details This function takes an input variant file (such as a VCF file) and generate a tab-delimited output file with many columns,
#'          each representing one set of annotations. Additionally, if the input is a VCF file, the program also generates a new output VCF file
#'          with the INFO field filled with annotation information.
#' @return CSV files from annovar format
#' @references ANNOVAR: functional annotation of genetic variants from high-throughput sequencing data
#' @seealso \url{http://annovar.openbioinformatics.org/en/latest/user-guide/}
#' @export
table.annovar= function(fns.annovar,
                        output.dir, 
                        sample.name,
                        annovar.db.dir,
                        ref = "hg38",
                        protocol,
                        protocol.type, 
                        nastring, 
                        run.cmd = TRUE,
                        mc.cores=1){
  
  out.dirs=file.path(output.dir, sample.name)
  lapply(1:length(out.dirs), function(a) dir.create(out.dirs[a], recursive = TRUE, showWarnings=FALSE))
  out.names=file.path(out.dirs, sample.name)
  
  message("[[",Sys.time(),"]] Run table_annovar -------")
  
  cmd=paste(table_annovar.pl, fns.annovar, annovar.db.dir, "--buildver", ref, "-out", out.names, "-remove", 
            "-protocol", protocol, "-operation", protocol.type, "-nastring ", nastring, "-csvout")
  
  print_message(cmd)
  
  if(run.cmd) mclapply(cmd, system,mc.cores=mc.cores)
  cat(cmd, file=file.path(output.dir, "annovar.run.log"), sep="\n", append = FALSE)
  
  out.fns=get.fns(output.dir, idx="multianno.csv$")
  out.fns
}



#' @title Vairant Effect Predictor (VEP)
#' @description A wrapper function to run VEP
#' @usage vep(fns.vcf, output.dir, sample.name, perl5.10.path="/usr/bin", vep.db.dir, run.cmd=TRUE, mc.cores=1)
#' @param fns.vcf Path to VCF files
#' @param output.dir Output directory
#' @param vep.db.dir Specify the cache directory to use.
#' @param perl5.10.path Absolute path to perl 5.10 version. VEP is a Perl based tool. We recommend version 5.10. (default="/usr/bin")
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details The VEP uses the coordinates and alleles in the VCF file to infer biological context for each variant including the location 
#'          of each mutation, its biological consequence (frameshift/ silent mutation), and the affected genes. Variants in the VCF files 
#'          are also matched to known variants from external mutation databases.
#' @return text file and html file included variant information
#' @references The Ensembl Variant Effect Predictor
#' @seealso \url{https://asia.ensembl.org/info/docs/tools/vep/script/vep_options.html}
#' @export
vep=function(fns.vcf, 
             output.dir,
             vep.db.dir,
             sample.name,
             perl5.10.path="/usr/bin",
             
             run.cmd=TRUE,
             mc.cores=1){
  
  out.dirs=file.path(output.dir, sample.name)
  lapply(1:length(out.dirs), function(a) dir.create(out.dirs[a], recursive = TRUE, showWarnings=FALSE))
  out.fns=file.path(out.dirs, paste0(sample.name, ".txt"))
  
  #Variant Effect Predictor
  message("[[", Sys.time(),"]] Run Variant Effect Predictor -------")
  
  cmd=paste(vep.path, "--offline", "-i", fns.vcf, "-o", out.fns, "--dir_cache", vep.db.dir, "--polyphen b --humdiv --sift b --biotype --force_overwrite")
  
  print_message(cmd)
  
  #Temporary set system pass
  path.o=Sys.getenv("PATH")
  Sys.setenv("PATH"=paste(perl5.10.path, path.o, sep=":"))
  
  if(run.cmd) mclapply(cmd, system,mc.cores=mc.cores)
  cat(cmd, file=file.path(output.dir, "vep.run.log"), sep="\n", append = FALSE)
  
  Sys.setenv("PATH"=path.o)
  out.fns
}

