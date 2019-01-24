
#' @title read.pileup.gz
#' @description A wrapper function to run samtools (mpileup)
#' @usage read.pileup.gz(ref.fa, fns.bam, sample.name, output.dir, mapQ=1, run.cmd=TRUE, mc.cores=1)
#' @param ref.fa Reference fasta file path
#' @param normal.bam BAM files of normal sample
#' @param tumor.bam BAM files of tumor sample
#' @param sample.name A character vector for the sample names
#' @param mapQ Mapping quality (default=1)
#' @param run.cmd  Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details Generates VCF, BCF or pileup for one or multiple BAM files. Alignment records are grouped by sample (SM) identifiers in @RG header 
#'          lines. If the sample identifiers are absent, each input file is regarded as one sample.
#' @seealso \url {http://www.htslib.org/doc/samtools.html}
#' @export
read.pileup.gz=function(ref.fa, 
                        fns.bam,
                        sample.name, 
                        mapQ=1,
                        output.dir,
                        run.cmd=TRUE, 
                        mc.cores=1){
  
  #output
  pileup.fns=file.path(output.dir, paste0(sample.name, ".pileup.gz"))
  
  #command line
  cmd=paste(samtools.path, "mpileup", "-f", ref.fa, "-Q", mapQ, fns.bam, "|", "gzip", ">", pileup.fns)
  print_message(cmd)
  message("[[",Sys.time(),"]] Run SAMtools pileup and gzip of normal samples ---- ")
  
  if(run.cmd) mclapply(cmd, system, mc.cores=mc.cores)
  cat(cmd, file=file.path(output.dir, "run.samtools.pileup.log"), sep="\n")
  message("[[",Sys.time(),"]] Done ---- ")
  
  #return
  pileup.fns
}



#' @title generate.GC
#' @description A wrapper function to run sequenza-utils.py in sequenza (GC-windows)
#' @usage generate.GC(window=1,000,000, output.dir, ref.fa, run.cmd=TRUE)
#' @param window A parameter value for -w in sequenza. Indicate a window size (in bases), to be used for the binning.
#'               The heterozygous positions and the positions carrying variant calls are not affected by binning.
#' @param output.dir Output directory
#' @param ref.fa Reference fasta file path
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @details Calculation GC contents from reference fasta file
#' @references Sequenza: allele-specific copy number and mutation profiles from tumor sequencing data.
#' @seealso {https://cran.r-project.org/web/packages/sequenza/vignettes/sequenza.pdf}
#' @export
generate.GC=function(window=1000000, 
                     output.dir, 
                     ref.fa, 
                     run.cmd=TRUE){
  
  options("scipen"=10)
  #output
  gc.fn=file.path(output.dir, paste0("gc", window, "Base.txt.gz"))
  
  #command line
  cmd=paste(sequenza.util, "GC-windows -w", window, ref.fa, "|", "gzip", ">", gc.fn)
  
  print_message(cmd)
  message("[[",Sys.time(),"]] Generating a genome-wide GC content file ---- ")
  
  if(run.cmd) system(cmd)
  cat(cmd, file=file.path(output.dir, "run.generateGC.log"), sep="\n")
  message("[[",Sys.time(),"]] Done ---- ")
  #return
  gc.fn
}



#' @title generate.seqz
#' @description A wrapper function to run sequenza-utils.py in sequenza (pileup2seqz, seqz-binning)
#' @usage pileup2seqz(gc.fn, normal.pileup.gz, window=1000000, tumor.pileup.gz, output.dir, run.cmd=T, mc.cores=1)
#' @param fn.gc output file of generate.GC function
#' @param normal.pileup.gz samtools pileup file of normal sample
#' @param tumor.pileup.gz samtools pileup file of tumor sample
#' @param window A parameter value for -w in sequenza. Indicate a window size (in bases), to be used for the binning.
#' @details A seqz file contains genotype information, alleles and mutation frequency, and other features. This file is used as input 
#'          for the R-based part of Sequenza.
#' @references Sequenza: allele-specific copy number and mutation profiles from tumor sequencing data.
#' @seealso {https://cran.r-project.org/web/packages/sequenza/vignettes/sequenza.pdf}

#' @export
pileup2seqz=function(gc.fn, 
                     normal.pileup.gz, 
                     window=1000000, 
                     tumor.pileup.gz, 
                     output.dir, 
                     run.cmd=TRUE, 
                     mc.cores=1){
  
  #normal sample, tumor sample
  pileup.gz.list=list(normal=normal.pileup.gz, tumor=tumor.pileup.gz)
  normal.sample=sub(".pileup.gz", "", basename(normal.pileup.gz))
  tumor.sample=sub(".pileup.gz", "", basename(tumor.pileup.gz))
  
  #output
  out.seqz=file.path(output.dir, paste0(normal.sample,"-",tumor.sample,".seqz.gz"))
  bin.out.seqz=file.path(output.dir, paste0(normal.sample,"-",tumor.sample,".bin.seqz.gz"))
  
  #command line
  cmd=paste(sequenza.util, "pileup2seqz", "-gc", gc.fn, "-n", pileup.gz.list$normal, "-t", pileup.gz.list$tumor, "|", "gzip", ">", out.seqz)
  print_message(cmd)
  message("[[",Sys.time(),"]] Generate a seqz file ---- ")
  
  if(run.cmd) mclapply(cmd, system, mc.cores=mc.cores)
  cat(cmd, file=file.path(output.dir, "run.pileup2seqz.log"), sep="\n", append=TRUE)
  message("[[",Sys.time(),"]] Done ---- ")
  
  # To reduce the size of the seqz file, using binning function
  cmd=paste(sequenza.util, "seqz-binning", "-w", window, "-s", out.seqz, "|", "gzip", ">", bin.out.seqz)
  message("[[",Sys.time(),"]] Reduce the size of the seqz file  ---- ")
  if(run.cmd) mclapply(cmd, system, mc.cores=mc.cores)
  cat(cmd, file=file.path(output.dir, "run.pileup2seqz.log"), sep="\n", append=TRUE)
  message("[[",Sys.time(),"]] Done ---- ")
  
  bin.out.seqz
}




#' @title seqz2rda
#' @description Saves seqz file in R data format.
#' @usage seqz2rda(cnv.dir)
#' @param cnv.dir output directory
#' @export
seqz2rda=function(cnv.dir){
  
  cnv.files=get.fns(input.dir=cnv.dir, idx="bin.seqz.gz$")
  sample.name=sub(".bin.seqz.gz", "", basename(cnv.files))
  
  # read seqz files and save Rdata
  for(i in 1:length(cnv.files)){
    
    chr.name=paste0("chr", c(1:22,"X"))
    
    mclapply(chr.name, function(a) {
      seqz=read.seqz(cnv.files[i], chr.name=a)
      seqz.name=paste0(sub(".gz", "", basename(cnv.files[i])), "_", a)
      assign(seqz.name, seqz)
      fn=file.path(cnv.dir, paste0(seqz.name, ".rda"))
      save(list=seqz.name, file=fn)
    })
  }
  message("[[",Sys.time(),"]] seqz files saved as R data in R object directory ----")
}




#' @title seqz2seg
#' @description Segmentation to estimate DNA copy number variation.
#' @usage seqz2seg(cnv.dir, window=1,000,000)
#' @param cnv.dir Output directory
#' @param window A parameter value for -w in sequenza. Indicate a window size (in bases), to be used for the binning.
#' @details Normalization of depth ratio and DNA segmentation
#' @references Sequenza: allele-specific copy number and mutation profiles from tumor sequencing data.
#' @seealso {https://cran.r-project.org/web/packages/sequenza/vignettes/sequenza.pdf}
#' @export
seqz2seg=function(cnv.dir, window=1000000){
  
  # 1. load seqz data
  options("scipen"=10)
  chr.name = paste0("chr", c(1:22, "X"))
  fns=get.fns(input.dir=cnv.dir, idx=".bin.seqz_chr")
  
  seqzList=mclapply(fns, function(a){
    seqz=get(load(a))
    rm(list = removeExt(basename(a)))
    seqz
  })
  names(seqzList)=sub(".rda", "" , basename(fns))
  
  # 2. segmentation
  #seqz.data_chr = seqzList$chr1
  segList = lapply(seqzList, function(seqz.data_chr){
    # Normalize coverage by GC-content
    gc.stats <- gc.norm(x = seqz.data_chr$depth.ratio, gc = seqz.data_chr$GC.percent)
    gc.vect <- setNames(gc.stats$raw.mean, gc.stats$gc.values)
    
    
    # Correct the coverage of the loaded chromosome:
    seqz.data_chr$adjusted.ratio = seqz.data_chr$depth.ratio / gc.vect[as.character(seqz.data_chr$GC.percent)]
    
    # Select the heterozygous positions
    seqz.het <- seqz.data_chr[seqz.data_chr$zygosity.normal != 'hom', ]
    
    # Detect breakpoints
    breaks <- find.breaks(seqz.het)
    
    # use heterozygous and homozygous position to measure segment values
    seg.s1 <- segment.breaks(seqz.data_chr, breaks = breaks)
    
    # Binning the values of depth ratio and B allele frequency
    seqz.r.win <- windowValues(x = seqz.data_chr$adjusted.ratio,
                               positions = seqz.data_chr$position,
                               chromosomes = seqz.data_chr$chromosome,
                               window = window, overlap = 1,
                               weight = seqz.data_chr$depth.normal)
    seqz.b.win <- windowValues(x = seqz.het$Bf,
                               positions = seqz.het$position,
                               chromosomes = seqz.het$chromosome,
                               window = window, overlap = 1,
                               weight = round(x = seqz.het$good.reads, digits = 0))
    
    list(seg.s1 = seg.s1, seqz.r.win = seqz.r.win, seqz.b.win = seqz.b.win)
  })
  save(segList, file=file.path(cnv.dir, "segmentAllList.rda"))
  message("[[",Sys.time(),"]] segments list saved in cnv directory ----")
  segList
}




#' @title ploidyNcellularity
#' @description Calculate ploidy and cellularity
#' @usage ploidyNcellularity(cnv.dir)
#' @param cnv.dir Output directory
#' @details Calculate ploidy and cellularity for each paired-sample and quantify the copy number
#' @references Sequenza: allele-specific copy number and mutation profiles from tumor sequencing data.
#' @seealso {https://cran.r-project.org/web/packages/sequenza/vignettes/sequenza.pdf}
#' @export
ploidyNcellularity=function(cnv.dir){
  
  segList=get(load(file.path(cnv.dir, "segmentAllList.rda")))
  chr.name=paste0("chr", c(1:22,"X"))
  uniq.sample.name=unique(sub("\\..*", "", names(segList)))
  
  message("[[",Sys.time(),"]] Calculation ploidy and cellularity ----")
  
  for(i in 1:length(segList)){
    segList[[i]]$seqz.r.win[[1]]$seqnames=names(segList[[i]]$seqz.r.win)
    segList[[i]]$seqz.b.win[[1]]$seqnames=names(segList[[i]]$seqz.b.win)
  }
  
  seg.s1.dfList=lapply(uniq.sample.name, function(a) data.frame(do.call(rbind, lapply(segList[paste0(a, ".bin.seqz_", chr.name)], function(b) b$seg.s1)), sample=a,row.names = NULL))
  seqz.baf.dfList=lapply(uniq.sample.name, function(a) data.frame(do.call(rbind, lapply(segList[paste0(a, ".bin.seqz_", chr.name)], function(b) b$seqz.b.win[[1]])), row.names = NULL))
  seqz.ratio.dfList=lapply(uniq.sample.name, function(a) data.frame(do.call(rbind, lapply(segList[paste0(a, ".bin.seqz_", chr.name)], function(b) b$seqz.r.win[[1]])), row.names = NULL))
  
  seg.filtered=lapply(seg.s1.dfList, function(seg) seg[(seg$end.pos - seg$start.pos) > 10e6,])
  CP=lapply(seg.filtered, function(a) baf.model.fit(Bf=a$Bf, depth.ratio = a$depth.ratio,
                                                    avg.depth.ratio = 1, ploidy = seq(0.5,3,0.05)))
  confint=lapply(CP, function(a) get.ci(a))
  ploidy=lapply(confint, function(a) a$max.ploidy)
  cellularity=lapply(confint, function(a) a$max.cellularity)
  
  cn.alleles=list()
  for(i in 1:length(seg.s1.dfList)){
    cn.alleles[[i]]=baf.bayes(Bf=seg.s1.dfList[[i]]$Bf, depth.ratio=seg.s1.dfList[[i]]$depth.ratio,
                              cellularity=cellularity[[i]], ploidy=ploidy[[i]], avg.depth.ratio=1)
    
  }
  
  seg=lapply(1:length(cn.alleles), function(a) cbind(seg.s1.dfList[[a]], cn.alleles[[a]], ploidy[[a]], cellularity[[a]]))
  names(seg)=uniq.sample.name
  save(seg, file=file.path(cnv.dir, "ploidyNcellularity.rda"))
  message("[[",Sys.time(),"]] Done ----")
  seg
}
