
#' @title cufflinks
#' @description A wrapper function to run Cufflinks for mRNA quantitation
#' @usage cufflinks(fns.bam, sample.name, output.dir, cufflinks_thread_number=4, cufflinks.gtf="G", ref.gtf, run.cmd=TRUE, mc.cores=1)
#' @param fns.bam Path to bam files
#' @param sample.name A character vector for the sample names
#' @param output.dir Output directory
#' @param cufflinks.thread.number A parameter value for -p in Cufflinks. A numeric value of the number of threads (default: 4)
#' @param cufflinks.gtf  If you set -G, Output will not include novel genes and isoforms that are assembled. (default: -G)
#' @param ref.gtf Path to reference gtf file
#' @param run.cmd  Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details Cufflinks algorithms for transcript assembly and expression quantification are much more accurate with paired-end reads. 
#' @return  mRNA quantification text files 
#' @import parallel
#' @references http://cole-trapnell-lab.github.io/cufflinks/papers/
#' @seealso \url{http://cole-trapnell-lab.github.io/cufflinks/}
#' @export
cufflinks=function(fns.bam, 
                   output.dir,
                   sample.name, 
                   cufflinks.thread.number=4, 
                   cufflinks.gtf=c("G", "g"), 
                   ref.gtf, 
                   run.cmd=TRUE, 
                   mc.cores=1){
  
  out.dirs=file.path(output.dir, sample.name)
  lapply(1:length(out.dirs), function(a) dir.create(out.dirs[a], recursive = TRUE, showWarnings=FALSE))
  sapply(out.dirs, dir.create,  showWarnings = FALSE)  
  
  cufflinks.gtf=paste0("-", cufflinks.gtf)
  cmd=paste(cufflinks.path, "-p", cufflinks.thread.number, "-o", out.dirs, cufflinks.gtf, ref.gtf, "-multi-read-correct", fns.bam)
  
  message("[[",Sys.time(),"]] Run cufflinks -------")
  print_message(cmd)
  
  if(run.cmd) mclapply(cmd, system, mc.cores=mc.cores)
  print_message(cmd)
  cat(cmd, file=file.path(output.dir, "run.cufflinks.log"), sep = "\n",append = FALSE)
  
  out.fns=list.files(out.dirs,".fpkm_tracking", full.names=TRUE)
  out.fns
}



#' @title htseq_count
#' @description A wrapper function to run htseq-count for mRNA or miRNA quantitation
#' @usage htseq_count(RNAtype="mRNA", fns.bam, sample.name, output.dir, Mode="intersection-nonempty", stranded="no", idattr="gene_id", htseq.r="pos", htseq.a=10, ref.gtf, mir.gff, run.cmd=TRUE, mc.cores=1)
#' @param fns.bam Path to input BAM or SAM files
#' @param sample.name A character vector for the sample names
#' @param output.dir Output directory
#' @param ref.gtf  Directoy stored at reference gtf file
#' @param mir.gff  Directoy stored at micro-RNA reference gff file
#' @param MODE A parameter value for -m in htseq-count. Mode to handle reads overlapping more than one feature (default:intersection-nonempty)
#' @param stranded A parameter value for -s in htseq-count. Whether the data is from a strand-specific assay (default:no)
#' @param idattr A parameter value for -i in htseq-count. GFF attribute to be used as feature ID (default:"gene_id")
#' @param htseq.r A parameter value for -r in htseq-count. Sorting order method (default:"pos")
#' @param htseq.a A parameter value for -a in htseq-count. Skip all reads with alignment quality lower than the given minimum value (default: 10)
#' @param run.cmd  Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details Counting reads in features. Given a file with aligned sequencing reads and a list of genomic features, a common task is to 
#'          count how many reads map to each feature.
#' @return Text file included read count information
#' @import parallel
#' @references HTSeq—a Python framework to work with high-throughput sequencing data
#' @seealso \url {https://htseq.readthedocs.io/en/release_0.9.1/}
#' @export
htseq_count=function(RNAtype=c("mRNA", "miRNA"),
                     fns.bam,
                     sample.name,
                     output.dir,
                     
                     #option
                     Mode=c("intersection-nonempty","union","intersection-strict"),
                     stranded=c("no","yes"),
                     idattr="gene_id",
                     htseq.r=c("pos","name"),
                     htseq.a=10,
                     
                     #ref
                     ref.gtf,
                     mir.gff,
                     
                     run.cmd = TRUE,
                     mc.cores = 1
){
  out.dirs=file.path(output.dir, sample.name)
  lapply(1:length(out.dirs), function(a) dir.create(out.dirs[a], recursive = TRUE, showWarnings=FALSE))
  out.fns=file.path(out.dirs,paste0(sample.name, ".count.txt"))
  
  message("[[",Sys.time(),"]] Run htseq-count----")
  
  if(RNAtype=="mRNA") {
    cmd=paste(samtools.path, "view -F 4", fns.bam, "|", htseq.path, "-m", Mode, "-i", idattr, "-r", htseq.r, "-s", stranded, "-", ref.gtf, ">", out.fns)
    print_message(cmd)
    if(run.cmd) mclapply(cmd, system, mc.cores = mc.cores)
    cat(cmd, file=file.path(output.dir, "htseq_count.run.log"), sep="\n", append=FALSE)
  }
  
  if(RNAtype=="miRNA") {
    
    cmd=paste(htseq.path, "-t miRNA -i Name -a", htseq.a, fns.bam, mir.gff, ">", out.fns)
    print_message(cmd)
    if(run.cmd) mclapply(cmd, system, mc.cores = mc.cores)
    cat(cmd, file=file.path(output.dir, "htseq_count.run.log"), sep="\n", append=FALSE)
  }
  out.fns
}



#' @title gtf2gr
#' @description Converts reference gtf file to GRanges form to execute FPKM estimation
#' @usage gtf2gr(ref.gtf, output.dir)
#' @param ref.gtf Directoy stored at reference gtf file (e.g. gencode.v22.gtf)
#' @param output.dir Output directory
#' @details To normalize the number of reads of each feature calculated in the previous step to the value of FPKM, convert the reference 
#'          gtf file to GRanges format.
#' @import limma
#' @import GenomicRanges
#' @export
gtf2gr=function(ref.gtf, output.dir){
    gtf=read.delim(ref.gtf, comment.char="#", header=FALSE, sep="\t")
    colnames(gtf)=c("seqnames","source", "feature", "start", "end", "score", "strand", "frame", "annot")
    
    # gene
    gtf.gene=gtf[gtf$feature == "gene",]
    annot=strsplit2(as.character(gtf.gene$annot), split = " |;")
    gtf.gene$gene_id=annot[,2]
    gtf.gene$gene_type=annot[,5]
    gtf.gene$gene_name=annot[,8]
    gtf.gene=gtf.gene[,-9]
    
    cols=c("seqnames", "start", "end", "strand")
    message("[[",Sys.time(),"]] Convert GTF file to GRanges form---- ")
    gtf.gr=GRanges(seqnames = gtf.gene$seqname, 
                   ranges = IRanges(gtf.gene$start, gtf.gene$end), 
                   strand = gtf.gene$strand, 
                   gtf.gene[,setdiff(colnames(gtf.gene), cols)])
    
    
    save(gtf.gr, file = file.path(output.dir, "gtf.gr.rda"))
    message("[[",Sys.time(),"]] GRanges is stored in the R object directory---- ")
    file.path(output.dir, "gtf.gr.rda")
}




#' @title gff2gr
#' @description Converts reference gff file to GRanges form 
#' @usage gff2gr(mir.gff, output.dir)
#' @param mir.gff Directoy stored at reference gff file 
#' @param output.dir Output directory
#' @import limma
#' @import GenomicRanges
#' @export
gff2gr=function(mir.gff, output.dir){
  gff=read.delim(mir.gff, comment.char="#", header=FALSE, sep="\t")
  colnames(gff)=c("seqnames","source", "feature", "start", "end", "score", "strand", "frame", "annot")
  
  # gene
  gff.gene=gff[gff$feature == "miRNA",]
  annot=strsplit2(as.character(gff.gene$annot), split = "=|;")
  gff.gene$ID=annot[,2]
  gff.gene$Alias=annot[,4]
  gff.gene$Name=annot[,6]
  gff.gene$Derives_from=annot[,8]
  gff.gene=gff.gene[,-9]
  
  cols=c("seqnames", "start", "end", "strand")
  message("[[",Sys.time(),"]] Convert GFF file to GRanges form---- ")
  gff.gr=GRanges(seqnames = gff.gene$seqname, 
                 ranges = IRanges(gff.gene$start, gff.gene$end), 
                 strand = gff.gene$strand, 
                 gff.gene[,setdiff(colnames(gff.gene), cols)])
  
  
  save(gff.gr, file = file.path(output.dir, "gff.gr.rda"))
  message("[[",Sys.time(),"]] GRanges is stored in the R object directory---- ")
  file.path(output.dir, "gff.gr.rda")
}



#' @title htseq.add.info
#' @description Add information to the htseq output file
#' @usage htseq.add.info(RNAtype="mRNA", count.fns, output.dir, output.dir, mc.cores=1)
#' @param RNAtype RNAtype (default="mRNA")
#' @param fns.count count file paths
#' @param output.dir Directory stored at FPKM conunt files
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details Adds information to the output file of htseq. (Gene name, chromosome, start position, end position, gene size, FPKM value)
#' @export
htseq.add.info=function(RNAtype=c("mRNA", "miRNA"),
                        count.fns,
                        output.dir,
                        mc.cores=1){
  
  sample.name=sub(".count.txt", "", basename(count.fns))
  out.dirs=file.path(output.dir, sample.name)
  lapply(1:length(out.dirs), function(a) dir.create(out.dirs[a], recursive = TRUE, showWarnings=FALSE))
  
  if(RNAtype=="mRNA"){
  # read counts
  gtf.gr=get(load(file.path(output.dir, "gtf.gr.rda")))
  countList=mclapply(count.fns, read.delim, header=FALSE, mc.cores=mc.cores)
  info.fns=file.path(out.dirs, sub(".count.txt", ".CountAddInfo.mRNA.txt", basename(count.fns)))
  # adding information
  countList=mclapply(countList, function(count){
    dat=count[1:(nrow(count)-5),]
    colnames(dat) = c("id", "counts")
    dat$gene_name=gtf.gr[match(as.character(dat[,1]), gtf.gr$gene_id)]$gene_name
    dat$seqnames=as.character(seqnames(gtf.gr[match(as.character(dat[,1]), gtf.gr$gene_id)]))
    dat$start=start(gtf.gr[match(as.character(dat[,1]), gtf.gr$gene_id)])
    dat$end=end(gtf.gr[match(as.character(dat[,1]), gtf.gr$gene_id)])
    dat$size=width(gtf.gr[match(as.character(dat[,1]), gtf.gr$gene_id)])
    dat$FPKM=as.numeric(c(dat$counts/dat$size/sum(dat$counts)*10^9))
    dat
  }, mc.cores=mc.cores)
  }
  
  
  if(RNAtype=="miRNA"){
  # read counts  
  gff.gr=get(load(file.path(output.dir, "gff.gr.rda")))
  countList=mclapply(count.fns, read.delim, header=FALSE, mc.cores=mc.cores)  
  info.fns=file.path(out.dirs, sub(".count.txt", ".CountAddInfo.miRNA.txt", basename(count.fns)))
  # adding information  
  countList=mclapply(countList, function(count){
    dat=count[1:(nrow(count)-5),]
    colnames(dat) = c("id", "counts")
    dat$Name=gff.gr[match(as.character(dat[,1]), gff.gr$Name)]$Name
    dat$seqnames=as.character(seqnames(gff.gr[match(as.character(dat[,1]), gff.gr$Name)]))
    dat$start=start(gff.gr[match(as.character(dat[,1]), gff.gr$Name)])
    dat$end=end(gff.gr[match(as.character(dat[,1]), gff.gr$Name)])
    dat$size=width(gff.gr[match(as.character(dat[,1]), gff.gr$Name)])
    dat
  }, mc.cores=mc.cores)
  }
  
  # write output
  for(i in 1:length(countList)) write.table(countList[[i]], file=info.fns[i], row.names=FALSE, sep="\t", quote=FALSE, append=FALSE)
  info.fns
}


