#' @title report
#' @description Reports the result of using SEQprocess()
#' @usage report(envList)
#' @param envList R environment list
#' @details Provides a report that summarizes the processing steps and visualized tables and plots for the processed results. The report file is automatically generated recording the workflows of the data processing steps, the options used in the processing, and the outcome results. 
#' @return pdf file include data processing result information
#' @import rmarkdown
#' @export 
report=function(envList){
  input.param.names=c("project.name", "fastq.dir", "output.dir", "mc.cores", "type","pipeline","qc","trim.method", "align.method", "bwa.method","rm.dup","realign","rseq.abundance.method", "cufflinks.gtf", "variant.call.method", "annotation.method", "ref","make.eSet")
  message("===Input Parameters===")
  
  input.params=envList[input.param.names]
  input.params=data.frame(name=names(input.params), parameter=as.character(input.params), row.names = NULL)
  print(input.params)
  analysis.pipeline= paste(toupper(envList[["type"]]),envList[["pipeline"]],sep = "-")
  
  message("===Output directory===")
  output.dir=envList$output.dir
  print(output.dir)
  
  message("===Result===")
  process.names=envList$proc.names
  print(process.names)
  
  # get report 
  proc.reports=lapply(process.names, function(proc) get.proc.report(proc, output.dir))
  names(proc.reports)=process.names
  print(proc.reports)
  
  #Number of samples
  Sample_ID=sub(fq1.idx,"",basename(unlist(envList[["input.fns"]]["fq1"])))
  
  # summary
  report.summary=c(list(input.params=input.params, process.names=process.names, analysis.pipeline=analysis.pipeline, Sample_ID=Sample_ID), proc.reports)
  save(report.summary, file = system.file("data/report.summary.rda",package = "SEQprocess"))
  library(rmarkdown)
  report.rmd.fn=system.file("extdata/report.rmd/SEQprocess_Report.Rmd",package = "SEQprocess")
  
  
  # set param
  qc.eval=is.element("qc", process.names)
  trim.eval=sum(grepl("trim|cutadapt", process.names))>0
  align.eval=sum(grepl("tophat2|bwa-mem|bwa-aln|bowtie2|star", process.names))>0
  rmdu.eval=sum(grepl("rmdu|rmdu_b", process.names))>0
  realign.eval=is.element("realign", process.names)
  annovar.eval=ifelse(length(report.summary$annot)!=0,sum(grepl("annot|vep", process.names))>0, FALSE)
  eset.eval=is.element("make.Set",process.names)
  
  
  preprocessing.eval=trim.eval&align.eval&rmdu.eval
  trimming.eval= if(!preprocessing.eval) trimming.eval=trim.eval else trimming.eval=FALSE
  alignment.eval= if(!preprocessing.eval) alignment.eval=align.eval else alignment.eval=FALSE
  rmdup.eval= if(!preprocessing.eval) rmdup.eval=rmdu.eval else rmdup.eval=FALSE
  if(eset.eval){
    eSet.eval=report.summary$make.Set$eSet.eval
    vset.eval=report.summary$make.Set$vset.eval
    cset.eval=report.summary$make.Set$cset.eval
  }else {
    eSet.eval=FALSE
    vset.eval=FALSE
    cset.eval=FALSE
  }
  
  render(report.rmd.fn, output_format ="pdf_document", output_file="SEQprocess_Report.pdf", output_dir=output.dir, params =  list(qc.eval=qc.eval, preprocessing.eval=preprocessing.eval, trimming.eval=trimming.eval, alignment.eval=alignment.eval, rmdup.eval=rmdup.eval, realign.eval=realign.eval,annovar.eval=annovar.eval,eset.eval=eset.eval ,eSet.eval=eSet.eval,vset.eval=vset.eval,cset.eval=cset.eval))
}

#' @title get.proc.report 
#' @description Writes report accoding to process
#' @usage get.proc.report(proc=c("qc", "trim","cutadapt", "bwa-mem", "bwa-aln", "tophat2","star","bowtie2","rmdu","rmdu_b", "realign", "gatk", "varscan2","mutect","muse","somaticsniper","annot","vep", "cufflinks","htseq","make.eSet"),output.dir)
#' @param proc Process name
#' @param output.dir Output directory
#' @details Creates report information according to processing result
#' @export 
get.proc.report=function(proc=c("qc", "trim","cutadapt", "bwa-mem", "bwa-aln", "tophat2","star","bowtie2","rmdu","rmdu_b", "realign", "gatk", "varscan2","mutect","muse","somaticsniper","annot","vep", "cufflinks","htseq","make.eSet"),output.dir){
  
  if(proc=="qc"){
    qc.dir=file.path(output.dir, "00_qc")
    qc.report=get.qc.report(qc.dir)
    
    report=qc.report
  }else if(proc=="trim"){
    trim.dir=file.path(output.dir,"01_trim")
    
    trim.report=get.trim.gal.report(trim.dir)       
    
    report=trim.report
  }else if(proc=="cutadapt"){
    
    trim.dir=file.path(output.dir, "01_trim")
    trim.report=get.trim_cut.report(trim.dir)       
    
    report=trim.report
  }else if(proc=="bwa-mem"){
    
    align.dir=file.path(output.dir, "02_align")
    align.report=get.collectmetrics.report(align.dir)
    
    report=align.report
  }else if(proc=="bwa-aln"){
    align.dir=file.path(output.dir, "02_align")
    align.report=get.single_end.metrics.report(align.dir)
    
    report=align.report
  }else if(proc=="tophat2"){
    tophat.dir=file.path(output.dir, "02_align")
    tophat.report=get.tophat.report(tophat.dir)
    
    report=tophat.report
    
  }else if(proc=="bowtie2"){       
    align.dir=file.path(output.dir, "02_align")
    align.report=get.collectmetrics.report(align.dir)
    
    report=align.report
    
  }else if(proc=="star"){
    align.dir=file.path(output.dir, "02_align")
    align.report=get.collectmetrics.report(align.dir)
    
    report=align.report
    
  }else if(proc=="rmdu"|proc=="rmdu_b"){
    rmdup.dir=file.path(output.dir, "03_rmdup")
    rmdup.report=get.collectmetrics.report(rmdup.dir)
    
    report=rmdup.report
    
  }else if(proc=="realign"){       
    realign.dir=file.path(output.dir, "04_realign")
    realign.report=get.collectmetrics.report(realign.dir)
    
    report=realign.report
    
  }else if(proc=="gatk"|proc=="varscan2"|proc=="mutect2"|proc=="muse"|proc=="somaticsniper"){ 
    report=NULL
    
  }else if(proc=="annot"){
    
    annot.dir=file.path(output.dir, "06_annot")
    annot.report=get.annovar.report(annot.dir)
    report=annot.report
  }else if(proc=="vep"){
    report=NULL
  }else if(proc=="cufflinks"|proc=="htseq"){
    
    report=NULL
    
  }else if(proc=="make.Set"){
    
    Robject.dir=file.path(output.dir, "09_Robject")
    Robject.report=get.Robject.report(Robject.dir)
    report=Robject.report
  }
  report
}
#' @title get.qc.report
#' @description Creates a data frame using a fastq file
#' @usage get.qc.report(qc.dir)
#' @param qc.dir Path to directory with fastQC output files
#' @details Adds information (ex.Q30) using the data frame combined with the fastqc file and output the result as a data frame
#' @return data frame of result summary
#' @export
get.qc.report=function(qc.dir){
  suppressPackageStartupMessages(library(fastqcr)) 
  fns= dir(qc.dir, "fastqc.zip$", full.names = TRUE)
  
  res.summary=data.frame(qc_stats(qc_aggregate(qc.dir, progressbar = FALSE)))[,-2]
  
  qual.scores=lapply(fns, function(a) data.frame(qc_read(a, modules = "per sequence quality scores", verbose = FALSE)))
  print(qual.scores)
  
  res.summary$Q30.perc=sapply(qual.scores, function(b) sum(b[which(b[,1]>30),2])/sum(b[,2])*100)
  
  colnames(res.summary)=c("Filename","GC%","Total Sequence","Sequence length","Phred Score(>30)(%)")
  
  list(report=res.summary)
}

#' @title get.trim.gal.report
#' @description Creates a data frame using the trimming output
#' @usage get.trim.gal.report(trim.dir)
#' @param trim.dir Path to directory with trimmed fastq files
#' @details Creates a data frame using the txt file, the output of Trim galore.
#' @return data frame of result summary
#' @export
get.trim.gal.report=function(trim.dir){
  
  fns=dir(trim.dir, "trimming_report.txt$", recursive = TRUE, full.names = TRUE)
  sample.names=sub("_trimming_report.txt", "", basename(fns))
  
  library(limma)
  trimList=lapply(fns, function(a) strsplit2(as.character(read.delim(a)[19:24,]), ":"))
  res.summary=t(sapply(trimList, function(a) gsub(" ", "", a[,2])))
  
  colnames(res.summary)=trimList[[1]][,1]
  report=data.frame(sample.names=sample.names, res.summary)
  report
}

#' @title get.trim_cut.report
#' @description Creates a data frame using the trimming output
#' @usage get.trim_cut.report(trim.dir)
#' @param trim.dir Path to directory with trimmed fastq files
#' @details Creates a data frame using the txt file, the output of Cutadapt.
#' @return data frame of result summary
#' @export
get.trim_cut.report=function(trim.dir){
  
  fns=dir(trim.dir, ".cutadapt.txt$", recursive = TRUE, full.names = TRUE)
  sample.names=sub(".cutadapt.txt$", "", basename(fns))
  
  library(limma)
  trimList=lapply(fns, function(a) strsplit2(as.character(read.delim(a)[5:10,]), ":"))
  res.summary=t(sapply(trimList, function(a) gsub(" ", "", a[,2])))
  
  colnames(res.summary)=trimList[[1]][,1]
  report=data.frame(sample.names=sample.names, res.summary)
  report
  
}

#' @title get.tophat.report
#' @description Creates a data frame using the tophat align summary
#' @usage get.tophat.report(tophat.dir)
#' @param tophat.dir Path to directory with bam files
#' @details Creates a data frame using the txt file, the output of Tophat.
#' @return data frame of result summary
#' @export 
get.tophat.report=function(tophat.dir){
  fns=dir(tophat.dir, "align_summary", recursive = TRUE, full.names = TRUE)
  
  summaryList=lapply(fns, function(a) as.character(read.delim(a, header=FALSE)$V1))
  
  #parser
  total.read=as.numeric(sapply(summaryList, function(a) sub(".* ", "", a[2]))) + as.numeric(sapply(summaryList, function(a) sub(".* ", "", a[6])))
  
  mapping.reads.r1=as.numeric(sapply(summaryList, function(a) sub(".* ", "", sub(" \\(.*", "", a[3]))))
  mapping.reads.r2=as.numeric(sapply(summaryList, function(a) sub(".* ", "", sub(" \\(.*", "", a[7]))))
  
  multi.map.reads.r1=as.numeric(sapply(summaryList, function(a) sub(".* ", "", sub(" \\(.*", "", a[4]))))
  multi.map.reads.r2=as.numeric(sapply(summaryList, function(a) sub(".* ", "", sub(" \\(.*", "", a[8]))))
  
  paired.map.reads=as.numeric(sapply(summaryList, function(a) sub(".* ", "", a[10])))
  
  sample.names=basename(dirname(fns))
  
  report=data.frame(sample.names=sample.names, total.read=total.read, mapping.reads.r1=mapping.reads.r1, mapping.reads.r2=mapping.reads.r2, multi.map.reads.r1=multi.map.reads.r1, multi.map.reads.r2=multi.map.reads.r2, paired.map.reads=paired.map.reads)
  report[is.na(report)]="."
  report  
}

#' @title get.collectmetrics.report
#' @description Creates a data frame using picard.collectmultiplemetrics()
#' @usage get.collectmetrics.report(bam.dir, collectmetrics.idx=".alignment_summary_metrics$")
#' @param bam.dir Path to directory with bam files
#' @param collectmetric.idx Index of files (default=".alignment_summary_metrics$")
#' @details Creates a data frame using the txt file, the output of  picard.collectmultiplemetrics().
#' @return data frame of result summary
#' @export
get.collectmetrics.report=function(bam.dir, collectmetrics.idx=".alignment_summary_metrics$", mc.cores=1){
  
  fns=dir(bam.dir, collectmetrics.idx, recursive = TRUE, full.names = TRUE)
  fns.bam=list.files(bam.dir, bam.idx, recursive = TRUE, full.names = TRUE)
  
  if(length(fns)==0) fns=picard.collectmetrics(fns.bam = fns.bam, ref.fa = ref.fa, mc.cores=mc.cores)
  
  statList=lapply(fns, function(a) read.delim(a, header=FALSE,comment.char = "#"))
  Total.reads=as.numeric(sapply(statList, function(a)  as.character(a$V2[4])))
  PF.reads=as.numeric(sapply(statList, function(a)  as.character(a$V3[4])))
  Mapped.reads.R1=as.numeric(sapply(statList, function(a)  as.character(a$V6[2])))
  Mapped.reads.R1.pct=as.numeric(sapply(statList, function(a)  as.character(a$V7[2])))
  Mapped.reads.R2=as.numeric(sapply(statList, function(a)  as.character(a$V6[3])))
  Mapped.reads.R2.pct=as.numeric(sapply(statList, function(a)  as.character(a$V7[3])))
  Paired.Mapped=as.numeric(sapply(statList, function(a)  as.character(a$V17[4])))
  Paired.Mapped.pct=as.numeric(sapply(statList, function(a)  as.character(a$V18[4])))
  Unmapped.reads.R1=as.numeric(sapply(statList, function(a)  as.character(a$V2[2])))-Mapped.reads.R1
  Unmapped.reads.R1.pct=1-Mapped.reads.R1.pct
  Unmapped.reads.R2=as.numeric(sapply(statList, function(a)  as.character(a$V2[3])))-Mapped.reads.R2
  Unmapped.reads.R2.pct=1-Mapped.reads.R2.pct
  Mapped.reads=as.numeric(sapply(statList, function(a)  as.character(a$V6[4])))
  Mapped.reads.pct=as.numeric(sapply(statList, function(a)  as.character(a$V7[4])))
  Unmapped.reads=as.numeric(sapply(statList, function(a)  as.character(a$V2[4])))-Mapped.reads
  Unmapped.reads.pct=1-Mapped.reads.pct
  
  sample.names=sub(bam.idx, "", basename(fns.bam))
  
  report=data.frame(sample.names=sample.names, Total.reads=Total.reads, PF.reads=PF.reads, Mapped.reads.R1=Mapped.reads.R1, Mapped.reads.R1.pct=Mapped.reads.R1.pct, Mapped.reads.R2=Mapped.reads.R2, Mapped.reads.R2.pct=Mapped.reads.R2.pct,Paired.Mapped=Paired.Mapped,Paired.Mapped.pct=Paired.Mapped.pct, Unmapped.reads.R1=Unmapped.reads.R1,Unmapped.reads.R1.pct=Unmapped.reads.R1.pct,Unmapped.reads.R2=Unmapped.reads.R2, Unmapped.reads.R2.pct=Unmapped.reads.R2.pct,Mapped.reads=Mapped.reads,Mapped.reads.pct=Mapped.reads.pct,Unmapped.reads=Unmapped.reads, Unmapped.reads.pct=Unmapped.reads.pct )
  report
}

#' @title get.single.end.metrics.report
#' @description Creates a data frame using picard.collectmultiplemetrics()
#' @usage get.single_end.metrics.report(sam.dir, collectmetrics.idx=".alignment_summary_metrics$")
#' @param sam.dir Path to directory with sam files
#' @param collectmetric.idx Index of file (default=".alignment_summary_metrics$")
#' @details Creates a data frame using the txt file, the output of  picard.collectmultiplemetrics(). Used for Single-end data(ex. miRSEQ).
#' @return data frame of the result summary
#' @export
get.single_end.metrics.report=function(sam.dir, collectmetrics.idx=".alignment_summary_metrics$", mc.cores=1){
  
  fns=dir(sam.dir, collectmetrics.idx, recursive = TRUE, full.names = TRUE)
  fns.sam=list.files(sam.dir, ".sam$", recursive = TRUE, full.names = TRUE)
  
  if(length(fns)==0) fns=picard.collectmetrics(fns.bam = fns.sam, ref.fa = ref.fa, mc.cores=mc.cores)
  
  statList=lapply(fns, function(a) read.delim(a, header=FALSE,comment.char = "#"))
  Total.reads=as.numeric(sapply(statList, function(a)  as.character(a$V2[2])))
  PF.reads=as.numeric(sapply(statList, function(a)  as.character(a$V3[2])))
  Mapped.reads=as.numeric(sapply(statList, function(a)  as.character(a$V6[2])))
  Mapped.reads.pct=as.numeric(sapply(statList, function(a)  as.character(a$V7[2])))
  Unmapped.reads=as.numeric(sapply(statList, function(a)  as.character(a$V2[2])))-Mapped.reads
  Unmapped.reads.pct=1-Mapped.reads.pct
  
  sample.names=sub(bam.idx, "", basename(fns.sam))
  
  report=data.frame(sample.names=sample.names, Total.reads=Total.reads, PF.reads=PF.reads,Mapped.reads=Mapped.reads,Mapped.reads.pct=Mapped.reads.pct,Unmapped.reads=Unmapped.reads, Unmapped.reads.pct=Unmapped.reads.pct )
  report
}

#' @title get.annovar.report
#' @description Creates a data frame using the annovar output files
#' @usage get.annovar.report(annot.dir)
#' @param annot.dir Path to directory with annovar output files
#' @details Provide data frame of annotation information
#' @return list of result summary
#' @export 
get.annovar.report=function(annot.dir){
  fns=dir(annot.dir, "annovar$", recursive = TRUE, full.names = TRUE)
  
  sample.names=sub(".annovar", "", basename(fns))
  
  # gatk QC
  library(data.table)
  gatkQC=lapply(fns, function(a) fread(a, sep="\t", select=12))
  gatk.qc=sapply(gatkQC, function(a) table(a$V12))
  
  # variant annot
  fns=dir(annot.dir, "multianno.csv$", recursive = TRUE ,full.names = TRUE)
  annoT=lapply(fns, function(a) fread(a, sep=",", select=c("Ref","Alt","Func.knownGene","Gene.knownGene","ExonicFunc.knownGene")))
  annot.pass=lapply(1:length(gatkQC), function(i) annoT[[i]][which(gatkQC[[i]]$V12=="PASS"),])
  
  sample.variant=lapply(annot.pass, nrow)
  total.variant=format(as.numeric(sum(unlist(sample.variant))),big.mark=",")
  min.variant=format(min(unlist(sample.variant)),big.mark = ",")
  max.variant=format(max(unlist(sample.variant)),big.mark=",")
  mean.variant=format(mean(unlist(sample.variant)),big.mark=",")
  #table(annot.pass[[1]][,"ExonicFunc.knownGene"])[-1]
  func.tbl=sapply(annot.pass, function(a) table(a$Func.knownGene))
  exonicFunc.tbl=sapply(annot.pass, function(a) table(a$ExonicFunc.knownGene)[-1])
  
  if(length(exonicFunc.tbl)!=0) colnames(exonicFunc.tbl)=sample.names
  if(length(exonicFunc.tbl)==0) exonicFunc.tbl=NULL
  if(length(func.tbl)!=0) colnames(func.tbl)=sample.names
  
  mut.type.tbl=sapply(annot.pass, function(a) table(mut.type.somatic(as.data.frame(a))))
  colnames(mut.type.tbl)=sample.names
  
  report=list(Func=func.tbl, ExonicFunc=exonicFunc.tbl, mut.type=mut.type.tbl,total.variant=total.variant,min.variant=min.variant,max.variant=max.variant,mean.variant=mean.variant)
  report
}

#' @title mut.type.somatic
#' @description Write mutation type 
#' @usage mut.type.somatic(df, ref="Ref",alt="Alt")
#' @param df data frame from get.annovar.report()
#' @param ref column name of the reference (default="Ref")
#' @param alt column name of the alteration (default="Alt")
#' @details In the data frame, enter the mutation type using reference and alteration
#' @return data frame with mutation type
#' @export
mut.type.somatic=function(df, ref="Ref",alt="Alt"){
  mut.type=paste(df[,ref],df[,alt], sep=">")
  mut.type[grep("-",mut.type)]="indel"
  mut.type[which(mut.type=="C>A"|mut.type=="G>T")]="C>A/G>T"
  mut.type[which(mut.type=="C>G"|mut.type=="G>C")]="C>G/G>C"
  mut.type[which(mut.type=="C>T"|mut.type=="G>A")]="C>T/G>A"
  mut.type[which(mut.type=="T>A"|mut.type=="A>T")]="T>A/A>T"
  mut.type[which(mut.type=="T>C"|mut.type=="A>G")]="T>C/A>G"
  mut.type[which(mut.type=="T>G"|mut.type=="A>C")]="T>G/A>C"
  mut.type
}

#' @title get.Robject.report
#' @description Reads rda file
#' @usage get.Robject.report(Robject.dir)
#' @param Robject.dir Path to Robject directory
#' @details Reads information according to data set.
#' @return list of data set summary
#' @export
get.Robject.report=function(Robject.dir){
  
  eset.fns=dir(Robject.dir, "eSet.rda$", full.names = TRUE)
  vset.fns=dir(Robject.dir, "filter_vSet.rda$", full.names = TRUE)
  cset.fns=dir(Robject.dir, "cSetList.rda$", full.names = TRUE)
  
  if(grepl("eSet.rda",basename(eset.fns))) {
    eset=get(load(fns[grep("eSet.rda",basename(eset.fns))]))
    total_feature=length(featureNames(eset))
    eSet.eval=TRUE
    expression.set=list(total_feature=total_feature)
  }else if(!grepl("eSet.rda",basename(eset.fns))) {
    eSet.eval=FALSE
    expression.set=NULL}
  
  if(grepl("filter_vSet.rda",basename(vset.fns))) {
    vset=get(load(fns[grep("filter_vSet.rda",basename(vset.fns))]))
    total_variant=length(featureNames(vset))
    vset.eval=TRUE
    variant.set=list(total_variant=total_variant)
  }else if(!grepl("filter_vSet.rda",basename(vset.fns))) {
    vset.eval=FALSE
    variant.set=NULL}
  
  if(grepl("cSetList.rda",basename(cset.fns))) {
    cset=get(load(fns[grep("cSetList.rda",basename(cset.fns))]))
    total_c.variants=length(featureNames(cset))
    cset.eval=TRUE
    c.variant.set=list(total_variant=total_variant)
  }else if(!grepl("cSetList.rda",basename(cset.fns))) {
    cset.eval=FALSE
    c.variant.set=NULL}
  
  report=list(Robject.dir=Robject.dir,expression.set=expression.set,eSet.eval=eSet.eval,variant.set=variant.set,vset.eval=vset.eval,c.variant.set=c.variant.set,cset.eval=cset.eval)
  report
}

#' @title qc_aggregate
#' @description Aggregates the information from the fastq file 
#' @usage qc_aggregate(qc.dir, pattern="fastqc.zip$")
#' @param qc.dir Path to directory with fastqc.zip files
#' @param pattern Index of file (default="fastqc.zip$")
#' @details Combine the fastqc results into one data.
#' @export 
qc_aggregate=function (qc.dir = ".", pattern="fastqc.zip$",progressbar = TRUE) 
{
  qc.files <- list.files(qc.dir, pattern = pattern, full.names = TRUE, recursive = FALSE)
  nfiles <- length(qc.files)
  if (nfiles == 0) 
    stop("Can't find any *fastqc.zip files in the specified qc.dir")
  res.summary <- NULL
  progressbar <- progressbar & nfiles > 3
  if (progressbar) {
    message("Aggregating FastQC Outputs \n")
    pb <- utils::txtProgressBar(max = nfiles, style = 3)
  }
  for (i in 1:nfiles) {
    qc <- qc_read(qc.files[i], modules = c("summary", "statistics", 
                                           "Sequence Duplication Levels"), verbose = FALSE)
    .summary <- qc$summary
    .statistics <- as.data.frame(qc$basic_statistics)
    rownames(.statistics) <- .statistics$Measure
    pct.dup <- round(100 - qc$total_deduplicated_percentage, 
                     2)
    .summary <- dplyr::mutate(.summary, tot.seq = rep(.statistics["Total Sequences", 
                                                                  2], nrow(.summary)), pct.gc = as.numeric(rep(.statistics["%GC", 
                                                                                                                           2], nrow(.summary))), seq.length = rep(.statistics["Sequence length", 
                                                                                                                                                                              2], nrow(.summary)), pct.dup = rep(pct.dup, nrow(.summary)))
    res.summary <- rbind(res.summary, .summary)
    if (progressbar) 
      utils::setTxtProgressBar(pb, i)
  }
  res.summary <- dplyr::select_(res.summary, "sample", "module", 
                                "status", "tot.seq", "seq.length", "pct.gc", "pct.dup")
  res.summary$sample <- gsub(".fastq.gz|.fastq", "", res.summary$sample, 
                             ignore.case = TRUE)
  if (progressbar) 
    close(pb)
  res.summary <- structure(res.summary, class = c("qc_aggregate", 
                                                  class(res.summary)))
  res.summary
}

