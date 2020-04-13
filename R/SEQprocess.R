

#' @title SEQprocess
#' @description Run the NGS data processing pipeline
#' @param fastq.dir If the user starts the process with fastq files, set the directory for the fastq files.
#' @param output.dir Output directory
#' @param argList The argument list used by the user in the shell
#' @param project.name User's project name
#' @param type Sequence data type
#' @param pipeline One of the six pipelines provided by SEQprocess
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @param report.mode Whether the process is finished and report is generated
#' @param config.fn Congifure file path
#' @param qc Whether quality check 
#' @param trim.method Set trimming method
#' @param align.method Set alignment method
#' @param build.transcriptome.idx (tophat) A transcriptome index and the associated data files (the original GFF file) can be thus reused for multiple TopHat runs with this option, so these files are only created for the first run with a given set of transcripts. (default=FALSE)
#' @param tophat.thread.number (tophat) A numeric value of the number of threads
#' @param bwa.method (bwa) Set bwa method
#' @param bwa.thread.number (bwa) A numeric value of the number of threads 
#' @param star.thread.number (STAR) A numeric value of the number of threads
#' @param rm.dup Set the remove duplicates method
#' @param realign Whether realignment
#' @param variant.call.method Set variant call method
#' @param annotation.method Set variant annotation method
#' @param ref (annovar) Set annovar reference version
#' @param rseq.abundance.method Set RNA quantification method
#' @param cufflinks.gtf (cufflinks) If you set "-G", Output will not include novel genes and isoforms that are assembled.
#' @param cufflinks.thread.number (cufflinks) A numeric value of the number of threads 
#' @param RNAtype (htseq) Choose mRNA or miRNA.
#' @param CNV Whether estimate copy number variation
#' @param make.eSet Make ExpressionSet R data(RNA expression, Copy number variation, Mutation)
#' @param eset2SE Convert ExpressionSet R data to SummarizedExperiment R data
#' @export
SEQprocess=function(
  fastq.dir=NULL, # fastq directory, if starting with fastQC set fastq directory
  output.dir=file.path(getwd(),"result", "SEQprocess_result"), # output directory
  argList=list(program.name="SEQprocess"),
  project.name="SEQprocess",
  
  type=c("WGS","WES","BarSEQ", "RSEQ","miRSEQ"), # flatform type
  pipeline=c("none", "GDC", "GATK", "BarSEQ", "Tuxedo", "miRSEQ"), # pipeline 
  mc.cores=1,
  run.cmd=TRUE,
  report.mode=FALSE,
  config.fn=system.file("data/config.R", package = "SEQprocess"),
  
  # QC method
  qc=TRUE, # fastqc
  trim.method=c("trim.galore", "cutadapt", "none"), # trimming
  
  ## alignment
  align.method=c("bwa", "bowtie2", "tophat2", "star", "none"), #align
  
  #tophat opt
  build.transcriptome.idx=FALSE, # make transcriptome index
  tophat.thread.number=4, # thread number
  
  #bwa opt
  bwa.method=c("mem", "aln"),
  bwa.thread.number=4,
  
  #star opt
  star.thread.number=8,
  
  ## rmdu
  rm.dup=c("MarkDuplicates","BARCODE","none"),
  
  ## realign
  realign=TRUE,
  
  ## variant call
  #RSEQ.variant.call=FALSE,
  variant.call.method=c("none", "gatk", "varscan2", "mutect2", "muse", "somaticsniper"),
  gatk.thread.number=4,
  
  ## annotation
  annotation.method=c("annovar", "vep", "none"),
  
  #annovar opt
  ref="hg38",
  
  ## expression
  rseq.abundance.method=c("none", "cufflinks", "htseq"),
  
  #cufflinks opt
  cufflinks.gtf=c("G", "g"),
  cufflinks.thread.number=4,
  
  #htseq opt
  RNAtype=c("mRNA", "miRNA"),
  
  ##copy number variation call
  CNV=FALSE,
  
  #make Set
  make.eSet=FALSE,
  eset2SummarizedExperiment=FALSE,
  mut.cnt.cutoff=8,
  
  #dirs 
  qc.dir=file.path(output.dir, "00_qc"),
  trim.dir=file.path(output.dir, "01_trim"),
  align.dir=file.path(output.dir, "02_align"),
  rmdup.dir=file.path(output.dir, "03_rmdup"),
  realign.dir=file.path(output.dir, "04_realign"),
  vcf.dir=file.path(output.dir, "05_vcf"),
  annot.dir=file.path(output.dir, "06_annot"),
  RNAquant.dir=file.path(output.dir, "07_RNAquant"),
  cnv.dir=file.path(output.dir, "08_cnv"),
  Robject.dir=file.path(output.dir, "09_Robject")
){
  # set default option
  type=match.arg(type)
  pipeline=match.arg(pipeline)
  trim.method=match.arg(trim.method)
  align.method=match.arg(align.method)
  bwa.method=match.arg(bwa.method)
  rseq.abundance.method=match.arg(rseq.abundance.method)
  RNAtype=match.arg(RNAtype)
  ref=match.arg(ref)
  cufflinks.gtf=match.arg(cufflinks.gtf)
  rm.dup=match.arg(rm.dup)
  variant.call.method=match.arg(variant.call.method)
  annotation.method=match.arg(annotation.method)
  
  for (i in 1:length(argList)) assign(names(argList)[i], argList[[i]])
  source(config.fn, local=.GlobalEnv)
  
  ## MuSE opt
  if(type=="WGS") MuSE.data.type="G" else if(type=="WES") MuSE.data.type="E"
  
  # get param
  # dseq-GDC
  if(type=="WGS"&pipeline=="GDC"|type=="WES"&pipeline=="GDC"){
    if(variant.call.method=="none") variant.call.method=="varscan2"
    annotation.method="vep" 
    # Other options are default
  }  
  
  # dseq-GATK
  if(type=="WGS"&pipeline=="GATK"|type=="WES"&pipeline=="GATK"){
    variant.call.method="gatk" 
    # Other options are default
  } 
  
  # dseq-barSeq
  if(type=="BarSEQ"&pipeline=="BarSEQ"){
    rm.dup="BARCODE"
    variant.call.method="gatk" 
    # Other options are default
  }  
  
  # rseq-GDC
  if(type=="RSEQ"&pipeline=="GDC"){
    trim.method="none"
    align.method="star"
    rseq.abundance.method="htseq"
    realign=FALSE
    rm.dup="none"
    annotation.method="none"
    # Other options are default
  }  
  
  # rseq-tuxedo
  if(pipeline=="Tuxedo"){
    type="RSEQ"
    trim.method="trim.galore"
    align.method="tophat2"
    rseq.abundance.method="cufflinks"
    realign=FALSE
    rm.dup="none" 
    annotation.method="none"
    # Other options are default
  }  
  
  # rseq-miRseq
  if(pipeline=="miRSEQ"){
    type="miRSEQ"
    trim.method="cutadapt"
    bwa.method="aln"
    rseq.abundance.method="htseq"
    RNAtype="miRNA"
    realign=FALSE
    rm.dup="none"
    annotation.method="none"
    # Other options are default
  }
  
  #preset 
  re.order=ifelse(variant.call.method!="none",TRUE, FALSE)
  #picard rmdu opt
  BARCODE_TAG=ifelse(rm.dup=="MarkDuplicates", FALSE, TRUE)
  #gatk opt
  unsafe=ifelse(type=="WGS"|type=="WES"|type=="BarSEQ", FALSE,TRUE)
  
  # Set process name, output directory, input files and sample name 
  proc.names=get.process.names(qc, trim.method, align.method, bwa.method , rm.dup, realign, variant.call.method, annotation.method, rseq.abundance.method, make.eSet)
  
  if(qc|trim.method!="none") input.fns=list(fq1=get.fns(fastq.dir,fq1.idx), fq2=get.fns(fastq.dir,fq2.idx))
  if(qc|trim.method!="none") sample.name=sub(fq1.idx, "", basename(input.fns[[1]]))
  
  # make directory 
  dir.create(qc.dir, recursive = TRUE, showWarnings=FALSE)
  dir.create(trim.dir,recursive = TRUE, showWarnings=FALSE)
  dir.create(align.dir,recursive = TRUE, showWarnings=FALSE)
  dir.create(rmdup.dir,recursive = TRUE, showWarnings=FALSE)
  dir.create(realign.dir,recursive = TRUE, showWarnings=FALSE)
  dir.create(vcf.dir,recursive = TRUE, showWarnings=FALSE)
  dir.create(annot.dir,recursive = TRUE, showWarnings=FALSE)
  if(rseq.abundance.method!="none") dir.create(RNAquant.dir,recursive = TRUE, showWarnings=FALSE)
  dir.create(cnv.dir, recursive=TRUE, showWarnings=FALSE)
  dir.create(Robject.dir,recursive=TRUE, showWarnings=FALSE)
  envList=as.list(environment())
  
  ########################################################
  ###################### System Run ######################
  ########################################################
  
  
  
  # Quality Check - FASTQC
  if(qc) {
    fns.fq1=input.fns[[1]]
    fns.fq2=input.fns[[2]]
    fastQC(fq1=fns.fq1, fq2=fns.fq2, output.dir=qc.dir, run.cmd, mc.cores)
  }
  
  
  
  # Trimming - Trim Galore or Cutadapt
  if(trim.method!="none") {
    trim.fns=run.trim(trim.method=trim.method, fq1=input.fns[[1]], fq2=input.fns[[2]], output.dir=trim.dir, trim.quality=trim.quality, 
                      trim.clip_R1=trim.clip_R1, trim.clip_R2=trim.clip_R2, m=m, adapt.seq=adapt.seq, sample.name=sample.name, 
                      run.cmd=run.cmd, mc.cores=mc.cores)
    fns.fq1=trim.fns[[1]]
    fns.fq2=trim.fns[[2]]
  }
  
  
  # Alignment - 1. BWA
  if(align.method=="bwa") {
    fns.bam=run.bwa(bwa.method, fq1=fns.fq1, fq2=fns.fq2, output.dir=align.dir, sample.name, ref.fa, bwa.idx, bwa.thread.number, RGSM, run.cmd, mc.cores, trim.dir)
    fns.bam
  }
  
  
  # Alignment - 2. STAR
  if(align.method=="star"){
    if(trim.method!="none") fastq.dir=trim.dir
    fns.bam=run.star(star.idx.dir, sample.name, fq1, fq2, fastq.dir, output.dir=align.dir, ref.fa, ref.gtf, sjdbOverhang, star.thread.number, outFilterMultimapScoreRange, outFilterMultimapNmax, outFilterMismatchNmax, alignIntronMax,  alignMatesGapMax, sjdbScore, alignSJDBoverhangMin, outFilterMatchNminOverLread, outFilterScoreMinOverLread, mc.cores, run.cmd, align.dir, RGSM)
    fns.bam
  }
  
  # Alignment - 3. Tophat2
  if(align.method=="tophat2"){
    fns.bam=run.tophat2(fq1, fq2, output.dir=align.dir, sample.name, ref.gtf, bowotie.idx, tophat.thread.number, build.transcriptome.idx, run.cmd, mc.cores, trim.dir, RGSM)
    fns.bam
  }
  
  # Alignment - 4. bowtie2
  if(align.method=="bowtie2"){
    fns.bam=run.bowtie2(fq1, output.dir=align.dir, sample.name, bowtie.idx, mc.cores, run.cmd)
    fns.bam
  }
  
  # Add Read groups & Reorder - Picard
  if(re.order){
    if(type=="WGS"|type=="WES"|type=="BarSEQ") fns.bam=get.fns(align.dir, idx=".bam$") else if(type=="RSEQ"|type=="miRSEQ") fns.bam=get.fns(align.dir, idx="accepted_hits.bam$")
    sample.name=basename(dirname(fns.bam))
    fns.bam=picard.addrg(fns.bam=fns.bam, output.dir=align.dir, RGSM=sample.name, run.cmd=run.cmd, mc.cores=mc.cores)
    fns.bam=picard.reorder(fns.bam=fns.bam, output.dir=align.dir, ref.fa=ref.fa, run.cmd=run.cmd, mc.cores=mc.cores)
    fns.bam
  }
  
  
  # Remove Duplicates - Picard
  if(rm.dup!="none") {
    fns.bam=get.fns(align.dir, idx=".rg.od.bam$")
    if(type=="WGS"|type=="WES"|type=="BarSEQ") sample.name=sub(bam.idx, "", basename(fns.bam)) else if(type=="RSEQ"|type=="miRSEQ") sample.name=basename(dirname(fns.bam))

    fns.bam=picard.rmdu(fns.bam=fns.bam, output.dir=rmdup.dir, type=type, sample.name=sample.name, BARCODE_TAG=BARCODE_TAG , run.cmd=run.cmd, mc.cores=mc.cores) 
    fns.bam
  }
  
  
  # Re-alingment 
  if(realign){
    fns.bam=run.realign(fns.bam, output.dir=realign.dir, sample.name, ref.fa=ref.fa, ref.gold_indels=ref.gold_indels, ref.dbSNP, unsafe=unsafe, gatk.thread.number=gatk.thread.number, run.cmd=run.cmd, mc.cores=mc.cores, rmdup.dir=rmdup.dir)
    fns.bam
  }
  
  
  # VCF
  if(variant.call.method=="gatk"){
    fns.vcf=run.gatk(fns.bam=rmdu.fns, output.dir=vcf.dir, sample.name=sample.name, ref.fa=ref.fa, ref.gold_indels=ref.gold_indels, ref.dbSNP, unsafe=unsafe, FS=FS, QD=QD, QUAL=QUAL, DP=DP, gatk.window=gatk.window, cluster=cluster, gatk.thread.number=gatk.thread.number, run.cmd=run.cmd, mc.cores=mc.cores, realign.dir=realign.dir)
    fns.vcf
  }
  
  if(variant.call.method=="varscan2"){
    fns.vcf=run.varscan(fns.bam=fns.bam, ref.dbSNP=ref.dbSNP, ref.gold_indels=ref.gold_indels, unsafe=unsafe, ref.fa=ref.fa, normal.bam=normal.bam, tumor.bam=tumor.bam, sample.name=sample.name, output.dir=vcf.dir, mapQ=mapQ, min_coverage_normal=min_coverage_normal, min_coverage_tumor=min_coverage_tumor, min_var_freq=min_var_freq, min_freq_for_hom=min_freq_for_hom, p_value=p_value, somatic_p_value=somatic_p_value, strand_filter=strand_filter, min_tumor_freq=min_tumor_freq, max_normal_freq=max_normal_freq, gatk.thread.number=gatk.thread.number, run.cmd=run.cmd, mc.cores=mc.cores, realign.dir=realign.dir)
    fns.vcf
  }
  
  if(variant.call.method=="muse"){
    fns.vcf=run.muse(fns.bam=fns.bam, output.dir=vcf.dir, tumor.bam=tumor.bam, normal.bam=normal.bam, sample.name=sample.name, ref.fa=ref.fa, ref.dbSNP=ref.dbSNP, ref.gold_indels=ref.gold_indels, MuSE.data.type=MuSE.data.type, unsafe=unsafe, gatk.thread.number=gatk.thread.number, run.cmd=run.cmd, mc.cores=mc.cores, realign.dir=realign.dir)
    fns.vcf
  }
  
  
  if(variant.call.method=="mutect2"){
    fns.vcf=run.mutect2(tumor.bam=tumor.bam, normal.bam=normal.bam, output.dir=vcf.dir, sample.name=sample.name, ref.fa=ref.fa, ref.dbSNP=ref.dbSNP, cosmic.vcf=cosmic.vcf, minN=minN, filteredrecordsmergetype=filteredrecordsmergetype, unsafe=unsafe, gatk.thread.number=gatk.thread.number, run.cmd=run.cmd, mc.cores=mc.cores, realign.dir=realign.dir)
    fns.vcf
  }
  
  if(variant.call.method=="somaticsniper"){
    fns.vcf=run.somaticsniper(fns.bam=fns.bam, output.dir=vcf.dir, sample.name=sample.name, ref.fa=ref.fa, ref.dbSNP=ref.dbSNP, ref.gold_indels=ref.gold_indels, unsafe=unsafe, mapQual=mapQual, LOH=LOH, Genotype=Genotype, somaticQual=somaticQual, somaticMutation=somaticMutation, Theta=Theta, Hap.number=Hap.number, Hap.diff=Hap.diff, out.format=out.format, gatk.thread.number=gatk.thread.number, run.cmd=run.cmd, mc.cores=mc.cores, realign.dir=realign.dir)
    fns.vcf
  }
  
  # annot
  if(annotation.method!="none"){
    fns.annot=run.annotation(annotation.method=annotation.method, fns.vcf=fns.vcf, output.dir=annot.dir, sample.name=sample.name, annovar.db.dir=annovar.db.dir, vep.db.dir=vep.db.dir, ref=ref, protocol=protocol, protocol.type=protocol.type, nastring=nastring, run.cmd=run.cmd, mc.cores=mc.cores, vcf.dir=vcf.dir)
    fns.annot
  }
  
  
  # RNAquant
  if(rseq.abundance.method!="none"){
    out.fns=run.rseq.abundance(rseq.abundance.method=rseq.abundance.method, fns.bam=fns.bam, output.dir=RNAquant.dir, sample.name=sample.name, cufflinks.thread.number=cufflinks.thread.number, cufflinks.gtf=cufflinks.gtf, ref.gtf=ref.gtf, RNAtype=RNAtype, Mode=Mode, stranded=stranded, idattr=idattr, htseq.r=htseq.r, htseq.a=htseq.a, mir.gff=mir.gff, run.cmd=run.cmd, mc.cores=mc.cores, align.dir=align.dir)
  }
  
  # cnv
  if(CNV&grepl("WGS|WES|BarSEQ", type)){
    run.cnv(ref.fa=ref.fa, fns.bam=fns.bam, rmdup.dir=rmdup.dir, sample.name=sample.name, mapQ=mapQ, output.dir=cnv.dir, window=window, cnv.dir=cnv.dir, run.cmd, mc.cores)
  } else if(CNV&grepl("RSEQ|miRSEQ", type)){
    message("CNV estimation will not run. It requires DSEQ data")
  }
  
  # Make ExpressionSet, VariantSet
  if(make.eSet){
    if(grepl("RSEQ|miRSEQ", type)){
      eset=make.eset(RNAquant.dir=RNAquant.dir, Robject.dir=Robject.dir)
      if(eset2SummarizedExperiment) eset2SE(eset=eset, vset=NULL, cset=NULL, Robject.dir=Robject.dir)
      } 
    if(length(get.fns(annot.dir, idx=".csv$")) > 0){
      vset=run.make.vset(annot.dir=annot.dir, Robject.dir=Robject.dir, vcf.dir=vcf.dir, ref.fa=ref.fa, unsafe=unsafe, minBaseQuality=minBaseQuality, minMappingQuality=minMappingQuality, mut.cnt.cutoff=mut.cnt.cutoff, run.cmd=run.cmd, mc.cores=mc.cores)
      if(eset2SummarizedExperiment) eset2SE(eset=NULL, vset=vset, cset=NULL, Robject.dir=Robject.dir)
      }
    if(CNV){
      cset=make.cset(cnv.dir=cnv.dir, refGene.gr=refGene.gr)
      if(eset2SummarizedExperiment) eset2SE(eset=NULL, vset=NULL, cset=cset, Robject.dir=Robject.dir)
      }
} 
  
  
  # Report
  if(report.mode) report(envList)
}


#' @export
run.trim=function(trim.method, fq1, fq2, output.dir, trim.quality, trim.clip_R1, trim.clip_R2, m, adapt.seq, sample.name, run.cmd, mc.cores){
  if(trim.method=="trim.galore"){
    out.fns=trim.gal(fq1=fq1, fq2=fq2, output.dir=output.dir, trim.quality=trim.quality, trim.clip_R1=trim.clip_R1, trim.clip_R2=trim.clip_R2, run.cmd=run.cmd, mc.cores=mc.cores)
    out.fns
  }
  if(trim.method=="cutadapt"){
    out.fns=cutadapt(fq1=fq1, output.dir=output.dir, sample.name=sample.name, adapt.seq=adapt.seq, m=m, run.cmd=run.cmd, mc.cores=mc.cores)
    out.fns
  }
}



#' @export
run.bwa=function(bwa.method, fq1, fq2, output.dir, sample.name, ref.fa, bwa.idx, bwa.thread.number, RGSM, run.cmd, mc.cores, trim.dir){
  fq1.fns=get.fns(trim.dir, fq1.idx)
  fq2.fns=get.fns(trim.dir, fq2.idx)
  sample.name=sub(fq1.idx, "", basename(fq1.fns))
  fns.bam=bwa(bwa.method=bwa.method, fq1=fq1.fns, fq2=fq2.fns, output.dir=output.dir, sample.name=sample.name, ref.fa=ref.fa, bwa.idx=bwa.idx, bwa.thread.number=bwa.thread.number, run.cmd=run.cmd, mc.cores=mc.cores)
}  


#' @export
run.star=function(star.idx.dir, sample.name, fq1, fq2, fastq.dir, output.dir, ref.fa, ref.gtf, sjdbOverhang, star.thread.number, outFilterMultimapScoreRange, outFilterMultimapNmax, outFilterMismatchNmax, alignIntronMax,  alignMatesGapMax, sjdbScore, alignSJDBoverhangMin, outFilterMatchNminOverLread, outFilterScoreMinOverLread, mc.cores, run.cmd, align.dir, RGSM){
  
  fq1.fns=get.fns(fastq.dir, fq1.idx)
  fq2.fns=get.fns(fastq.dir, fq2.idx)
  sample.name=sub(fq1.idx, "", basename(fq1.fns))
  
  #build.star.idx(star.idx.dir=star.idx.dir, sample.name=sample.name, align.dir=align.dir, ref.fa=ref.fa, ref.gtf=ref.gtf, sjdbOverhang=sjdbOverhang, star.thread.number=star.thread.number, fasta.idx=TRUE, SJ.idx=FALSE, run.cmd=run.cmd)
  
  #STAR(star.idx.dir=star.idx.dir, output.dir=output.dir, sample.name=sample.name, fq1=fq1.fns, fq2=fq2.fns, star.thread.number=star.thread.number, outFilterMultimapScoreRange=outFilterMultimapScoreRange, outFilterMultimapNmax=outFilterMultimapNmax, outFilterMismatchNmax=outFilterMismatchNmax, alignIntronMax=alignIntronMax,  alignMatesGapMax=alignMatesGapMax, sjdbScore=sjdbScore, alignSJDBoverhangMin=alignSJDBoverhangMin, outFilterMatchNminOverLread=outFilterMatchNminOverLread, outFilterScoreMinOverLread=outFilterScoreMinOverLread, sjdbOverhang=sjdbOverhang, SJ.detect=TRUE, SJ.align=FALSE, run.cmd=run.cmd, mc.cores=mc.cores)
  
  #build.star.idx(star.idx.dir=star.idx.dir, sample.name=sample.name, align.dir=align.dir, ref.fa=ref.fa, ref.gtf=ref.gtf, sjdbOverhang=sjdbOverhang, star.thread.number=star.thread.number, fasta.idx=FALSE, SJ.idx=TRUE, run.cmd=run.cmd)
  
  fns.bam=STAR(star.idx.dir=star.idx.dir, output.dir=output.dir, sample.name=sample.name, fq1=fq1.fns, fq2=fq2.fns, star.thread.number=star.thread.number, outFilterMultimapScoreRange=outFilterMultimapScoreRange, outFilterMultimapNmax=outFilterMultimapNmax, outFilterMismatchNmax=outFilterMismatchNmax, alignIntronMax=alignIntronMax, alignMatesGapMax=alignMatesGapMax, sjdbScore=sjdbScore, alignSJDBoverhangMin=alignSJDBoverhangMin, outFilterMatchNminOverLread=outFilterMatchNminOverLread, outFilterScoreMinOverLread=outFilterScoreMinOverLread, sjdbOverhang=sjdbOverhang, SJ.detect=FALSE, SJ.align=TRUE, run.cmd=run.cmd, mc.cores=mc.cores)
  
  fns.bam
}


#' @export
run.tophat2=function(fq1, fq2, output.dir, sample.name, ref.gtf, bowotie.idx, tophat.thread.number, build.transcriptome.idx, run.cmd, mc.cores, trim.dir, RGSM){
  fq1.fns=get.fns(trim.dir, fq1.idx)
  fq2.fns=get.fns(trim.dir, fq2.idx)
  sample.name=sub(fq1.idx, "", basename(fq1.fns))
  
  fns.bam=tophat2(fq1=fq1.fns, fq2=fq2.fns, output.dir=output.dir, sample.name=sample.name, ref.gtf=ref.gtf, bowtie.idx=bowtie.idx, tophat.thread.number=tophat.thread.number, build.transcriptome.idx=build.transcriptome.idx, run.cmd=run.cmd, mc.cores=mc.cores)
  
  fns.bam
}


#' @export
run.bowtie2=function(fq1, output.dir, sample.name, bowtie.idx, mc.cores, run.cmd){
  fq1.fns=get.fns(trim.dir, fq1.idx)
  sample.name=sub(fq1.idx, "", basename(fq1.fns))
  
  fns.bam=bowtie2(fq1=fq1.fns, output.dir=output.dir, sample.name=sample.name, bowtie.idx=bowtie.idx, mc.cores=mc.cores, run.cmd=run.cmd)
  fns.bam
}  

#' @export
run.realign=function(fns.bam, output.dir, sample.name, ref.fa=ref.fa, ref.gold_indels=ref.gold_indels, ref.dbSNP, unsafe=unsafe, gatk.thread.number=gatk.thread.number, run.cmd=run.cmd, mc.cores=mc.cores, rmdup.dir=rmdup.dir){
  
  fns.bam=get.fns(rmdup.dir, bam.idx)
  sample.name=sub(bam.idx, "", basename(fns.bam))
  fns.interval=gatk.targetcreator(fns.bam=fns.bam, output.dir=output.dir, sample.name=sample.name, ref.fa=ref.fa, ref.gold_indels=ref.gold_indels, unsafe=unsafe, gatk.thread.number=gatk.thread.number, run.cmd=run.cmd, mc.cores=mc.cores)
  out.fns=gatk.realign(fns.bam=fns.bam, fns.interval=fns.interval, output.dir=output.dir, sample.name=sample.name, ref.fa=ref.fa, ref.gold_indels=ref.gold_indels, unsafe=unsafe, run.cmd=run.cmd, mc.cores=mc.cores)
  
  out.fns
}


#' @export
run.gatk=function(fns.bam, output.dir, sample.name, ref.fa=ref.fa, ref.gold_indels=ref.gold_indels, ref.dbSNP, unsafe=unsafe, FS, QD, QUAL, DP, gatk.window, cluster, gatk.thread.number, run.cmd=run.cmd, mc.cores=mc.cores, realign.dir=realign.dir){
  
    fns.bam=get.fns(realign.dir, bam.idx)
    sample.name=sub(bam.idx, "", basename(fns.bam))
    fns.grp=gatk.baserecalibrator(fns.bam=fns.bam, output.dir=output.dir, sample.name=sample.name, ref.fa=ref.fa, ref.dbSNP=ref.dbSNP, ref.gold_indels=ref.gold_indels, unsafe=unsafe, gatk.thread.number=gatk.thread.number, run.cmd=run.cmd, mc.cores=mc.cores)
    fns.bam=gatk.applyBQSR(fns.bam=fns.bam, fns.grp=fns.grp, output.dir=output.dir, sample.name=sample.name, ref.fa=ref.fa, unsafe=unsafe, gatk.thread.number=gatk.thread.number, run.cmd=run.cmd, mc.cores=mc.cores)
    
    fns.vcf=gatk.haplotypecaller(fns.bam=fns.bam, output.dir=output.dir, sample.name=sample.name, ref.fa=ref.fa, unsafe=unsafe, run.cmd=run.cmd, mc.cores=mc.cores)
    out.fns=gatk.variantfilter(fns.vcf=fns.vcf, output.dir=output.dir, ref.fa=ref.fa, FS=FS, QD=QD, QUAL=QUAL, DP=DP, gatk.window=gatk.window, cluster= cluster, run.cmd=run.cmd, mc.cores=mc.cores)
    out.fns
}  

#' @export
run.varscan=function(fns.bam, ref.dbSNP, ref.gold_indels, unsafe, ref.fa, normal.bam, tumor.bam, sample.name, output.dir, mapQ, min_coverage_normal, min_coverage_tumor, min_var_freq, min_freq_for_hom, p_value, somatic_p_value, strand_filter, min_tumor_freq, max_normal_freq, gatk.thread.number, run.cmd, mc.cores, realign.dir){
    
    fns.bam=get.fns(realign.dir, bam.idx)  
    sample.name=sub(bam.idx, "", basename(fns.bam))
    fns.grp=gatk.baserecalibrator(fns.bam=fns.bam, output.dir=output.dir, sample.name=sample.name, ref.fa=ref.fa, ref.dbSNP=ref.dbSNP, ref.gold_indels=ref.gold_indels, unsafe=unsafe, gatk.thread.number=gatk.thread.number, run.cmd=run.cmd, mc.cores=mc.cores)
    fns.bam=gatk.applyBQSR(fns.bam=fns.bam, fns.grp=fns.grp, output.dir=output.dir, sample.name=sample.name, ref.fa=ref.fa, unsafe=unsafe, gatk.thread.number=gatk.thread.number, run.cmd=run.cmd, mc.cores=mc.cores)
    
    normal.bam=fns.bam[grep(n.sample, fns.bam)]
    tumor.bam=fns.bam[grep(t.sample, fns.bam)]
    fn.pileup=multiple.reads.pileup(ref.fa=ref.fa, normal.bam=normal.bam, tumor.bam=tumor.bam, sample.name=sample.name, output.dir=output.dir, run.cmd=run.cmd, mc.cores=mc.cores, mapQ=mapQ)
    fns.vcf=varscan(fn.pileup=fn.pileup, output.dir=output.dir, sample.name=sample.name, min_coverage_normal=min_coverage_normal, min_coverage_tumor=min_coverage_tumor, min_var_freq=min_var_freq, min_freq_for_hom=min_freq_for_hom, somatic_p_value=somatic_p_value, strand_filter=strand_filter, run.cmd=run.cmd, mc.cores=mc.cores)
    fns.vcf=get.fns(output.dir, idx="indel.vcf$|snp.vcf$")
    processSomatic(fns.vcf=fns.vcf, min_tumor_freq=min_tumor_freq, output.dir=output.dir, max_normal_freq=max_normal_freq, p_value=p_value, run.cmd=run.cmd, mc.cores=mc.cores)
                              
}


#' @export
run.muse=function(fns.bam, output.dir, tumor.bam, normal.bam, sample.name, ref.fa, ref.dbSNP, ref.gold_indels, MuSE.data.type, unsafe, gatk.thread.number, run.cmd, mc.cores, realign.dir){
  
  fns.bam=get.fns(realign.dir, bam.idx)  
  sample.name=sub(bam.idx, "", basename(fns.bam))
  fns.grp=gatk.baserecalibrator(fns.bam=fns.bam, output.dir=output.dir, sample.name=sample.name, ref.fa=ref.fa, ref.dbSNP=ref.dbSNP, ref.gold_indels=ref.gold_indels, unsafe=unsafe, gatk.thread.number=gatk.thread.number, run.cmd=run.cmd, mc.cores=mc.cores)
  fns.bam=gatk.applyBQSR(fns.bam=fns.bam, fns.grp=fns.grp, output.dir=output.dir, sample.name=sample.name, ref.fa=ref.fa, unsafe=unsafe, gatk.thread.number=gatk.thread.number, run.cmd=run.cmd, mc.cores=mc.cores)
  
  normal.bam=fns.bam[grep(n.sample, fns.bam)]
  tumor.bam=fns.bam[grep(t.sample, fns.bam)]
  out.txt=MuSE.call(tumor.bam=tumor.bam, normal.bam=normal.bam, output.dir=output.dir, sample.name=sample.name, ref.fa=ref.fa, run.cmd=run.cmd, mc.cores=mc.cores)
  MuSE.sump(MuSE.txt=out.txt, output.dir=output.dir, MuSE.data.type=MuSE.data.type, ref.dbSNP=ref.dbSNP, ref.gold_indels=ref.gold_indels, run.cmd=run.cmd, mc.cores=mc.cores)
}

#' @export
run.mutect2=function(tumor.bam, normal.bam, output.dir, sample.name, ref.fa, ref.dbSNP, cosmic.vcf, minN, filteredrecordsmergetype, unsafe, gatk.thread.number, run.cmd, mc.cores, realign.dir){
  
  fns.bam=get.fns(realign.dir, bam.idx)  
  sample.name=sub(bam.idx, "", basename(fns.bam))
  fns.grp=gatk.baserecalibrator(fns.bam=fns.bam, output.dir=output.dir, sample.name=sample.name, ref.fa=ref.fa, ref.dbSNP=ref.dbSNP, ref.gold_indels=ref.gold_indels, unsafe=unsafe, gatk.thread.number=gatk.thread.number, run.cmd=run.cmd, mc.cores=mc.cores)
  fns.bam=gatk.applyBQSR(fns.bam=fns.bam, fns.grp=fns.grp, output.dir=output.dir, sample.name=sample.name, ref.fa=ref.fa, unsafe=unsafe, gatk.thread.number=gatk.thread.number, run.cmd=run.cmd, mc.cores=mc.cores)
  
  normal.bam=fns.bam[grep(n.sample, fns.bam)]
  tumor.bam=fns.bam[grep(t.sample, fns.bam)]
  
  out.vcf=gatk.mutect2.normal(normal.bam=normal.bam, sample.name=sample.name, ref.dbSNP=ref.dbSNP, cosmic.vcf=cosmic.vcf, output.dir=output.dir, run.cmd=run.cmd, mc.cores=mc.cores)
  pon.vcf=gatk.combinevariants(ref.fa=ref.fa, normal.vcf=out.vcf, minN=minN, filteredrecordsmergetype=filteredrecordsmergetype, output.dir=output.dir)
  out.vcf=mutect2(output.dir=output.dir, ref.fa=ref.fa, tumor.bam=tumor.bam, normal.bam=normal.bam, pon.vcf=pon.vcf, cosmic.vcf=cosmic.vcf, ref.dbSNP=ref.dbSNP, contamination_fraction_to_filter=contamination_fraction_to_filter)
  out.vcf
}


#' @export
run.somaticsniper=function(fns.bam, output.dir, sample.name, ref.fa, ref.dbSNP, ref.gold_indels, unsafe, mapQual, LOH, Genotype, somaticQual, somaticMutation, Theta, Hap.number, Hap.diff, out.format, gatk.thread.number, run.cmd, mc.cores, realign.dir){
  
  fns.bam=get.fns(realign.dir, bam.idx)  
  sample.name=sub(bam.idx, "", basename(fns.bam))
  fns.grp=gatk.baserecalibrator(fns.bam=fns.bam, output.dir=output.dir, sample.name=sample.name, ref.fa=ref.fa, ref.dbSNP=ref.dbSNP, ref.gold_indels=ref.gold_indels, unsafe=unsafe, gatk.thread.number=gatk.thread.number, run.cmd=run.cmd, mc.cores=mc.cores)
  fns.bam=gatk.applyBQSR(fns.bam=fns.bam, fns.grp=fns.grp, output.dir=output.dir, sample.name=sample.name, ref.fa=ref.fa, unsafe=unsafe, gatk.thread.number=gatk.thread.number, run.cmd=run.cmd, mc.cores=mc.cores)
  
  normal.bam=fns.bam[grep(n.sample, fns.bam)]
  tumor.bam=fns.bam[grep(t.sample, fns.bam)]
  
  out.vcf=somaticsniper(tumor.bam=tumor.bam, normal.bam=normal.bam, output.dir=output.dir, sample.name=sample.name, ref.fa=ref.fa, mapQual=mapQual, LOH=LOH, Genotype=Genotype, somaticQual=somaticQual, somaticMutation=somaticMutation, Theta=Theta, Hap.number=Hap.number, Hap.diff=Hap.diff, out.format=out.format, run.cmd=run.cmd, mc.cores=mc.cores)
  out.vcf
}



#' @export
run.annotation=function(annotation.method, fns.vcf, output.dir, sample.name, annovar.db.dir, vep.db.dir, ref, protocol, protocol.type, nastring, run.cmd, mc.cores, vcf.dir=vcf.dir){
  fns.vcf=get.fns(vcf.dir, vcf.idx)
  sample.name=sub(vcf.idx, "", basename(fns.vcf))
  if(annotation.method=="annovar"){
    fns.annovar=vcf2annovar(fns.vcf=fns.vcf, output.dir=output.dir, sample.name=sample.name, run.cmd=run.cmd, mc.cores=mc.cores)
    out.fns=table.annovar(fns.annovar=fns.annovar, output.dir=output.dir, sample.name=sample.name, annovar.db.dir=annovar.db.dir, ref=ref, protocol=protocol, protocol.type=protocol.type, nastring=nastring, run.cmd=run.cmd, mc.cores=mc.cores)
    out.fns
  }
  
  if(annotation.method=="vep"){
    out.fns=vep(fns.vcf=fns.vcf, output.dir=output.dir, vep.db.dir=vep.db.dir, sample.name=sample.name, run.cmd=run.cmd, mc.cores=mc.cores)
    out.fns
  }
  out.fns
}


#' @export
run.rseq.abundance=function(rseq.abundance.method, fns.bam, output.dir, sample.name, cufflinks.thread.number, cufflinks.gtf, ref.gtf, RNAtype, Mode, stranded, idattr, htseq.r, htseq.a, mir.gff, run.cmd, mc.cores, align.dir){
  fns.bam=get.fns(align.dir, idx=".Aligned.sortedByCoord.out.bam$|.sam$|accepted_hits.bam$")
  if(rseq.abundance.method=="cufflinks"){
    if(!sum(basename(dirname(fns.bam))==sample.name)!=0) sample.name=basename(dirname(fns.bam))
    out.fns=cufflinks(fns.bam, output.dir=output.dir, sample.name=sample.name, cufflinks.thread.number=cufflinks.thread.number, cufflinks.gtf=cufflinks.gtf, ref.gtf=ref.gtf, run.cmd=run.cmd, mc.cores=mc.cores)
  }
  
  if(rseq.abundance.method=="htseq"&RNAtype=="mRNA"){
    out.fns=htseq_count(RNAtype=RNAtype, fns.bam, sample.name=sample.name, output.dir=output.dir, Mode=Mode, stranded=stranded, idattr=idattr, htseq.r=htseq.r, htseq.a=htseq.a, ref.gtf=ref.gtf, mir.gff=mir.gff, run.cmd=run.cmd, mc.cores=mc.cores)
    gtf.gr=gtf2gr(ref.gtf, output.dir=output.dir)
    info.fns=htseq.add.info(RNAtype=RNAtype, count.fns=out.fns, output.dir=output.dir, mc.cores=mc.cores)
    info.fns
  }
  
  if(rseq.abundance.method=="htseq"&RNAtype=="miRNA"){
    out.fns=htseq_count(RNAtype=RNAtype, fns.bam, sample.name=sample.name, output.dir=output.dir, Mode=Mode, stranded=stranded, idattr=idattr, htseq.r=htseq.r, htseq.a=htseq.a, ref.gtf=ref.gtf, mir.gff=mir.gff, run.cmd=run.cmd, mc.cores=mc.cores)
    gff.gr=gff2gr(mir.gff, output.dir=output.dir)
    info.fns=htseq.add.info(RNAtype=RNAtype, count.fns=out.fns, output.dir=output.dir, mc.cores=mc.cores)
  }
}


#' @export
run.cnv=function(ref.fa, fns.bam, rmdup.dir, sample.name, mapQ, output.dir, window, cnv.dir, run.cmd, mc.cores){
  
  fns.bam=get.fns(rmdup.dir, idx="rmdu.bam$")
  fn.pileup=read.pileup.gz(ref.fa=ref.fa, fns.bam=fns.bam, sample.name=sample.name, mapQ=mapQ, output.dir=output.dir, run.cmd=run.cmd, mc.cores=mc.cores)
  
  gc.fn=generate.GC(window=window, output.dir=output.dir, ref.fa=ref.fa, run.cmd=run.cmd)
  
  normal.pileup.gz=get.fns(output.dir, idx=n.sample)
  tumor.pileup.gz=get.fns(output.dir, idx=t.sample)
  bin.out.seqz=pileup2seqz(gc.fn=gc.fn, normal.pileup.gz=normal.pileup.gz, window=window, tumor.pileup.gz=tumor.pileup.gz, output.dir=output.dir, run.cmd=run.cmd, mc.cores=mc.cores)
  
  seqz2rda(cnv.dir=cnv.dir)
  seqz2seg(cnv.dir, window=window)
  ploidyNcellularity(cnv.dir)
  
}

#' @export
run.make.vset=function(annot.dir, Robject.dir, vcf.dir, ref.fa, unsafe, minBaseQuality, minMappingQuality, mut.cnt.cutoff, run.cmd, mc.cores){
  make.vset(annot.dir=annot.dir, Robject.dir=Robject.dir, mc.cores=mc.cores)
  gatk.depthOFcoverage(vcf.dir=vcf.dir, annot.dir=annot.dir, Robject.dir=Robject.dir, ref.fa=ref.fa, unsafe=unsafe, minBaseQuality=minBaseQuality, minMappingQuality=minMappingQuality, run.cmd=run.cmd, mc.cores=mc.cores)
  vset.preprocess(Robject.dir=Robject.dir, mut.cnt.cutoff=mut.cnt.cutoff)
}


