

#' @title make.eset
#' @description For the expression data are transformed to a file with extension .eSet
#' @usage make.eset(RNAquant.dir, Robject.dir, mc.cores=1)
#' @param RNAquant.dir Cufflinks or htseq-count output directory
#' @param Robject.dir Output directory 
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @import Biobase
#' @import parallel
#' @details Quantify mRNA and miRNA and store them in ExpressionSet format for convenient analysis
#' @export
make.eset=function(RNAquant.dir, Robject.dir, mc.cores=1){
  
  message("[[",Sys.time(),"]] Start creating an ExpressionSet.---- ")
  
  #### Make eset from cufflinks
  if(length(which(grepl("tracking", dir(RNAquant.dir, recursive=TRUE, full.names=TRUE))))>0){
    
    # Read Cufflinks fpkm tracking files
    cuff=get.fns(input.dir=RNAquant.dir, idx="genes.fpkm_tracking")
    sample.name=basename(dirname(cuff))
    cuff.list=mclapply(cuff, read.delim, mc.cores=mc.cores)
    
    # ID
    id.list=lapply(cuff.list, function(a) paste(a$tracking_id, a$locus, sep="_"))
    uniq.id=as.character(unique(unlist(id.list)))
    
    # list of cufflinks fpkm tracking files
    cuff.list=lapply(1:length(cuff.list), function(a) {rownames(cuff.list[[a]])=id.list[[a]]; cuff.list[[a]][uniq.id,]})
    names(cuff.list)=sample.name
    
    # assayData
    expr=do.call(cbind, lapply(cuff.list, function(a) a$FPKM))
    rownames(expr)=rownames(cuff.list[[1]])
    
    # featureData
    fdata=cuff.list[[1]][,1:9]
  }
  
  ##### make eset from htseq-count, mRNA
  if(length(which(grepl("CountAddInfo.mRNA", dir(RNAquant.dir, recursive=TRUE, full.names=TRUE))))>0){
    # read htseq fpkm text files
    htseq=get.fns(input.dir=RNAquant.dir, idx=".CountAddInfo.mRNA.txt$")
    sample.name=sub(".CountAddInfo.mRNA.txt", "", basename(htseq))
    htseq.list=mclapply(htseq, read.delim, mc.cores=mc.cores)
    
    # ID
    id.list=lapply(htseq.list, function(a) a$id)
    uniq.id=unique(unlist(id.list))
    
    # list of cufflinks fpkm tracking files
    htseq.list=lapply(1:length(htseq.list), function(a) {rownames(htseq.list[[a]])=id.list[[a]]; htseq.list[[a]][uniq.id,]})
    names(htseq.list)=sample.name
    
    # assayData
    expr=do.call(cbind, lapply(htseq.list, function(a) a$FPKM))
    rownames(expr)=rownames(htseq.list[[1]])
    
    # featureData
    fdata=htseq.list[[1]][,1:7]
  }
  
  ##### make eset from htseq-count, miRNA
  if(length(which(grepl("CountAddInfo.miRNA", dir(RNAquant.dir, recursive=TRUE, full.names=TRUE))))>0){
    htseq=get.fns(input.dir=RNAquant.dir, idx=".CountAddInfo.miRNA.txt")
    sample.name=sub(".CountAddInfo.miRNA.txt", "", basename(htseq))
    htseq.list=lapply(htseq, read.delim)
    
    # ID
    id.list=lapply(htseq.list, function(a) a$id)
    uniq.id=unique(unlist(id.list))
    
    # list of cufflinks fpkm tracking files
    htseq.list=lapply(1:length(htseq.list), function(a) {rownames(htseq.list[[a]])=id.list[[a]]; htseq.list[[a]][uniq.id,]})
    names(htseq.list)=sample.name
    
    # assayData
    expr=do.call(cbind, lapply(htseq.list, function(a) a$counts))
    rownames(expr)=rownames(htseq.list[[1]])
    
    # featureData
    fdata=htseq.list[[1]][,1:7]
  }
  
  # phenotypeData
  pdata=data.frame(sampleID=colnames(expr), row.names=colnames(expr))
  library(Biobase)
  metadata=data.frame(labelDescription=colnames(pdata), row.names=colnames(pdata))
  pData=new("AnnotatedDataFrame", data=as.data.frame(pdata), varMetadata=metadata)
  fmetadata=data.frame(labelDescription = colnames(fdata), row.names=colnames(fdata))
  fData=new("AnnotatedDataFrame", data=as.data.frame(fdata), varMetadata=fmetadata)
  
  eset=new("ExpressionSet", exprs=expr, phenoData=pData, featureData=fData)
  eset=eset[which(!rowSums(exprs(eset))==0)]
  save(eset, file=file.path(Robject.dir, "eSet.rda"))
  message("[[",Sys.time(),"]] ExpressionSet is saved in your R object directory---- ")
  eset
}


#' @title refGene2gr
#' @description Convert gene reference file to GRanges form
#' @usage refGene2gr(refGene.path, cnv.dir)
#' @param refGene.path Path to refGene.txt
#' @param cnv.dir Copy number variation directory
#' @export
refGene2gr=function(refGene.path, cnv.dir){
  refGene=read.table(refGene.path)
  colnames(refGene)=c("bin","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","score","name2","cdsStartStat","cdsEndStat","exonFrames")
  refGene.GRanges=with(refGene, GRanges(seqnames=refGene[,3], ranges=IRanges(refGene[,5],refGene[,6]), strand=refGene[,4], refGene[,setdiff(colnames(refGene), c("chrom","strand", "txStart" ,"txEnd"))]))
  save(refGene.GRanges, file=file.path(cnv.dir, "refGene.gr.rda"))
}


#' @title make.cset
#' @description The copy number variants data are transformed to a file with extension .cSet
#' @param cnv.dir sequenza output directory
#' @param Robject.dir Ouptut directory
#' @usage make.cset(cnv.dir, Robject.dir)
#' @export
make.cset=function(cnv.dir, Robject.dir){
  seg=get(load(file.path(cnv.dir, "ploidyNcellularity.rda")))
  
  update.seqnames=function(seqnames){
    seqnames=as.character(seqnames)
    seqnames=sub("chr", "",seqnames, ignore.case=T)
    seqnames=sub("23", "X", seqnames, ignore.case=T)
    seqnames=sub("24", "Y", seqnames, ignore.case=T)
    seqnames=paste0("chr", seqnames)
    seqnames
  }
  
  seg2grList=function(segList, score="depth.ratio", chr="chromosome", start="start.pos", end="end.pos"){
    seg.grList=NULL
    for (i in 1:length(segList)){
      
      seg.i=segList[[i]]
      seg.i=seg.i[which(!is.na(seg.i[,score])),]
      seg.seqnames=update.seqnames(seg.i[,"chromosome"])
      gr=with(seg.i, GRanges(seqnames=seg.seqnames,  ranges=IRanges(seg.i[,start], seg.i[,end]), strand="*", score=seg.i[,score], seg.i[,setdiff(colnames(seg.i), c(chr, start, end, score))]))
      seg.grList <- c(seg.grList, gr)
    }
    names(seg.grList)=names(segList)
    seg.grList
  }
  
  seg.grList=seg2grList(seg, score = "depth.ratio")
  
  getScores=function(interval, signal, score="seg.mean"){
    overlaps <- findOverlaps(interval, signal)
    score(interval) = NA
    score(interval)[queryHits(overlaps)] = data.frame(signal)[,score][subjectHits(overlaps)]
    score(interval)
  }
  
  updateSymbol=function(keys, separator = "[;|\\s+|/+]", unapproved.ids="^LOC|^KIAA|orf|^OTTHUMG|^DKFZ", checkGeneSymbol=F,species = "Hs", rm.llid=T){
    
    run.checkGeneSymbols=function(keys,separator=separator,unmapped.as.na=F){
      library(HGNChelper)
      keys=checkGeneSymbols(keys,unmapped.as.na=unmapped.as.na)[,3]
      keys=get.1stkey(keys,separator=separator)
    }
    
    #for one to many 
    keys.el=get.1stkey(keys,separator=separator)
    if(species=="Hs" & checkGeneSymbol ) keys.el=run.checkGeneSymbols(keys.el,separator=separator)
    
    # handle upapproved ids 
    el = strsplit(as.character(keys), separator, perl = T)
    
    if(rm.llid) {
      library(limma)
      ll.st=grepl(unapproved.ids, keys.el)
      keys.el = alias2SymbolTable(keys.el, species = species)
      ll.st=which(ll.st | is.na(keys.el))
      if(length(ll.st)>0) {
        llid = unlist(lapply(el[ll.st], function(a) a[which(!grepl(unapproved.ids, a))][1]))
        keys.el[ll.st]=llid
      }
      
    }else{
      ll.st=which(is.na(keys.el) |is.null(keys.el) )
      if(species=="Hs" & checkGeneSymbol) keys.el[ll.st]=run.checkGeneSymbols(keys.el[ll.st],unmapped.as.na=TRUE)
    }
    
    return(keys.el)
  }
  
  get.1stkey=function(keys,separator = "[;|\\s+|/+]"){
    keys.el=sub(paste0(separator, ".*"),"",keys)
    return(keys.el)
  }
  
  seg2cset.symbol=function(interval, seg.grList, symbol="name2", score="score"){
    
    # seg.grList to cset
    expr=sapply(1:length(seg.grList),  function(a) getScores(interval = interval, signal = seg.grList[[a]], score=score)) # getScore
    colnames(expr)=names(seg.grList)
    expr.agg=aggregate(expr, list(updateSymbol(as.character(data.frame(interval)[,symbol]))), mean, rm.na=T)
    rownames(expr.agg)=expr.agg$Group.1
    expr.agg=expr.agg[,-1]
    
    # cset
    fdata=data.frame(interval)[match(rownames(expr.agg), as.character(data.frame(interval)[,symbol])),]
    rownames(fdata)=rownames(expr.agg)
    pdata=data.frame(sampleID=colnames(expr.agg), row.names=colnames(expr.agg))
    library(Biobase)
    metadata=data.frame(labelDescription=colnames(pdata), row.names=colnames(pdata))
    pData=new("AnnotatedDataFrame", data=as.data.frame(pdata), varMetadata=metadata)
    fmetadata=data.frame(labelDescription = colnames(fdata), row.names=colnames(fdata))
    fData=new("AnnotatedDataFrame", data=as.data.frame(fdata), varMetadata=fmetadata)
    
    eset=new("ExpressionSet", exprs=expr.agg, phenoData=pData, featureData=fData)
    eset=eset[!rowSums(is.na(exprs(eset)))==length(sampleNames(eset)),]
    
    eset
  }
  
  refGene.gr=get(load(file.path(cnv.dir, "refGene.gr.rda")))
  csetList=mclapply(c("score", "Bf", "A", "B", "CNt"), function(a) seg2cset.symbol(interval = refGene.gr, seg.grList = seg.grList, score = a), mc.cores=5)
  names(csetList) = c("cset.dp", "cset.bf", "cset.cn_A", "cset.cn_B", "cset.cn_t")
  for(i in 1:length(csetList)) {
    csetList[[i]]=csetList[[i]][order(fData(csetList[[i]])$start)]
    csetList[[i]]=csetList[[i]][!is.na(fData(csetList[[i]])$seqnames)]
  }
  
  save(csetList, file=file.path(Robject.dir, "cSetList.rda"))
  message("[[",Sys.time(),"]] ExpressionSet is saved in your R object directory---- ")
  
  csetList
}


#' @title make.vset
#' @description The mutations data are transformed to a file with extension .vSet
#' @param annot.dir ANNOVAR output directory
#' @param Robject.dir Ouptut directory
#' @param mut.cnt.cutoff Depth filter
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @import Biobase
#' @import parallel
#' @export 
make.vset=function(annot.dir, Robject.dir, mc.cores=1){
  
  annot.fns=get.fns(input.dir=annot.dir, idx=".csv$")
  anno=mclapply(annot.fns, read.csv, mc.cores=mc.cores)
  sample.name=basename(dirname(annot.fns))
  
  for(i in 1:length(sample.name)) rownames(anno[[i]]) = paste(anno[[i]]$Chr, anno[[i]]$Start, anno[[i]]$End, anno[[i]]$Ref, anno[[i]]$Alt, sep = "_")
  idList=lapply(anno, rownames)
  uniq.id=unique(unlist(idList))
  
  ## vdat
  vdat=matrix(0, nrow = length(uniq.id), ncol = length(sample.name))
  rownames(vdat)=uniq.id
  colnames(vdat)=sample.name
  
  for(i in 1:length(sample.name)) vdat[rownames(anno[[i]]),i] = 1
  
  ## fdat
  fdata=matrix(nrow=length(uniq.id), ncol=ncol(anno[[1]]))
  rownames(fdata)=uniq.id
  colnames(fdata)=colnames(anno[[1]])
  for(i in 1:length(sample.name)) fdata[idList[[i]],] = as.matrix(anno[[i]][,colnames(fdata)])
  
  # pdat
  pdata=data.frame(sampleID=colnames(vdat), row.names=colnames(vdat))
  library(Biobase)
  metadata=data.frame(labelDescription=colnames(pdata), row.names=colnames(pdata))
  pData=new("AnnotatedDataFrame", data=as.data.frame(pdata), varMetadata=metadata)
  fmetadata=data.frame(labelDescription = colnames(fdata), row.names=colnames(fdata))
  fData=new("AnnotatedDataFrame", data=as.data.frame(fdata), varMetadata=fmetadata)
  
  vset=new("ExpressionSet", exprs=vdat, phenoData=pData, featureData=fData)
  vset=vset[featureNames(vset)[!grepl("KI|random|Un|GL", featureNames(vset))]]
  save(vset, file=file.path(Robject.dir, "vSet.rda"))
  message("[[",Sys.time(),"]] VariantSet is saved in your R object directory---- ")
  vset
}


#' @title gatk.depthOFcoverage
#' @description Calculate the read depth of the position with single nucleotide variations
#' @param vcf.dir Output of variant call step (directory of vcf files)
#' @param annot.dir Output of annotation step (directory of .annovar files) 
#' @param ref.fa Reference fasta file path
#' @param unsafe A parameter value for -U ALLOW_N_CIGAR_READS in GATK. This parameter must be TRUE in RNA-seq data. 
#' @param minBaseQuality Minimum base quality (default=1)
#' @param minMappingQuality Minimum mapping quality (default=1)
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to dedicate. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details When creating a vSet, use read depth to determine whether a mutation exists. GATK DepthOfCoverage uses the interval 
#'          bed file to calculate the depth of the position.
#' @export 
gatk.depthOFcoverage=function(vcf.dir, annot.dir, Robject.dir, ref.fa, unsafe, minBaseQuality=1, minMappingQuality=1, run.cmd=TRUE, mc.cores=1){
  
  # BAM file list
  bam.fns=get.fns(vcf.dir, idx=".recal.bam$")
  bam.list=paste("-I", bam.fns)
  sample.name=sub(bam.idx, "", basename(bam.list))
  
  bam.fns.list=""
  for(i in 1:length(bam.fns)) bam.fns.list=paste(bam.fns.list, bam.list[i])
  
  # annotation csv file 
  csv.fns=get.fns(annot.dir, idx=".csv$")
  message("[[",Sys.time(),"]] Read annovar files --------")
  annot.list=mclapply(csv.fns, read.csv, mc.cores=mc.cores)
  
  # id list
  vset=get(load(file.path(Robject.dir, "vSet.rda")))
  uniq.id=featureNames(vset)
  
  # chr, start, end, ref, alt
  chr=strsplit2(as.character(uniq.id), split="_")[,1]
  start.pos=as.numeric(strsplit2(as.character(uniq.id), split="_")[,2])
  end.pos=as.numeric(strsplit2(as.character(uniq.id), split="_")[,3])
  ref=strsplit2(as.character(uniq.id), split="_")[,4]
  alt=strsplit2(as.character(uniq.id), split="_")[,5]
  
  # variant region -> bed file
  bed=data.frame(chr=as.character(chr), start=as.numeric(start.pos)-1, end=as.numeric(end.pos), ref=as.character(ref), alt=as.character(alt))
  bed.path=file.path(annot.dir, "target.bed")
  message("[[",Sys.time(),"]] List of intervals saved in annot.dir  --------")
  write.table(bed, file=bed.path, quote=F, row.names=FALSE, col.names = F)
  
  out.fns=file.path(annot.dir, "BaseCounts")
  
  # run depthOFcoverage
  unsafe=TRUE
  cmd.add=ifelse(unsafe, "-U ALLOW_N_CIGAR_READS", "")
  cmd=paste("java -jar", GATK.path, "-T DepthOfCoverage", "-R", ref.fa, "-o", out.fns, bam.fns.list, "-L", bed.path, "--printBaseCounts", "--minBaseQuality", minBaseQuality, "--minMappingQuality", minMappingQuality, cmd.add)
  
  message("[[",Sys.time(),"]] Run GATK DepthOfCoverage --------")
  if(run.cmd) system(cmd)
  cat(cmd, file=file.path(annot.dir, "run.depthOFcoverage.log"), sep="\n", append=FALSE)
  message("[[",Sys.time(),"]] Done --------")
  
  # read depthOFcoverage output file
  BaseCounts=read.delim(out.fns)
  
  # BaseCount => mutDepth, wildDepth, allDepth
  BaseCounts=BaseCounts[match(paste(bed$chr, bed$start+1, sep=":"),BaseCounts$Locus),]
  rownames(BaseCounts)=paste(bed$chr, bed$start+1, bed$ref, bed$alt, sep="_")
  BaseCounts$Chr=chr
  BaseCounts$Start=start.pos
  BaseCounts$End=end.pos
  BaseCounts$Ref=ref
  BaseCounts$Alt=alt
  
  Depth=BaseCounts[,grep("base_counts|Depth_for|Chr|Start|End|Ref|Alt", colnames(BaseCounts))]
  
  ATGC=list()
  
  for(i in 1:length(sample.name)){
    ATGC[[i]]=strsplit2(as.character(Depth[,(2*i)]), split = " |:")
  }
  names(ATGC)=sample.name
  
  for(i in 1:length(sample.name)){
    Depth=data.frame(Depth, ATGC[[1]][,c(2,4,6,8)])
    colnames(Depth)[(ncol(Depth)-3)]=paste(as.character(sample.name[i]), "A", sep="_")
    colnames(Depth)[(ncol(Depth)-2)]=paste(as.character(sample.name[i]), "C", sep="_")
    colnames(Depth)[(ncol(Depth)-1)]=paste(as.character(sample.name[i]), "G", sep="_")
    colnames(Depth)[ncol(Depth)]=paste(as.character(sample.name[i]), "T", sep="_")
  }
  
  wtDP=matrix(NA, nrow = nrow(BaseCounts), ncol = length(sample.name))
  for(i in 1:length(sample.name)){
    mat=ATGC[[i]]
    
    wtDP[Depth$Ref=="A",i]=mat[Depth$Ref=="A",2]
    wtDP[Depth$Ref=="C",i]=mat[Depth$Ref=="C",4]
    wtDP[Depth$Ref=="G",i]=mat[Depth$Ref=="G",6]
    wtDP[Depth$Ref=="T",i]=mat[Depth$Ref=="T",8]
    
  }
  
  mutDP=matrix(NA, nrow = nrow(BaseCounts), ncol = length(sample.name))
  for(i in 1:length(sample.name)){
    mat=ATGC[[i]]
    
    mutDP[Depth$Alt=="A",i]=mat[Depth$Alt=="A",2]
    mutDP[Depth$Alt=="C",i]=mat[Depth$Alt=="C",4]
    mutDP[Depth$Alt=="G",i]=mat[Depth$Alt=="G",6]
    mutDP[Depth$Alt=="T",i]=mat[Depth$Alt=="T",8]
  }
  
  for(i in 1:length(sample.name)){
    allDP=Depth[,grep("Depth_for" , colnames(Depth))]
  }
  
  colnames(mutDP)=paste(sample.name, "mutDP", sep="_")
  colnames(wtDP)=paste(sample.name, "wtDP", sep="_")
  colnames(allDP)=paste(sample.name, "allDP", sep="_")
  
  Depth.df=data.frame(mutDP, wtDP, allDP)
  
  fData(vset)=data.frame(fData(vset), Depth.df)
  
  save(vset, file=file.path(Robject.dir, "vSet.rda"))
  message("[[",Sys.time(),"]] Depth information added to vSet --------")
}

#' @title vset.preprocess
#' @description Use the depth of the variants position to determine the presence or absence of the mutation. 
#' @usage vset.preprocess(Robject.dir, mut.cnt.cutoff)
#' @param annot.dir Output of annotation step (directory of .annovar files) 
#' @param Robject.dir Ouptut directory
#' @param mut.cnt.cutoff Criterion of depth (default=8)
#' @param mc.cores The number of cores to dedicate. Must be at least one(default=1), and parallelization requires at least two cores.
#' @import Biobase
#' @details vSet before preprocessing is simply denoted as 1 if there is a mutation and 0 otherwise. After the depth of the mutated position
#'          is calculated, the variants position with a depth that does not meet a certain criterion is not included in the analysis.
#' @export 
vset.preprocess=function(annot.dir, Robject.dir, mut.cnt.cutoff=8, mc.cores=1){
  
  annovar.out=dir(file.path(annot.dir), pattern=".annovar$", recursive=TRUE, full.names=TRUE)
  annoList=mclapply(annovar.out, read.delim, header=FALSE, mc.cores=mc.cores)
  for(i in 1:length(annoList)){
    rownames(annoList[[i]])=paste(annoList[[i]]$V1, annoList[[i]]$V2, annoList[[i]]$V3, as.character(annoList[[i]]$V4), as.character(annoList[[i]]$V5), sep="_")
    annoList[[i]]=annoList[[i]][rownames(annoList[[i]])[!grepl("KI|random|Un|GL", rownames(annoList[[i]]))],]
  }
  
  vset=get(load(file.path(Robject.dir, "vSet.rda")))
  sample.name=sampleNames(vset)
  gatkQC.mat=matrix(NA, nrow=nrow(vset), ncol=ncol(vset))
  rownames(gatkQC.mat)=featureNames(vset)
  colnames(gatkQC.mat)=paste0(sample.name, "_gatkQC")
  for(i in 1:length(sample.name)) gatkQC.mat[rownames(annoList[[i]]),i]=as.character(annoList[[i]]$V12)
  gatkQC.mat=data.frame(gatkQC.mat)
  
  fData(vset)=data.frame(fData(vset), gatkQC.mat)
  head(fData(vset))
  fData(vset)$indel=FALSE
  fData(vset)$indel[nchar(as.character(fData(vset)$Ref))>1 | nchar(as.character(fData(vset)$Alt))>1 | as.character(fData(vset)$Ref)=="-" | as.character(fData(vset)$Alt)=="-"]=TRUE
  
  ## vdat
  gatkQC=fData(vset)[,paste0(sample.name, "_gatkQC")]
  mutDP=fData(vset[!fData(vset)$indel])[,paste0(sample.name, "_mutDP")]
  wtDP=fData(vset[!fData(vset)$indel])[,paste0(sample.name, "_wtDP")]
  mutDP=apply(fData(vset)[,paste0(sample.name, "_mutDP")],2, as.numeric)
  wtDP=apply(fData(vset)[,paste0(sample.name, "_wtDP")],2, as.numeric)
  exprs(vset)[which(!fData(vset)$indel & mutDP>0 & mutDP<mut.cnt.cutoff | !fData(vset)$indel & wtDP<mut.cnt.cutoff & mutDP==0)]=NA
  exprs(vset)[which(!fData(vset)$indel & wtDP>=mut.cnt.cutoff & mutDP==0)]=0
  exprs(vset)[which(!fData(vset)$indel & mutDP>=mut.cnt.cutoff & gatkQC=="PASS")]=1
  
  fil.vset=vset
  fil.vset=fil.vset[!rowSums(is.na(exprs(fil.vset)))==length(sampleNames(fil.vset)),]
  save(fil.vset, file=file.path(Robject.dir, "filter_vSet.rda"))
}  


#' @title eset2SE
#' @description Convert ExpressionSet to SummarizedExperiment
#' @param eset RNA abundance ExpressionSet
#' @param vset Variant ExpressionSet
#' @param cset CNV ExpressionSet
#' @param Robject.dir Output directory
#' @import GenomicRanges 
#' @import SummarizedExperiment 
#' @details SEQprocess also provides SummarizedExperiment data format for comportable data analysis and management
#' @export
eset2SE=function(eset=NULL, vset=NULL, cset=NULL, Robject.dir){
  library(GenomicRanges)
  library(SummarizedExperiment)
  
  if(!is.null(eset)){
    counts=exprs(eset)
    
    chrname=sapply(as.character(fData(eset)$locus), function(a) strsplit(a, split=":")[[1]][1])
    start.pos=as.numeric(sapply(as.character(fData(eset)$locus), function(a) strsplit(a, split=":|-")[[1]][2]))
    end.pos=as.numeric(sapply(as.character(fData(eset)$locus), function(a) strsplit(a, split=":|-")[[1]][3]))
    rowRanges=GRanges(chrname,
                      IRanges(start=start.pos, end=end.pos),
                      feature_id=featureNames(eset), fData(eset))
    colData=DataFrame(pData(eset))
    
    SE=SummarizedExperiment(assays=list(counts=counts),
                            rowRanges=rowRanges, colData=colData)
    save(SE, file=file.path(Robject.dir, "exprSE.rda"))
  }
  
  if(!is.null(vset)){
    counts=exprs(vset)
    
    chrname=as.character(fData(vset)$Chr)
    start.pos=as.numeric(as.character(fData(vset)$Start))
    end.pos=as.numeric(as.character(fData(vset)$End))
    
    rowRanges=GRanges(chrname,
                      IRanges(start=start.pos, end=end.pos),
                      strand="*",
                      feature_id=featureNames(vset), fData(vset))
    colData=DataFrame(pData(vset))
    
    SE=SummarizedExperiment(assays=list(counts=counts),
                            rowRanges=rowRanges, colData=colData)
    save(SE, file=file.path(Robject.dir, "variantSE.rda"))
  }
  
  if(!is.null(cset)){
    counts=lapply(cset, function(a) exprs(a))
    
    chrname=fData(cset[[1]])$seqnames
    start.pos=as.numeric(fData(cset[[1]])$start)
    end.pos=as.numeric(fData(cset[[1]])$end)
    pos.width=as.numeric(fData(cset[[1]])$width)
    colnames(fData(cset[[1]]))[1:5]=toupper(colnames(fData(cset[[1]]))[1:5])
    rowRanges=GRanges(chrname, IRanges(start=start.pos, end=end.pos, width=pos.width), 
                             strand="*", feature_id=featureNames(cset[[1]]), fData(cset[[1]]))
    colData=DataFrame(pData(cset[[1]]))
      
    SE=lapply(1:length(counts), function(a) SummarizedExperiment(assays=list(counts=counts), rowRanges=rowRanges, colData=colData))
    
    save(SE, file=file.path(Robject.dir, "cnvSE.rda"))
  }
}

