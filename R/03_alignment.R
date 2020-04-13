
#' @title bwa
#' @description A wrapper function to run BWA.
#' @usage bwa(bwa.method, fq1, fq2, output.dir, sample.name, ref.fa, bwa.idx, bwa_thread_number=4, run.cmd=TRUE, mc.cores=1)
#' @param bwa.method bwa algorithms of mem and aln can be used(mem: for paired-end data, aln: for single-end data)
#' @param fq1 Path to read1 fastq files
#' @param fq2 Path to read2 fastq files (bwa-mem only)
#' @param output.dir Output directory
#' @param sample.name A character vector for the sample names
#' @param ref.fa Path to reference fasta file
#' @param bwa.idx Path to bwa index files
#' @param bwa.thread.number A parameter value for -t in BWA. A numeric value of the number of threads (default: 4)
#' @param run.cmd  Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details BWA is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome. 
#'          "bwa" can be run with option either of BWA-mem or BWA-aln.
#' @return Aligned BAM files
#' @import parallel
#' @references Fast and accurate short read alignment with Burrows-Wheeler transform
#' @seealso \url{http://bio-bwa.sourceforge.net/bwa.html}
#' @export
bwa=function(bwa.method=c("mem", "aln"),
             fq1, fq2,
             output.dir,
             sample.name,
             ref.fa,
             bwa.idx,
             bwa.thread.number=4,
             run.cmd=TRUE,
             mc.cores=1){
  
  out.dirs=file.path(output.dir, sample.name)
  lapply(1:length(out.dirs), function(a) dir.create(out.dirs[a], recursive = TRUE, showWarnings=FALSE))
  fn.bam=file.path(out.dirs, paste0(sample.name, ".bam"))
  fn.sai=file.path(out.dirs, paste0(sample.name, ".sai"))
  fn.sam=file.path(out.dirs, paste0(sample.name, ".sam"))
  #mem
  
  if(bwa.method=="mem"){
    
    cmd=paste(bwa.path, "mem", "-t", bwa.thread.number, bwa.idx, fq1, fq2, "|", samtools.path, "view", "-bS", "-", "-o", fn.bam)
    
    message("[[",Sys.time(),"]] Run bwa mem---- ")
    print_message(cmd)
    
    if(run.cmd) mclapply(cmd, system,mc.cores=mc.cores)
    cat(cmd, file=file.path(output.dir, "run.bwa-mem.log"), sep="\n", append = FALSE)
    
  }
  
  #aln
  if(bwa.method=="aln"){
    fq2=""
    cmd = paste0(bwa.path, " aln ", "-t ", bwa.thread.number, " -0 ", ref.fa, " ", fq1, " > ", fn.sai)
    
    # run
    message("[[",Sys.time(),"]] Run bwa aln---- ")
    print_message(cmd)
    
    if(run.cmd) mclapply(cmd, system,mc.cores=mc.cores)
    cat(cmd, file=file.path(output.dir, "run.bwa-aln.log"), sep="\n")
    
    # bwa samse
    cmd= paste0(bwa.path, " samse ", ref.fa, " ", fn.sai, " ", fq1, " > ", fn.sam)
    message("[[",Sys.time(),"]] Run bwa samse---- ")
    message(cmd)
    if(run.cmd) mclapply(cmd, system,mc.cores=mc.cores)
    cat(cmd, file=file.path(output.dir, "run.bwa-samse.log"), sep="\n", append = FALSE)
  }
  
  dir(out.dirs, pattern="bam$|sam$", recursive=TRUE, full.names=TRUE)
}



#' @title tophat2
#' @description A wrapper function to run tophat2.
#' @usage tophat2(fq1, fq2, output.dir, sample.name, ref.gtf, bowtie.idx, tophat_thread_number=4, build.transcriptome.idx=FALSE, run.cmd=TRUE, mc.cores=1)
#' @param fq1 Path to read1 fastq files
#' @param fq2 Path to read2 fastq files
#' @param output.dir Output directory
#' @param sample.name A character vector for the sample names
#' @param ref.gtf Path to reference gtf file
#' @param bowtie.idx Path to directory with bowtie indexes and a prefix for the bowtie indexes
#' @param tophat_thread_number A parameter value for -p in tophat2. A numeric value of the number of threads (default: 4) 
#' @param build.transcriptome.idx A parameter value for --transcriptome-index in tophat2. A transcriptome index and the associated data 
#'                                files (the original GFF file) can be thus reused for multiple TopHat runs with this option, so these 
#'                                files are only created for the first run with a given set of transcripts. (default=FALSE)
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details TopHat is a program that aligns RNA-Seq reads to a genome to identify exon-exon splice junctions. It is built on 
#'          the ultrafast short read mapping program Bowtie.
#' @return Aligned BAM files
#' @import parallel
#' @references TopHat2: accurate alignment of transcriptomes in the presence of insertions, deletions and gene fusions
#' @seealso \url{https://ccb.jhu.edu/software/tophat/manual.shtml}
#' @export
tophat2=function(fq1, 
                 fq2, 
                 output.dir,
                 sample.name, 
                 ref.gtf, 
                 bowtie.idx,
                 tophat.thread.number=4, 
                 build.transcriptome.idx=FALSE,
                 run.cmd=TRUE, 
                 mc.cores=1){
  
  #filepath
  out.dirs=file.path(output.dir, sample.name)
  
  #build transcriptome index file 
  if(build.transcriptome.idx){
    message("[[",Sys.time(),"]] Build tophat transcript ---- ")
    cmd.trancscript.build=paste0(tophat2.path,  " --b2-very-sensitive -p ", tophat.thread.number, " -G ", ref.gtf," --transcriptome-index=", transcriptome.idx, " ", 
                                 "-o", out.dirs, bowtie.idx, " ", fq1[1], " ",fq2[1])
    
    system(cmd.trancscript.build)
    
    cmd.trancscript.build=paste0(tophat2.path,  " --b2-very-sensitive -p ", tophat.thread.number, " -G ", ref.gtf," --transcriptome-index=", transcriptome.idx, " ",
                                 "-o", out.dirs, bowtie.idx, " ", fq1[-1], " ", fq2[-1])
    print_message(cmd.trancscript.build)
    
    mclapply(cmd.trancscript.build, system, mc.cores=mc.cores)
  }
  
  #command line
  if(build.transcriptome.idx==FALSE){
    message("[[",Sys.time(),"]] Run tophat---- ")
    cmd = paste(tophat2.path, "--b2-very-sensitive --no-coverage-search -p", tophat.thread.number, "-G", ref.gtf, "-o", out.dirs, bowtie.idx, fq1, fq2)
    
    print_message(cmd)
    if(run.cmd) mclapply(cmd, system, mc.cores=mc.cores)
    cat(cmd, file=file.path(output.dir, "tophat.run.log"), sep="\n", append = FALSE)
  }
  
  print_message(cmd)
  
  out.fns = list.files(out.dirs, "accepted_hits.bam$", full.names = TRUE)
  out.fns
}




#' @title bowtie2
#' @description A wrapper function to run bowtie2.
#' @usage bowtie2(fq1, output.dir, sample.name, bowtie.idx, mc.cores=1, run.cmd=TRUE)
#' @param fq1 Path to read1 fastq files
#' @param output.dir Output directory
#' @param sample.name A character vector for the sample names
#' @param bowtie.idx Path to bowtie index files
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @details Bowtie 2 is an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences. 
#'          Bowtie2 is used only for single-end sequencing data.
#' @return Aligned SAM files
#' @import parallel
#' @references Fast gapped-read alignment with Bowtie 2
#' @seealso \url{http://bowtie-bio.sourceforge.net/bowtie2/index.shtml}
#' @export
bowtie2=function(fq1, 
                 output.dir, 
                 sample.name,
                 bowtie.idx, 
                 mc.cores=1, 
                 run.cmd=TRUE){
  
  out.dirs=file.path(output.dir, sample.name)
  lapply(1:length(out.dirs), function(a) dir.create(out.dirs[a], recursive = TRUE, showWarnings=FALSE))
  
  fn.sam=file.path(out.dirs, paste0(sample.name, ".sam"))
  
  cmd=paste0(bowtie2.path, " -x ", bowtie.idx, " -U ", fq1, " -S ", fn.sam)
  
  message("[[",Sys.time(),"]] Alignment start using bowtie2---- ")
  print_message(cmd)
  if(run.cmd) mclapply(cmd, system,mc.cores=mc.cores)
  cat(cmd, file=file.path(output.dir, "run.bowtie2.log"), sep="\n", append = FALSE)
  
  fn.sam
}


#' @title STAR
#' @description A wrapper function to run STAR (runMode genomeGenerate)
#' @usage build.star.idx(star.idx.dir=file.path(reference.dir, "STAR.idx"), sample.name, align.dir, ref.fa, ref.gtf, sjdbOverhang=100, star_thread_number=8, fasta.idx=TRUE, SJ.idx=FALSE, run.cmd=TRUE)
#' @param star.idx.dir Directory of STAR index files
#' @param sample.name A character vector for the sample names
#' @param ref.fa Reference fasta file path
#' @param ref.gtf Reference gtf file path (e.g., gencode.gtf)
#' @param sjdbOverhang A parameter value for the --sjdbOverhang in STAR. Length of the donor/acceptor sequence on each side of the junctions, ideally=(mate_length-1) (default=100)
#' @param star_thread_number A parameter value for --runThreadN in STAR. A numeric value of the number of threads (default: 8) 
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @param fasta.idx Indexing reference fasta file (when first indexing => TRUE)
#' @param SJ.idx Indexing splicing junction (when second indexing => TRUE)
#' @details Indexing reference fasta file and splicing junction site from fastq files.
#' @return STAR reference and splicing junction indexing
#' @import parallel
#' @references STAR: ultrafast universal RNA-seq aligner
#' @seealso \url {https://github.com/alexdobin/STAR}
#' @export
build.star.idx=function(star.idx.dir=file.path(reference.dir, "STAR.idx"),
                        sample.name,
                        align.dir,
                        ref.fa,
                        ref.gtf,
                        sjdbOverhang=100,
                        star.thread.number=8,
                        fasta.idx=TRUE,
                        SJ.idx=FALSE,
                        run.cmd=TRUE){
  
  dir.create(star.idx.dir, recursive=TRUE, showWarnings=FALSE)
  tmp.dir=file.path(star.idx.dir, "tmp")
  SJ.out.tab=file.path(align.dir, sample.name ,paste0(sample.name, ".SJ.out.tab"))
  idx1=ifelse(fasta.idx, paste("--sjdbGTFfile", ref.gtf), "")
  idx2=ifelse(SJ.idx, paste("--sjdbFileChrStartEnd", SJ.out.tab), "")
  cmd=paste(STAR.path, "--runMode genomeGenerate", "--genomeDir", star.idx.dir, "--genomeFastaFiles", ref.fa, "--sjdbOverhang", sjdbOverhang, "--runThreadN", star.thread.number, "--outTmpDir", tmp.dir, idx1, idx2)
  print_message(cmd)
  
  message("[[",Sys.time(),"]] Build STAR index----")
  if(run.cmd) system(cmd)
  cat(cmd, file.path(star.idx.dir, "run.staridx.log"), append = FALSE)
  
  star.idx.dir
}



#' @title STAR
#' @description  A wrapper function to run STAR.
#' @usage STAR(STAR.idx, fq1, fq2, sample.name, output.dir, run.cmd=TRUE, mc.cores=1)
#' @param star.idx.dir Directory of STAR index
#' @param output.dir Output directory
#' @param sample.name A character vector for the sample names
#' @param fq1 path to read1 fastq files 
#' @param fq2 path to read2 fastq files
#' @param star_thread_number A parameter value for --runThreadN in STAR. A numeric value of the number of threads (default: 8) 
#' @param outFilterMultimapScoreRange A parameter value for --outFilterMultimapScoreRange in STAR. The score range below the maximum score for multimapping alignments (default: 1) 
#' @param outFilterMultimapNmax A parameter value for --outFilterMultimapNmax in STAR. Read alignments will be output only if the read maps fewer than this value, otherwise no alignments will be output (default: 20) 
#' @param outFilterMismatchNmax A parameter value for --outFilterMismatchNmax in STAR. Alignment will be output only if it has fewer mismatches than this value (default: 10) 
#' @param alignIntronMax A parameter value for --alignIntronMax in STAR. Maximum intron length (default: 500,000)
#' @param alignMatesGapMax A parameter value for --alignMatesGapMax in STAR. Maximum genomic distance between mates (default: 1,000,000)
#' @param sjdbScore A parameter value for --sjdbScore in STAR. Extra alignment score for alignmets that cross database junctions (default: 2)
#' @param alignSJDBoverhangMin  A parameter value for --alignSJDBoverhangMin in STAR. Minimum overhang for annotated junctions (default: 1)
#' @param outFilterMatchNminOverLread A parameter value for --outFilterMatchNminOverLread in STAR. Float: outFilterMatchNmin normalized to read length (sum of mates’ lengths for paired-end reads) (default: 0.33)
#' @param outFilterScoreMinOverLread A parameter value for --outFilterScoreMinOverLread in STAR. Float: outFilterScoreMin normalized to read length (sum of mates’ lengths for paired-end reads) (default: 0.33)
#' @param sjdbOverhang A parameter value for --sjdbOverhang in STAR. >=0: Length of the donor/acceptor sequence on each side of the junctions, if =0, splice junction database is not used (default: 100)
#' @param SJ.detect First align, detection of splicing junction (default=TRUE)
#' @param SJ.align Second align, mapping reads to fastq files (default=FALSE)
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details Spliced Transcripts Alignment to a Reference (STAR), which was designed to specifically address many of the challenges of 
#'          RNA-seq data mapping, and uses a novel strategy for spliced alignments.
#' @return Aligned BAM files
#' @import parallel
#' @references STAR: ultrafast universal RNA-seq aligner
#' @seealso \url {https://github.com/alexdobin/STAR}
#' @export
STAR=function(star.idx.dir=file.path(reference.dir, "STAR.idx"),
              output.dir,
              sample.name,
              fq1, 
              fq2, 
              star.thread.number=8, 
              outFilterMultimapScoreRange=1, 
              outFilterMultimapNmax=20, 
              outFilterMismatchNmax=10,
              alignIntronMax=500000, 
              alignMatesGapMax=1000000,
              sjdbScore=2, 
              alignSJDBoverhangMin=1, 
              outFilterMatchNminOverLread=0.33, 
              outFilterScoreMinOverLread=0.33, 
              sjdbOverhang=100,
              SJ.detect=TRUE, SJ.align=FALSE,
              run.cmd=TRUE, mc.cores=1){
  
  out.dirs=file.path(output.dir, sample.name)
  lapply(1:length(out.dirs), function(a) dir.create(out.dirs[a], recursive = TRUE, showWarnings=FALSE))
  prefix=file.path(out.dirs, paste0(sample.name, "."))
  
  fq1.list=NULL
  fq2.list=NULL
  
  for(i in 1:length(fq1)){
    if(i==1){
      fq1.list=fq1[1]
    }else fq1.list=paste(fq1.list,fq1[i],sep=",")
  }
  
  for(i in 1:length(fq1)){
    if(i==1){
      fq2.list=fq2[1]
    }else fq2.list=paste(fq2.list,fq1[i],sep=",")
  }
  
  fastq.list=paste(fq1, fq2, sep=" ")
  pre.align=ifelse(SJ.detect, paste("--outSAMtype None --outSAMmode None"), "")
  post.align=ifelse(SJ.align, paste("--limitBAMsortRAM 0 --outSAMattributes NH HI NM MD AS XS --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate --outSAMheaderHD @HD VN:1.4"), "")
  
#  cmd=paste(STAR.path, 
#            "--genomeDir", star.idx.dir,
#            "--readFilesIn", fastq.list,
#            "--runThreadN", star.thread.number,
#            "--outFileNamePrefix", prefix,
#            "--outFilterMultimapScoreRange", outFilterMultimapScoreRange,
#            "--outFilterMultimapNmax", outFilterMultimapNmax,
#            "--outFilterMismatchNmax", outFilterMismatchNmax,
#            "--alignIntronMax", alignIntronMax,
#            "--alignMatesGapMax", alignMatesGapMax, 
#            "--sjdbScore", sjdbScore,
#            "--alignSJDBoverhangMin", alignSJDBoverhangMin,
#            "--genomeLoad NoSharedMemory",
#            "--outFilterMatchNminOverLread", outFilterMatchNminOverLread,
#            "--outFilterScoreMinOverLread", outFilterScoreMinOverLread, 
#            "--sjdbOverhang", sjdbOverhang,
#            pre.align, post.align)
  
  cmd=paste(STAR.path, 
            "--genomeDir", star.idx.dir,
            "--readFilesIn", fastq.list,
            "--runThreadN", star.thread.number,
            "--outFileNamePrefix", prefix,
            "--genomeLoad NoSharedMemory",
            "--sjdbOverhang", sjdbOverhang,
            pre.align, post.align)
  
  message("[[",Sys.time(),"]] run STAR ---- ")
  print_message(cmd)
  if(run.cmd) mclapply(cmd, system, mc.cores=mc.cores)
  cat(cmd, file=file.path(output.dir, "run.star.log"), sep = "\n", append=FALSE)
  
  if(SJ.detect) out.fns=dir(output.dir, pattern="SJ.out.tab", full.names=TRUE) else if(SJ.align) out.fns=dir(output.dir, pattern="Aligned.sorotedByCoord", full.names=TRUE)
  out.fns
}




#' @title picard.addrg
#' @description A wrapper function to run Picard (AddOrReplaceReadGroups)
#' @usage picard.addrg(fns.bam, output.dir, sample.name, RGLB="lC", RGPL="Illumina", RGPU="runbarcode", RGSM, SORT_ORDER="coordinate", VALIDATION_STRINGENCY="LENIENT", CREATE_INDEX=TRUE, tmp.dir=file.path(output.dir, "tmp"), run.cmd=TRUE, mc.cores=1)
#' @param fns.bam Path to BAM files
#' @param output.dir Output directory
#' @param sample.name A character vector for the sample names
#' @param RGLB A parameter value for RGLB in picard. A character value of Read Group library (default="LC")
#' @param RGPL A parameter value for RGPL in picard. A character value of Read Group platform (default="Illumina")
#' @param RGPU A parameter value for RGPU in picard. A character value of Read Group platform unit (default="runbarcode")
#' @param RGSM A parameter value for RGSM in picard.  character value of Read Group sample name (default=sample.name)
#' @param SORT_ORDER A parameter value for SO in picard. Sort order, a character value of sorting method (default="coordinate")
#' @param VALIDATION_STRINGENCY A parameter value for VALIDATION_STRINGENCY in picard. A character value of validation stringency (default="LENIENT")
#' @param CREATE_INDEX A parameter value for CREATE_INDEX in picard. A character value whether to create .bam index files (default="true")
#' @param tmp.dir Temporary directory path (default=./tmp)
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details This tool enables the user to replace all read groups in the INPUT file with a single new read group and assign 
#'          all reads to this read group in the output BAM files.
#' @return BAM files added read groups
#' @import parallel
#' @seealso \url{http://broadinstitute.github.io/picard/}
#' @export
picard.addrg=function(fns.bam,
                      output.dir,
                      RGLB="LC", RGPL="Illumina", RGPU="runbarcode", RGSM, 
                      SORT_ORDER="coordinate",
                      VALIDATION_STRINGENCY="LENIENT",
                      CREATE_INDEX=TRUE,
                      tmp.dir=file.path(output.dir, "tmp"),
                      run.cmd=TRUE,
                      mc.cores=1){
  
  sample.name=sub(bam.idx,"", basename(fns.bam))
  out.dirs=file.path(output.dir, sample.name)
  out.fns=sub(".bam$", ".rg.bam", fns.bam)
  
  # command 
  cmd=paste0("java -jar ", picard.path, " AddOrReplaceReadGroups", " I=", fns.bam, " O=", out.fns, " SO=", SORT_ORDER, 
             " RGLB=", RGLB, " RGPL=", RGPL, " RGPU=", RGPU, " RGSM=", RGSM, " VALIDATION_STRINGENCY=", VALIDATION_STRINGENCY,
             " CREATE_INDEX=", CREATE_INDEX, " TMP_DIR=", tmp.dir)
  
  # run
  message("[[",Sys.time(),"]] Run picard addRG---- ")
  print_message(cmd)
  if(run.cmd) mclapply(cmd, system, mc.cores=mc.cores) # run
  cat(cmd, file=file.path(output.dir, "run.picard.addrg.log"), sep="\n", append = FALSE) # log
  out.fns
}




#' @title picard.reorder
#' @description A wrapper function to run Picard (ReordrSam)
#' @usage picard.reorder(fns.bam, output.dir, sample.name, ref.fa, ALLOW_INCOMPLETE_DICT_CONCORDANCE=FALSE, ALLOW_CONTIG_LENGTH_DISCORDANCE=FALSE, CREATE_INDEX=TRUE, run.cmd=TRUE, mc.cores=1)
#' @param fns.bam Path to BAM files
#' @param output.dir Output directory
#' @param sample.name A character vector for the sample names
#' @param ref.fa Reference fasta file (eg. GRCh38.fa)
#' @param ALLOW_INCOMPLETE_DICT_CONCORDANCE  A parameter value for ALLOW_INCOMPLETE_DICT_CONCORDANCE in picard. Logical, allow discordant contig (default=FALSE)
#' @param ALLOW_CONTIG_LENGTH_DISCORDANCE  A parameter value for ALLOW_CONTIG_LENGTH_DISCORDANCE in picard. Logical, allow contig of different length (default=FALSE)
#' @param CREATE_INDEX  A parameter value for CREATE_INDEX in picard. Logical, whether to create .bam index files (default=TRUE)
#' @param tmp.dir Temporary directory (default= ./tmp)
#' @param run.cmd Whether to execute the command line (default=TRUE)
#' @param mc.cores The number of cores to use. Must be at least one(default=1), and parallelization requires at least two cores.
#' @details ReorderSam reorders reads in a BAM file to match the contig ordering in a provided reference file, as determined by exact
#'          name matching of contigs. Reads mapped to contigs but absent in the new reference are dropped. Runs substantially faster if the input 
#'          is an indexed BAM file.
#' @return Reordered BAM files (e.g., .rg.od.bam)
#' @import parallel
#' @seealso \url{http://broadinstitute.github.io/picard/}
#' @export
picard.reorder= function(fns.bam, 
                         output.dir, 
                         ref.fa,
                         ALLOW_INCOMPLETE_DICT_CONCORDANCE = FALSE, 
                         ALLOW_CONTIG_LENGTH_DISCORDANCE = FALSE, 
                         CREATE_INDEX = TRUE,
                         tmp.dir = file.path(output.dir, "tmp"),
                         run.cmd = TRUE,
                         mc.cores=1){
  sample.name=sub(bam.idx, "", basename(fns.bam))
  out.dirs=file.path(output.dir, sample.name)
  out.fns = sub(".bam$", ".od.bam", fns.bam)
  
  # command
  cmd = paste0("java -jar ", picard.path, " ReorderSam", " INPUT=", fns.bam, " OUTPUT=", out.fns, " REFERENCE=", ref.fa,  
               " ALLOW_INCOMPLETE_DICT_CONCORDANCE=", ALLOW_INCOMPLETE_DICT_CONCORDANCE, 
               " ALLOW_CONTIG_LENGTH_DISCORDANCE=", ALLOW_CONTIG_LENGTH_DISCORDANCE,
               " CREATE_INDEX=", CREATE_INDEX,
               " TMP_DIR=", tmp.dir)
  
  # run
  message("[[",Sys.time(),"]] Reorder Sam---- ")
  print_message(cmd)
  if(run.cmd) mclapply(cmd, system, mc.cores=mc.cores)
  cat(cmd, file=file.path(output.dir, "run.picard.reorder.log"), sep = "\n", append = FALSE)
  
  out.fns
}
