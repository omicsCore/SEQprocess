get.option=function(){
  library(optparse)
  option_list = list(
    make_option(c("--fastq.dir"), type="character", default=NULL, 
                help="Directory path of fastq files [default = NULL]", action="store_true"),
    
    make_option(c("--output.dir"), type="character", default=NULL, 
                help="Output directory path [default = Working directory/result/project.name", action="store_true"),
    
    make_option(c("--project.name"), type="character", default=NULL, 
                help="Project name [default = SEQprocess]", action="store_true"),
    
    make_option(c("--type"), type="character", default=NULL, 
                help="Sequencing data type [default = WGS", action="store_true"),
    
    make_option(c("--pipeline"), type="character", default=NULL, 
                help="Select pre-customized six data processing pipeline  [default = none", action="store_true"),
    
    make_option(c("--mc.cores"), type="numeric", default=1, 
                help="Number of cores
                If the number of cpu is more than one, you can change the settings of the file to use the number you want.
                [default = 1]", action="store_true"),
    
    make_option(c("--run.cmd"), type="logical", default=TRUE, 
                help="Whether run SEQprocess [default = TRUE]", action="store_true"),
    
    make_option(c("--report.mode"), type="logical", default=FALSE,
                help="Whether to create a report or not. If report mode option is TRUE, SEQprocess do not run and make report file [default = FALSE]", action="store_true"),
    
    make_option(c("--config.fn"), type="character", default="config.fn",
                help="Configure file contained program path, reference files and index files [default = /data/config.path.R]", action="store_true"),
    
    make_option(c("--qc"), type="logical", default=TRUE, 
                help="Quality check of fastq sequence files [default = TRUE]", action="store_true"),
    
    make_option(c("--trim.method"), type="logical", default=TRUE, 
                help="Select trimming off low quality bases, and cleaning up adapter sequences method. 
                [default = trim.galore]", action="store_true"),
    
    make_option(c("--align.method"), type="character", default="bwa", 
                help="[default] bwa = Mapping fastq files to reference genome using BWA-mem
                bowtie2 = Mapping fastq files to reference genome using bowtie2
                tophat2 = Mapping fastq files to reference genome using bowtie2, tophat2
                star = Mapping fastq files to reference genome using STAR
                none = Do not run any alignment tools", action="store_true"),
    
    make_option(c("--build.transcriptome.idx"), type="logical", default=4, 
                help="A transcriptome index and the associated data files (the original GFF file) can be thus reused for multiple TopHat
                runs with this option, so these files are only created for the first run with a given set of transcripts.
                [default = FALSE]", action="store_true"),
    
    make_option(c("--tophat.thread.number"), type="numeric", default=4, 
                help="A numeric value of the number of threads
                [default = 4]", action="store_true"),
    
    make_option(c("--bwa.method"), type="character", default="mem", 
                help="[default] mem = For paired-end data, 
                aln: For single-end data", action="store_true"),
    
    make_option(c("--bwa.thread.number"), type="numeric", default=4, 
                help="A numeric value of the number of threads
                [default = 4]", action="store_true"),
    
    make_option(c("--star.thread.number"), type="numeric", default=8, 
                help="A numeric value of the number of threads
                [default = 8]", action="store_true"),
    
    make_option(c("--rm.dup"), type="character", default="MarkDuplicates", 
                help="[default] MarkDuplicates = Remove duplicants of bam files.
                BARCODE = Remove duplicate barcodes of Barcode sequencing data
                none = Do not remove duplicate barcodes and reads
                [default = MarkDuplicates]", action="store_true"),
    
    make_option(c("--realign"), type="logical", default=TRUE,
                help="Realign the indel position to increase accuracy when calling variants
                [default = TRUE]", action="store_true"),
    
    make_option(c("--variant.call.method"), type="character", default="gatk", 
                help="[default] gatk = Variant calling from sequencing data using gatk
                varscan2 = Variant call using VarScan2
                mutect2 = Variant call using MuTect2
                muse = Variant call using MuSE
                somaticsniper = Variant call using SomaticSniper
                none = Do not call variants", action="store_true"),
    
    make_option(c("--gatk.thread.number"), type="numeric", default=4,
                help="A numeric value of the number of threads
                [default = 4]", action="store_true"),
    
    make_option(c("--annotation.method"), type="character", default="annovar", 
                help="[default] annovar = Information annotate about variants with chromosome, start position, end position, reference nucleotide and alteration nucleotides.
                vep = Use the variant effect predictor according to the GDC pipeline
                none = Do not annotate information about variants", action="store_true"),
    
    make_option(c("--ref"), type="character", default="hg38", 
                help="Genome reference version [default = hg38]", action="store_true"),
    
    make_option(c("--rseq.abundance.method"), type="character", default="none", 
                help="[default] none = Do not estimate RNA abundant
                cufflinks = Assembles transcriptomes from RNA-Seq data and quantifies their expression.
                htseq = Count how many reads map to each feature and convert reads count to FPKM value.", action="store_true"),
    
    make_option(c("--cufflinks.gtf"), type="character", default="G", 
                help="[defafult] G = Do not include novel genes and isoforms that are assembled.
                g = Include novel genes and isoforms that are assembled", action="store_true"),
    
    make_option(c("--cufflinks.thread.number"), type="numeric", default=4, 
                help="A numeric value of the number of threads
                [default = 4]", action="store_true"),
    
    make_option(c("--RNAtype"), type="character", default="mRNA", 
                help="Select RNA type when running htseq-count
                [default] mRNA
                miRNA", action="store_true"), 
    
    make_option(c("--CNV"), type="logical", default=FALSE, 
                help="DNA copy number estimation using Sequenza R package
                [default = FALSE]" , action="store_true"),
    
    make_option(c("--make.eSet"), type="logical", default=FALSE, 
                help="The processed data can be transformed to an ‘ExpressionSet’ compatible R data type, which can facilitate subsequent data analysis 
                [default = FALSE]", action="store_true"),
    
    make_option(c("--eset2SummarizedExperiment"), type="logical", default=FALSE, 
                help="ExpressionSet is converted to SummarizedExperiment class
                [default = FALSE]", action="store_true"),
    
    make_option(c("--mut.cnt.cutoff"), type="numeric", default=8,
                help="Standard depth to determine the presence of a mutation
                [default = 8]", action="store_true"),
    
    make_option(c("--qc.dir"), type="character", default="./output.dir/00_qc",
                help="Output directory of quality check step", action="store_true"),
    
    make_option(c("--trim.dir"), type="character", default="./output.dir/01_trim",
                help="Output directory of trim step", action="store_true"),
    
    make_option(c("--align.dir"), type="character", default="./output.dir/02_align",
                help="Output directory of alignment step", action="store_true"),
    
    make_option(c("--rmdup.dir"), type="character", default="./output.dir/03_rmdup",
                help="Output directory of remove duplicates step", action="store_true"),
    
    make_option(c("--realign.dir"), type="character", default="./output.dir/04_realign",
                help="Output directory of realignment step", action="store_true"),
    
    make_option(c("--vcf.dir"), type="character", default="./output.dir/05_vcf",
                help="Output directory of variant calling step", action="store_true"),
    
    make_option(c("--annot.dir"), type="character", default="./output.dir/06_annot",
                help="Output directory of variant annotation step", action="store_true"),
    
    make_option(c("--RNAquant.dir"), type="character", default="./output.dir/07_RNAquant",
                help="Output directory of RNA quantitation step", action="store_true"),
    
    make_option(c("--cnv.dir"), type="character", default="./output.dir/08_cnv",
                help="Output directory of copy number estimation step", action="store_true"),
    
    make_option(c("--Robject.dir"), type="character", default="./output.dir/09_Robject",
                help="The directory in which to store the R data type", action="store_true")
    )
  option_list=option_list
  opt_parser=OptionParser(usage="usage: %prog --project.name [Project name] --fastq.dir [fastq files path] --output.dir [Output file path] --type [Data type] --pipeline [Select pre-customized pipeline]......", option_list=option_list)
  opt = parse_args(opt_parser)
}


