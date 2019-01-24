

# WGS GATK Variant call pipeline example
library(parallel)

##################################
#### configure paths & names #####
##################################
fq1.idx=".1.fastq$|_R1.fastq$|.1_val_1.fq$|_1.fq$|.1_val_1.fq$|.R1_val_1.fq$"
fq2.idx=".2.fastq$|_R2.fastq$|.2_val_2.fq$|_2.fq$|.2_val_2.fq$|.R2_val_2.fq$"
bam.idx=".rg.od.bam$|.rmdu.bam$|.realign.bam$|.recal.bam$|.Aligned.sortedByCoord.out.bam$|.sam$|accepted_hits.bam$"
vcf.idx=".f.vcf$|.muse_variants.vcf$|.snp.vcf$|.indel.vcf$|_variants.vcf$"


##################################
######### 1.QualtiyCheck #########
##################################
# FASTQC
fastqc.path="/data/program/bin/fastqc"


##################################
########### 2.Trimming ###########
##################################
# Trim_galore
trim_galore.path="/data/program/bin/trim_galore"
trim.quality=30
trim.clip_R1=13
trim.clip_R2=13

# Cutadapt
cutadapt.path="/data/program/bin/cutadapt"
m=17 # minimum read length
adapt.seq="TGGAATTCTCGGGTGCCAAGG"


##################################
########## 3. Alignment ##########
##################################
# BWA-mem & BWA-aln
bwa.path="/data/program/bin/bwa"

# samtools 
samtools.path="/data/program/bin/samtools"

# Picard AddOrReplaceReadGroups & Picard ReorderSam
picard.path="/data/program/picard/picard/build/libs/picard.jar"
RGLB="LC"
RGPL="Illumina"
RGPU="runbarcode"
SORT_ORDER="coordinate"
VALIDATION_STRINGENCY="LENIENT"
ALLOW_INCOMPLETE_DICT_CONCORDANCE=FALSE 
ALLOW_CONTIG_LENGTH_DISCORDANCE=FALSE
CREATE_INDEX=TRUE


##################################
###### 4. Remove Duplicates ######
##################################
# Picard MarkDupliecates
REMOVE_DUPLICATES=TRUE


##################################
######## 5. Re-Alignment #########
##################################
# GATK
GATK.path="/data/program/gatk/3.7/GenomeAnalysisTK.jar"

# HaplotypeCaller
genotyping_mode="DISCOVERY"
output_mode="EMIT_VARIANTS_ONLY"
stand_call_conf_number=30
FS=30.0
QD=2.0
FS=trimws(format(round(FS, 1), nsmall = 1))
QD=trimws(format(round(QD, 1), nsmall = 1))
QUAL=50
DP=5
gatk.window=35
cluster=3


##################################
######## 7. Annotation ###########
##################################
# ANNOVAR
annovar.db.dir=system.file("extdata", "humandb", package="SEQprocess")
vcf2annovar.pl="/data/program/annovar/annovar/convert2annovar.pl"
format = "vcf4"
coverage = 0

table_annovar.pl="/data/program/annovar/annovar/table_annovar.pl"
ref = "hg38"
protocol = "knownGene,cytoBand,exac03"
protocol.type = "g,r,f"
nastring = "."


##################################
###### 9. Reference FASTA ########
##################################
# Version = GRCH38
ref.fa=system.file("extdata", "reference.fa", package="SEQprocess")


##################################
###### 10. Reference index #######
##################################
# bwa
bwa.idx=system.file("extdata", "reference.fa", package="SEQprocess")


##################################
####### 12. Reference VCF ########
##################################
# Variant Calling Format
ref.gold_indels=system.file("extdata", "indel.vcf", package="SEQprocess")
cosmic.vcf=system.file("extdata", "cosmic.vcf", package="SEQprocess")
ref.dbSNP=system.file("extdata", "snp.vcf", package="SEQprocess")
