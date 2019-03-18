###### dseq-GATK
SEQprocess(project.name="Dseq-GATK-variantcall", 
           fastq.dir="[FastQC directory]", 
           output.dir="[Output directory]",
           config.fn=system.file("data", "WGS.GATK_ex.config.R", package="SEQprocess"), 
           type="WGS", 
           pipeline="GATK",
           mc.cores=1, 
           run.cmd=TRUE)
