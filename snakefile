# Snakemake htseq-counts pipeline for unpaired gene
# snakemake --cores 10 --use-conda --resources mem_mb=20000 --retries 3

# Important information
#   Pipeline is engineered for GRCh38 reference genome
#   Input files have the following nomenclature:
#   fastq_files : {sample_name}.fastq.gz


configfile: "workflow/config.yaml"

# Get working directory
from os import getcwd
working_directory = getcwd()

# Get prefix of fastq files
from os import listdir
from sys import exit
sample_names = []
for name in listdir(config["fastq_dir"]):
    if ".fastq.gz" in name:
        
         sample_names.append(name[:-9])

    else:
        exit()

#sample_names = list(set(sample_names))

rule all:
    input:
        "multiqc_report.html"

rule fastp:
    input:
        fastq_files=expand("fastq/{sample}.fastq.gz", sample=sample_names)
 
    output:
        clean_fastq=expand("clean_fastq/clean_{sample}.fastq.gz", sample=sample_names),
        report=expand("results/quality_control/{sample}.fastp.html", sample=sample_names),
        json=expand("results/quality_control/{sample}.fastp.json", sample=sample_names)
    conda:
        "envs/fastp.yaml"
    threads:
        config["threads"]
    shell:
        "fastp -i {input.fastq_files} "
        "-o {output.clean_fastq} "
        "-h {output.report} -j {output.json} "
        "--detect_adapter_for_pe --qualified_quality_phred 15 "
        "--unqualified_percent_limit 40 --length_required 15 "
        "--thread {threads}"

rule bwa:
    input:
        clean_fastq=expand("clean_fastq/clean_{sample}.fastq.gz", sample=sample_names)
    params:
        genome="resources/genome/GRCh38.primary_assembly.genome.fa.gz"
    output:
        temp("sam/{sample}.sam")
    conda:
        "envs/bwa.yaml"
    threads:
        config["threads"]
    shell:
        "bwa mem -t {threads} {params.genome} {input.clean_fastq} > "
        "{output}"


rule sortsam:
    input:
        "sam/{sample}.sam"
    output:
        temp("bam/{sample}.sorted.bam")
    conda:
        "envs/samtools.yaml"
    threads:
         config["threads"]
    shell:
        "samtools sort -l 0 -o {output} -O BAM --threads {threads} {input}" 

rule markduplicates:
    input:
        "bam/{sample}.sorted.bam"
    output:
        bam="bam/{sample}.sorted.removed_duplicates.bam",
        metrics="results/quality_control/markduplicates/{sample}.txt"
    conda:
        "envs/gatk4.yaml"
    shell:
        "gatk MarkDuplicates -I {input} -O {output.bam} "
        "-M {output.metrics} --REMOVE_DUPLICATES "
        "--REMOVE_SEQUENCING_DUPLICATES"

rule bamindex:
    input:
        "bam/{sample}.sorted.removed_duplicates.bam"
    output:
        "bam/{sample}.sorted.removed_duplicates.bam.bai"
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools index {input} {output}"

rule htseq:
    input:
        bam="bam/{sample}.sorted.removed_duplicates.bam.bai",
        annotation="resources/genome/gencode.v43.annotation.gtf.gz"
    output:
        "counts/{sample}_counts.txt"
    conda:
        "envs/htseq.yaml"
    threads: 4
    shell:
         "gunzip -c {input.annotation} > tempfile; "
         "htseq-count -f bam -r pos -s no -i gene_id -p 4 {input.bam} tempfile > tempfile2; "
         "echo '{input.bam}' | tr ' ' '\t' | cat - tempfile2 > {output}"

 
rule multiqc:
    input:
        expand("counts/{sample}_counts.txt", sample=sample_names)
    output:
        "multiqc_report.html"
    conda:
        "envs/multiqc.yaml"
    shell:
        "multiqc --force --ignore resources ."

onsuccess:
    print("Workflow finished, no error")
