import os
import glob

SRA,FRR = glob_wildcards("rawReads/{sra}_{frr}.fastq")

rule all:
    input:
	expand("rawQC/{sra}_{frr}_fastqc.{extension}",sra=SRA, frr=FRR, extension=["zip","html"]),
	expand("starAligned/{sra}Aligned.sortedByCoord.out.bam",sra=SRA)

rule rawFastqc:
	"""
	Quality control trough FASTQC
	"""
    input:
        rawread="rawReads/{sra}_{frr}.fastq"
    output:
        zip="rawQC/{sra}_{frr}_fastqc.zip",
	html="rawQC/{sra}_{frr}_fastqc.html"
    threads:
        1
	params:
		path="rawQC/"
	shell:
        """
	fastqc {input.rawread} --threads {threads} -o {params.path}
        """

rule trimmomatic:
	"""
	Trimming with Trimmomatic
	"""
    input:
        read1="rawReads/{sra}_1.fastq",
	read2="rawReads/{sra}_2.fastq"
    output:
        forwardPaired="trimmedReads/{sra}_1P.fastq",
	reversePaired="trimmedReads/{sra}_2P.fastq"
    threads:
        4
	params:
		basename="trimmedReads/{sra}.fastq",
		log="trimmedReads/{sra}.log"
	shell:
        """
	trimmomatic PE -threads {threads} {input.read1} {input.read2} -baseout {params.basename} ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
        """

rule trimFastqc:
	"""
	Quality control trough FASTQC
	"""
    input:
        trimreadFQC=rules.trimmomatic.output.forwardPaired
        trimreadRQC=rules.trimmomatic.output.reversePaired
    output:
        zip="rawQC/{sra}_{frr}_fastqc.zip",
	html="rawQC/{sra}_{frr}_fastqc.html"
    threads:
        1
	params:
		path="rawQC/"
	shell:
        """
	fastqc {input.rawread} --threads {threads} -o {params.path}
        """

rule starindex:
	"""
	Alignment with reference genome through STAR
	"""
    input:
        genomefasta="rawReads/Caenorhabditis_elegans.WBcel235.dna.toplevel"
    output:
        genome="starAligned/{sra}GENOME_INDEX",
    threads:
        6
	shell:
        """
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output.genome} --genomeFastaFiles {input.genomefasta}
        """

rule staralign:
	"""
	Alignment with reference genome through STAR
	"""
    input:
        read1=rules.trimmomatic.output.forwardPaired,
	read2=rules.trimmomatic.output.reversePaired
    output:
        bam="starAligned/{sra}Aligned.sortedByCoord.out.bam",
	log="starAligned/{sra}Log.final.out"
    threads:
        45
	params:
		prefix="starAligned/{sra}"
	shell:
        """
	STAR --runThreadN {threads} --genomeDir {rules.starindex.output.genome} --genomeLoad LoadAndKeep --readFilesIn {input.read1} {input.read2} --outFilterIntronMotifs RemoveNoncanonical --outFileName {params.prefix} --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 5000000000 --outReadsUnmapped Fastx
        """

rule picard:
	"""
	Duplicates removal through Picard
	"""
    output:
        bamPC="starAligned/{sra}Aligned.sortedByCoordPICARD.out.bam",
	log="starAligned/{sra}metrics.PICARD.txt"
    shell:
        """
        picard MarkDuplicates --REMOVE_DUPLICATES true -I {rules.staralign.output.bam} -O {output.bamPC} -M {output.metrics}
        """

rules samindex:
    output:
        samF="starAligned/{sra}Aligned.sortedByCoordPICARD.SAMindex.bai"
    shell:
        """
        samtools index {rules.picard.output.bamPC} -o {output.samF}
        samtools flagstats {output.samF}
        samtools coverage {output.samF}
        """
