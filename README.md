<div align="center">
<h1>
  RNA seq
  </h1>
  </div>

# EN-US
## About
Workflow for complete RNAseq analysis from .gz raw data.

[Repository](https://github.com/capuccino26/RNAseq)

# Workflow
The steps 1 to 7 consists of wet lab procedures, while steps 7 to 14 consists of bioinformatics steps, which will be automated in the present project using python/R integration. The project is still under development, and the finished steps will be marked down in the following workflow.

- [X] **1. Total RNA Extraction**
- [X] **2. mRNA isolation**

  2.1. Northern Blotting for degradation control
- [X] **3. Fragmentation**
- [X] **4. cDNA library convertion**
- [X] **5. Addition of Adapters**
- [X] **6. Amplification through PCR**

  6.1. Quality control of products concentration and length
- [ ] **7. Sequencing**

  7.1. Results in the FASTQ files compacted as .gz

  7.2. Quality control trough *[FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)*

  7.3. Contamination analysis through comparison of common reference genomes through BLAST, BWA and kraken2

  7.4. GC content analysis compared to the expected in the reference genome
- [ ] **8. Trimming with *[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)***

  8.1. Quality control through *[FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)* and comparsion with results before trimming
- [ ] **9. Alignment with reference genome through *[STAR](https://github.com/alexdobin/STAR)***
- [ ] **10. Duplicates removal through *[Picard](https://broadinstitute.github.io/picard/)***
- [ ] **11. Indexing through *[SAMTools](http://www.htslib.org/)***
- [ ] **12. Hits counting through *[FeatureCounts](http://www.bioconductor.org/packages/release/bioc/html/Rsubread.html)***

  12.1. Coverage analysis

  12.2. Normalization for RPKM (reads per kilobase per million of reads mapped) through *[gatk](https://github.com/broadinstitute/gatk/releases)*
- [ ] **13. Quantification through *[BedTools](https://github.com/arq5x/bedtools2)***
- [X] **14. Differential expression analysis through *[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)***
