<div align="center">
<h1>
  RNA seq
  </h1>
  </div>

# EN-US
## About
Workflow for complete RNAseq analysis from .gz raw data.

[Repository](https://github.com/capuccino26/RNAseq)

## Workflow
The steps 1 to 7 consists of wet lab procedures, while steps 7 to 13 consists of bioinformatics steps, which will be automated in the present project using command line/java/python/R integration. The project is still under development, and the finished steps will be marked down in the following workflow.

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
- [X] **8. Trimming with *[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)***

  8.1. Quality control through *[FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)* and comparison with results before trimming
- [X] **9. Alignment with reference genome through *[STAR](https://github.com/alexdobin/STAR)***
- [X] **10. Duplicates removal through *[Picard](https://broadinstitute.github.io/picard/)***
- [X] **11. Indexing through *[SAMTools](http://www.htslib.org/)***
- [X] **12. Hits counting through *[FeatureCounts](http://www.bioconductor.org/packages/release/bioc/html/Rsubread.html)***

  12.1. Coverage analysis

  12.2. Normalization for RPKM (reads per kilobase per million of reads mapped) through *[gatk](https://github.com/broadinstitute/gatk/releases)* if needed.
  
- [X] **13. Differential expression analysis through *[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)***

<div align="center">
  <img src="https://i.imgur.com/UGuCNlo.png">
</div>

## Disclaimer
The data used here as example are not mine, please refer to the proper [guidelines](https://holab-hku.github.io/R-workshop/rna-seq-analysis-in-r.html) from Dr Joshua Ho's at the [Bioinformatics and Digital Health Laboratory](https://holab-hku.github.io/).


# PT-BR
## Sobre
Fluxo de trabalho da an??lise completa de RNAseq com dados brutos no formato .gz.

[Reposit??rio](https://github.com/capuccino26/RNAseq)

## Fluxo de trabalho
As etapas 1 a 7 consistem de t??cnicas de an??lises laboratoriais, enquanto as etapas 7 a 13 consistem de etapas de bioinform??tica, que ser??o automatizada no presente projeto atrav??s da integra????o de linha de comando/java/python/R. O projeto ainda est?? em desenvolvimento, de forma que as etapas finalizadas ser??o marcadas no seguinte fluxo de trabalho.

- [X] **1. Extra????o de RNA total**
- [X] **2. Isolamento de mRNA**

  2.1. An??lise de degrada????o por *Northern Blotting*
- [X] **3. Fragmenta????o**
- [X] **4. Convers??o em biblioteca de cDNA**
- [X] **5. Adi????o dos adaptadores**
- [X] **6. Amplifica????o por PCR**

  6.1. Controle de qualidade dos produtos de PCR quanto ?? concentra????o e tamanho dos fragmentos
- [ ] **7. Sequenciamento**

  7.1. Resultados no formato FASTQ compactados como .gz

  7.2. Controle de qualidade atrav??s do *[FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)*

  7.3. An??lise de contamina????o pela compara????o com genomas de refer??ncia comuns atrav??s do BLAST, BWA e kraken2

  7.4. Compara????o do conte??do de GC das sequ??ncias analisadas com o esperado do genoma de refer??ncia
- [X] **8. Remo????o dos adaptadores e sequ??ncias de baixa qualidade atrav??s do *[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)***

  8.1. Controle de qualidade atrav??s do *[FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)* e compara????o com os dados anteriores ao tratamento com *Trimmomatic*
- [X] **9. Alinhamento com o genoma de refer??ncia atrav??s do *[STAR](https://github.com/alexdobin/STAR)***
- [X] **10. Remo????o de duplicatas pelo *[Picard](https://broadinstitute.github.io/picard/)***
- [X] **11. Indexa????o atrav??s do *[SAMTools](http://www.htslib.org/)***
- [X] **12. Contagem de associa????es positivas atrav??s do *[FeatureCounts](http://www.bioconductor.org/packages/release/bioc/html/Rsubread.html)***

  12.1. An??lise de cobertura

  12.2. Normaliza????o por *Reads* por Kilobase por Milh??o de *Reads* Mapeados (*reads per kilobase per million of reads mapped - RPKM*) atrav??s do *[gatk](https://github.com/broadinstitute/gatk/releases)* se necess??rio.

- [X] **13. An??lise de express??o diferencial atrav??s do *[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)***

<div align="center">
  <img src="https://i.imgur.com/R5h2CIs.png">
</div>

## Aviso
Os dados aqui utilizados como exemplo n??o s??o de minha autoria, refira-se ??s [diretrizes adequadas](https://holab-hku.github.io/R-workshop/rna-seq-analysis-in-r.html) do Dr Joshua Ho no [*Bioinformatics and Digital Health Laboratory*](https://holab-hku.github.io/).
