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
The steps 1 to 7 consists of wet lab procedures, while steps 7 to 14 consists of bioinformatics steps, which will be automated in the present project using command line/java/python/R integration. The project is still under development, and the finished steps will be marked down in the following workflow.

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

  8.1. Quality control through *[FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)* and comparison with results before trimming
- [ ] **9. Alignment with reference genome through *[STAR](https://github.com/alexdobin/STAR)***
- [ ] **10. Duplicates removal through *[Picard](https://broadinstitute.github.io/picard/)***
- [ ] **11. Indexing through *[SAMTools](http://www.htslib.org/)***
- [ ] **12. Hits counting through *[FeatureCounts](http://www.bioconductor.org/packages/release/bioc/html/Rsubread.html)***

  12.1. Coverage analysis

  12.2. Normalization for RPKM (reads per kilobase per million of reads mapped) through *[gatk](https://github.com/broadinstitute/gatk/releases)*
- [ ] **13. Quantification through *[BedTools](https://github.com/arq5x/bedtools2)***

  The RAW data is used for quantification, not the normalized through RPKM, which is used for statistical purposes only.
  
- [X] **14. Differential expression analysis through *[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)***

<div align="center">
  <img src="https://i.imgur.com/sKLJN2b.png">
</div>

## Disclaimer
The data used here as example are not mine, please refer to the proper [guidelines](https://holab-hku.github.io/R-workshop/rna-seq-analysis-in-r.html) from Dr Joshua Ho's at the [Bioinformatics and Digital Health Laboratory](https://holab-hku.github.io/).


# PT-BR
## Sobre
Fluxo de trabalho da análise completa de RNAseq com dados brutos no formato .gz.

[Repositório](https://github.com/capuccino26/RNAseq)

## Fluxo de trabalho
As etapas 1 a 7 consistem de técnicas de análises laboratoriais, enquanto as etapas 7 a 14 consistem de etapas de bioinformática, que serão automatizada no presente projeto através da integração de linha de comando/java/python/R. O projeto ainda está em desenvolvimento, de forma que as etapas finalizadas serão marcadas no seguinte fluxo de trabalho.

- [X] **1. Extração de RNA total**
- [X] **2. Isolamento de mRNA**

  2.1. Análise de degradação por *Northern Blotting*
- [X] **3. Fragmentação**
- [X] **4. Conversão em biblioteca de cDNA**
- [X] **5. Adição dos adaptadores**
- [X] **6. Amplificação por PCR**

  6.1. Controle de qualidade dos produtos de PCR quanto à concentração e tamanho dos fragmentos
- [ ] **7. Sequenciamento**

  7.1. Resultados no formato FASTQ compactados como .gz

  7.2. Controle de qualidade através do *[FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)*

  7.3. Análise de contaminação pela comparação com genomas de referência comuns através do BLAST, BWA e kraken2

  7.4. Comparação do conteúdo de GC das sequências analisadas com o esperado do genoma de referência
- [ ] **8. Remoção dos adaptadores e sequências de baixa qualidade através do *[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)***

  8.1. Controle de qualidade através do *[FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)* e comparação com os dados anteriores ao tratamento com *Trimmomatic*
- [ ] **9. Alinhamento com o genoma de referência através do *[STAR](https://github.com/alexdobin/STAR)***
- [ ] **10. Remoção de duplicatas pelo *[Picard](https://broadinstitute.github.io/picard/)***
- [ ] **11. Indexação através do *[SAMTools](http://www.htslib.org/)***
- [ ] **12. Contagem de associações positivas através do *[FeatureCounts](http://www.bioconductor.org/packages/release/bioc/html/Rsubread.html)***

  12.1. Análise de cobertura

  12.2. Normalização por *Reads* por Kilobase por Milhão de *Reads* Mapeados (*reads per kilobase per million of reads mapped - RPKM*) através do *[gatk](https://github.com/broadinstitute/gatk/releases)*
- [ ] **13. Quantificação pelo *[BedTools](https://github.com/arq5x/bedtools2)***

  Os dados brutos são utilizados na quantificação, não os dados normalizados por RPKM, que são utilizados apenas para finalidades estatísticas.
  
- [X] **14. Análise de expressão diferencial através do *[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)***

<div align="center">
  <img src="https://i.imgur.com/Ygr3F67.png">
</div>

## Aviso
Os dados aqui utilizados como exemplo não são de minha autoria, refira-se às [diretrizes adequadas](https://holab-hku.github.io/R-workshop/rna-seq-analysis-in-r.html) do Dr Joshua Ho no [*Bioinformatics and Digital Health Laboratory*](https://holab-hku.github.io/).
