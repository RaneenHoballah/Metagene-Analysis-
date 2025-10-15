# Metagene-Analysis-
# Metagene analysis pipeline for ChIP-seq data developed during my MSc thesis in Bioinformatics.
# 🧬 Metagene Analysis of ChIP-seq Data for MSCC Thesis in Bioinformatics

This repository contains the workflows and scripts developed for my MSc thesis in Bioinformatics, under the supervision of Dr. Saverio Brogna at the University of Birmingham.  
The project focuses on performing metagene analysis using ChIP-seq data to explore differences in transcriptional regulation between **protein-coding** and **noncoding (lncRNA)** genes in *Schizosaccharomyces pombe*.

---

## 📖 Background

Gene expression in eukaryotic organisms involves both **protein-coding genes** and **noncoding RNA (ncRNA)** genes.  
While protein-coding genes produce stable proteins required for cellular function, noncoding RNAs—including **long noncoding RNAs (lncRNAs)**—do not encode proteins but share key molecular features with mRNAs, such as 5′ capping, 3′ polyadenylation, and splicing.

Despite these similarities, the **role of lncRNAs** in gene expression and regulation remains poorly understood.  
Understanding their transcriptional behavior can shed light on key aspects of **gene regulation, development, and cellular differentiation** in eukaryotes.

---

## 🔬 Research Objective

This project investigates whether the transcription of **protein-coding genes** differs from that of **lncRNAs** by analyzing **RNA Polymerase II (RNAPII)** ChIP-seq data.

### 🧠 Hypothesis
> The transcription of protein-coding genes differs from the transcription of long noncoding RNAs (lncRNAs).

---

## ⚙️ Methods Overview

ChIP-seq datasets specific to **RNA Polymerase II (RNAPII)** were analyzed to examine transcriptional patterns across gene bodies.  
The analysis pipeline includes:

1. **Data Acquisition**
   - Downloading publicly available ChIP-seq datasets from [Larochelle et al., 2018](https://doi.org/10.1038/s41467-018-06758-z) and data provided by Dr. Brogna’s laboratory.

2. **Preprocessing**
   - Quality control and filtering using **Trimmomatic**
   - Sequence alignment with **Bowtie**
   - File handling and conversion using **SAMtools** 

3. **Data Conversion**
   - Converting alignment files (`.bam`) into Python-compatible digital formats for analysis.
   - Generating **pickle (`.pickle`) files** to store processed data efficiently.  
     These pickle files are then used as input to **start the metagene analysis**.

4. **Metagene Analysis**
   - Generating average transcription profiles (metagenes) across gene bodies for coding and noncoding genes using **Python**.
   - Comparing profiles between gene categories to identify transcriptional differences.

---

## 🧰 Tools and Libraries

- **Computational Environment:** BlueBEAR HPC
- **Languages:** Python, Bash, and Linux 
- **Key Tools:**  
  `SRA Toolkit`, `SAMtools`, `Trimmomatic`, `Bowtie`
- **Python Libraries:**  
  `pandas`, `numpy`, `matplotlib`


