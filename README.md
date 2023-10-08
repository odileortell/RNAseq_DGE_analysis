# Differential gene expression analysis of RNA-seq data
This repository contains the code and data for differential gene expression analysis of RNA-seq data.

**README.md**

**Data**

The GC6-74 RNA-seq dataset used in this analysis is available from the NCBI Gene Expression Omnibus (GEO) under accession number GSE94438.

**Methods**

I wrote bash scripts for quality control, trimming, indexing and alignment of reads to the reference genome. The differential gene expression analysis was performed using R code.

1. **Quality control:** FastQC was used to assess the quality of the RNA-seq reads.
2. **Trimming:** Adapters were trimmed using trimmomatic (ver.0.36), java (ver.17.0.2) or HTStream (ver.1.3.1).
3. **Indexing and alignment:** The RNA-seq reads were aligned to the GENCODE human primary genome (release 43) using STAR (ver.2.7.7).
4. **Count matrix generation:** The number of reads mapping to each gene was counted.
5. **Differential gene expression analysis:** The edgeR and limma R package was used to perform differential gene expression analysis between male and female TB patients in R (ver.4.2.1).

**Usage**

To reproduce the results of this analysis, you will need to have the following software installed:

* FastQC
* trimmomatic (ver.0.36)
* java (ver.17.0.2)
* HTStream (ver.1.3.1)
* STAR (ver.2.7.7)
* R (ver.4.2.1)
* edgeR

Once you have installed the required software, you can follow these steps:

1. Download the scripts to your local machine.
2. In the terminal, run the bash scripts in the directory of your data.
3. In R, install the required R packages
4. Run the R script to perform the differential gene expression analysis

**Output**

The output of the differential gene expression analysis includes a table of differentially expressed genes, an MDS plot, a volcano plot, venn diagram, and a heatmap.

