# MODA
Multi-Omics Data Analysis (MODA) is an analitical platform for presicion medicine in Python and R, the MODApy and MODAr family. It allows easy and friendly analysis the Whole Exome Sequencing (WES) data or clinical applications such as SNP identification and priorization and SNP database generation from analyzed patients. In this sence it allows building a home made data base for local and private use.

## Current facilities
The platform can be deployed in a local or remote server acceced trhough a Shiny web interface trhough RStudio Server and Python. The platform allow retrieving raw sequencing data from sequencong providers such as Macrogen(r) and others, by providing approrpiate url to the data. The patients will be downloaded and orgnanized into a local data base.
* SNP hunting: It apply GATK best practices for mutation detection and annotation trhough an optimized pipeline.
* SNP population based database: It allows building a local data base of SNP, with Zigocity information, allele frequency and SNP population frequency. Such type of database is important when analyzing local variant frequencies that may allow filtering or population characterization.
* Single, duos and trios WES analysis and comparison by means of User defined Gene sets implementing state of the art processing pipelines applying good practices protocols. A gene set is a list of genes pathological or pehnotipical associated by the user. It can be used for the study of mendelian/hereditary diseases. The results are saved in an Excel(r) file with clinvar, dbGap, dbSNP annotation plus local allee frequency and population based frequency. The files also hold url links to genecards, snpDB and Varsome web platforms for each detected variant to simplify acces to uptodate information.

## Under development
* Gene Fusion RNAseq-based detection: It is based on [Arriba](https://arriba.readthedocs.io/en/latest/) software,  the winner of the [DREAM SMC-RNA Challenge](https://www.synapse.org/SMC_RNA) in addition with specific visualization capabilities exclusive for MODA family. 
* Gene expression analysis: Gene and exon quantification
* Immune tumor microenvironment deconvolution of RNAseq samples trhough [MIXTURE]


## Current services
* Local and remote Instalation of the MODA family in your local or remote server for institutional use
* We provided WES analysis and advice
* We provide Gene-Fusion, gene expression and immune deconvolution analysis

Please contact for further details efernandez at cidie . ucc . edu . ar

Authors: 

Juan Carlos Vázquez, UTN FRC - UCC

Elmer Andrés Fernández (PhD) - CIDIE - UCC - CONICET 
