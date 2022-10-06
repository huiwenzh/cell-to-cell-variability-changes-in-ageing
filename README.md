# Cell-to-cell variability changes in complex biological processes like differentiationa and ageing


## This repository include the scripts to generate for paper: 
We used the **Tabula Muris Senis** dataset to perform a systematic evaluation of 14 cell-to-cell variability metrics that are generic or transcriptomic-data specific and demonstrated the significant impact of cell-to-cell variability changes during the B lymphocytes differentiaion processes and ageing. 

### Data accesses
**Tabula muris sensis Marrow datasets**

The raw data was downloaded from [link](https://figshare.com/articles/dataset/Tabula_Muris_Senis_Data_Objects/12654728) for both FACS-sorted SMARTSeq2 sequencing and 10X droplet-based sequencing data. In brief, we selected cell types with different samples sizes in two sequencing platforms and similar samples within one sequencing platform. More preprocessing steps please refer to our [paper].

**Large B dataset**

One external dataset was used to quantify the impact of immense sample sizes, which can be found [here](https://www.10xgenomics.com/resources/datasets/cd-19-plus-b-cells-1-standard-1-1-0).

**Techinical spike-ins dataset**

Spike-in dataset was used as a negative control in the evaluation of metrics performance, which can be found [GSE54695](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54695). In theory, there is no genuine biological variability derived from this experiment, so this dataset represents a limited level of variability compared to other gene expression data.

**Simulation datasets**

We simulated *four* datasets with different cell-cell variation at gene levels and different degrees of dropout rate from [splatter](https://bioconductor.org/packages/release/bioc/html/splatter.html). We used the parameters from HSC droplet datasets with 400 known highly variable genes. 

### Normalisation 
**Generic metrics** such as *SD*, *IQR*, *MAD*, *CV*, *FF* requires normalised values as input so that CPM-normalisation was applied to keep consistency. 
**Transcriptomic-data-specific metrics** were processed with the in-built normalisation methods so raw data was used as input.

