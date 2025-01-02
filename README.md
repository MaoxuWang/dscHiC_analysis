# Droplet-based high-throughput 3D genome structure mapping of single cells with simultaneous transcriptomics
This repository stores the preprocessing and analysis code used in the **"dscHiC (Honggui Wu, Maoxu Wang et al.)"** paper.

## Pipeline
The subfolder [pipeline] contains:
### 1. Cell barcode Calling
To process dscHi-C data, dscHiCtools <https://github.com/MaoxuWang/dscHiCtools> were developed.
Note that the cell barcode is on the i2 (R1/R2/i2) (or R2 if R1/R2/R3 file format)
- dscHi-C: 1 - 16 bp 
- dscHi-C-multiome: 9 - 24 bp

dscHiCtools firstly identify the true cell barcodes from index reads against whitelist barcodes that allows one mismatch. Barcodes with no or multiple hits will be discarded. Called cell barcodes are added to read name.

### 2. Mapping
Reads then were mapped to reference genome by bwa-mem (v0.7.17) in 5SP mode.

### 3. Contacts extraction
Notably, we did minor modifications to Hickit <https://github.com/MaoxuWang/hickit> to adapt the single-cell features. Briefly, we extracted the cell barcode to each contact information and cell barcodeS were also considered in contacts duplication rather than distance (set parameters --dup-dist=100) alone.

### 4. Imputation
FastHigashi <https://github.com/ma-compbio/Fast-Higashi> and scHiCluster <https://github.com/zhoujt1994/scHiCluster> methods were used in downstream analysis.

## Analysis
All downstream analyses were completed with R language.

The subfolder [analysis] contains:
- cal_hic_feature.R: Function extensively uesd including creation of Hi-C assay in Seurat
- cell_calling.R: identify the true cell barcodes by contacts distribution.
- cell_cycle_phasing.R: assign cell cycle phase to each single cell by contact matrix
- contacts_distribution.R: QC
- marker_AB.heatmap.R: plot function for visualization
- prepareInput.R: preparation datasets for imputation 
- scAB_embedding.R: embedding by scA/B matrix


## Reference
Honggui Wu, Maoxu Wang, Yinghui Zheng, X. Sunney Xie. "Droplet-based high-throughput 3D genome structure mapping of single cells with simultaneous transcriptomics"
