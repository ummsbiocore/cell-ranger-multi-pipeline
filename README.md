Cell Ranger (v9.0.1) is a set of analysis pipelines that process Chromium single-cell RNA-seq output to align reads, generate feature-barcode matrices and perform clustering and gene expression analysis. 

Steps:
  1. The cellranger bclConvert demultiplexes the Illumina sequencer's base call files (BCLs) for each flow cell directory into FASTQ files. 
  2. The cellranger multi is required to analyze Multiplexed data.
  
Please check following websites for detailed information: 
- `https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/multi`
- BCLConvert v4.3 userGuide: `https://help.dragen.illumina.com/product-guides/dragen-v4.3/bcl-conversion`
- 10X Demultiplexing Guide: `https://www.10xgenomics.com/support/jp/software/cell-ranger-atac/latest/analysis/inputs/direct-demultiplexing-with-illumina-software`
- Illumina BCL Convert SampleSheet Guide: `https://support-docs.illumina.com/SW/BCL_Convert/Content/SW/BCLConvert/SampleSheets_swBCL.htm`

For running BCL Convert for gene expression libraries, here is an example sample sheet CSV that can be used:
```
[Header]
FileFormatVersion,2

[BCLConvert_Settings]
CreateFastqForIndexReads,0

[BCLConvert_Data]
Lane,Sample_ID,index,index2
4,PBMC_60k_Rep1,GCGGGTAAGT,CTTAGTGCTA
```

**Suggested pipelines based on library type: **

| library type | Cell Ranger Pipeline |
| -------- | -------- |
| 3' Gene Expression     | count     |
| 3' Gene Expression + Antibody/CRISPR Guide Capture     | count     |
| Antibody Capture only     | count     |
| 3' Gene Expression + Cell Multiplexing (+ Antibody/CRISPR Guide Capture)     | multi     |
| Flex Gene Expression (+ Antibody/CRISPR Guide Capture)    | multi     |
| 5' V(D)J only     | vdj     |
| 5' Gene Expression only    | count     |
| Antibody Capture (+ Gene Expression)     | count     |
| CRISPR (FB) + Gene Expression     | count     |
| 5' Gene Expression + V(D)J (+ FB)    | multi     |
| 5' Gene Expression + V(D)J + Antigen Capture (BEAM) (+ Antibody Capture)     | multi     |



** Inputs for Cellranger Multi**
- **lanes_10x**: Optional: You can use this column if you want to execute cell ranger multi separately. We will group samples that share the same lanes_10x value and run them separately. Most users do not need to include this column."
- **physical_library_id**: Optional. Library type. Note: by default, the library type is detected automatically based on specified feature_types (recommended). Users typically do not need to include the physical_library_id column in the CSV file.


