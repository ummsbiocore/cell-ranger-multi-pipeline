Cell Ranger (v7.0.0) is a set of analysis pipelines that process Chromium single-cell RNA-seq output to align reads, generate feature-barcode matrices and perform clustering and gene expression analysis. 

Steps:
  1. The cellranger mkfastq demultiplexes the Illumina sequencer's base call files (BCLs) for each flow cell directory into FASTQ files. 
  2. The cellranger multi is required to analyze 3' Cell Multiplexing data. Otherwise, users can continue to use cellranger count.
  2. Cellranger count takes FASTQ files performs alignment, filtering, barcode counting, and UMI counting. It uses the Chromium cellular barcodes to generate feature-barcode matrices, determine clusters, and perform gene expression analysis.  cellranger count also processes Feature Barcoding data alongside Gene Expression reads.
  
Please check following web sites for detailed information: 
- `https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count`  
- `https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/multi`


** Inputs for Cellranger Multi**
- **lanes_10x**: Optional: You can use this column if you want to execute cell ranger multi separately. We will group samples that share the same lanes_10x value and run them separately. Most users do not need to include this column."
- **physical_library_id**: Optional. Library type. Note: by default, the library type is detected automatically based on specified feature_types (recommended). Users typically do not need to include the physical_library_id column in the CSV file.


**Inputs for Cellranger Count**

1. sample: Sample name as specified in the sample sheet supplied to cellranger mkfastq. Allowable characters in sample names are letters, numbers, hyphens, and underscores.
2. lane of library: Optional lane id.
3. library types: Select one of these options:"Gene Expression","Antibody Capture","CRISPR Guide Capture","Multiplexing Capture","VDJ-T"


**Cellranger Documentation:**
- <a class="link-underline" target="_blank" href="https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/bcl2fastq-direct"> Demultiplexing FASTQs with bcl2fastq</a>
- <a class="link-underline" target="_blank" href="https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count"> Single-Library Analysis with cellranger count</a>
- <a class="link-underline" target="_blank" href="https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/multi"> Cell Multiplexing with cellranger multi</a>

