**CellRanger Multi** (v9.0.1) pipeline is designed to process Chromium single-cell RNA-seq output to align reads, generate feature-barcode matrices and perform clustering and gene expression analysis. BclConvert demultiplexes the Illumina sequencer's basecall files (BCLs) for each flow cell directory into FASTQ files or the FASTQ files can be given directly.  

## Table of Contents
- [Features](#features)
- [Inputs](#inputs)
- [Outputs](#outputs)
- [Example Datasets](#example-datasets)
- [References and Additional Documentation](#references-and-additional-documentation)

## Features

* **Support for Multiplexing and Immune Profiling**: 	Optional input for cmo_set, feature_reference, and VDJ_reference.  
* **Cloud 	Integration**: Supports reading input data directly from S3 or GCP with automatic handling.  
* **Sample Metadata Parsing**: Integrates sample information via the Metadata input.  
* **Modular Design:** All tools are integrated into the pipeline as modules allowing turning some of the tools on, or off  
* **Scalability**: Adaptable to various compute environments and data sizes.  
* **QC Integration:** Optional FastQC and MultiQC reporting modules (if enabled).
  
Please check following websites for detailed information: 
- `https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/multi`
- BCLConvert v4.3 userGuide: `https://help.dragen.illumina.com/product-guides/dragen-v4.3/bcl-conversion`
- 10X Demultiplexing Guide: `https://www.10xgenomics.com/support/jp/software/cell-ranger-atac/latest/analysis/inputs/direct-demultiplexing-with-illumina-software`
- Illumina BCL Convert SampleSheet Guide: `https://support-docs.illumina.com/SW/BCL_Convert/Content/SW/BCLConvert/SampleSheets_swBCL.htm`

### Key Use cases:

  1. 3' Gene Expression + Cell Multiplexing (+ Antibody/CRISPR Guide 		Capture)  
  2. Flex Gene Expression (+ Antibody/CRISPR Guide Capture)  
  3. 5' Gene Expression + V(D)J (+ FB)  
  4. 5' Gene Expression + V(D)J + Antigen Capture (BEAM) (+ Antibody 		Capture)

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


## **Inputs**
### **Required**

**Run BCL-Convert**

* ***Description**: 	Choose whether to demultiplex BCL files to FASTQ.
* ***Format**: 	DROPDOWN
* ***Options**: 	"yes", "no"
* ***Default**: 	"no"

**Genome Build**

* ***Description**: 	Choose genome reference, e.g., human, mouse, or custom.
* ***Format**: 	DROPDOWN
* ***Options**: 	human_hg38_gencode_v32_cellranger_v6, 	human_hg38_cellranger_GRCh38-2024-A, 	mouse_mm10_gencode_vm23_cellranger_v6, 	mouse_GRCm39_cellranger_GRCm39-2024-A, 	human-mouse_cellranger_GRCh38-GRCm39-2024-A, 	zebrafish_GRCz11plus_ensembl, 	d_melanogaster_BDGP6_32_ensembl_105_cellranger_v6, 	d_melanogaster_flybase_r6_45_cellranger_v6, custom
* ***Default**: 	-

**Run FastQC**

* ***Description**: 	Enable FastQC for raw data quality checking.
* ***Format**: 	DROPDOWN
* ***Options**: 	"yes", "no"
* ***Default**: 	"yes"


**Run Cellranger Multi**

* ***Description**: 	Run and change settings for Cell Ranger Multi.
* ***Format**: 	DROPDOWN
* ***Options**: 	 "yes", "no"
* ***Default**: 	"yes"

**libraries**

* ***Description**: 
	* fastq_id: A name to identify fastq file., 
	* group: Optional. Use this column to run Cell Ranger Multi separately for different groups of files. Fastq files that share the same group value will be grouped together and processed together.
	* feature types: should be chosen from a dropdown list that includes:
  1. *Gene 	Expression*  
  2. *Antibody Capture*  
  3. *CRISPR Guide Capture*  
  4. *Multiplexing Capture*  
  5. *VDJ*  
  6. *VDJ-T*  
  7. *VDJ-T-GD*  
  8. *VDJ-B*  
  9. *Antigen Capture*  
  10. *Custom*	  
* ***Format**: 	TABLE
* ***Example**:

| *fastq_id* | *group* | *feature_types* |
| :---- | :---- | :---- |
| *sample1_fastq1* | *sample1* | *(dropdown)* |
| *sample1_fastq2* | *sample1* | *(dropdown)* |
| *sample2_fastq1* | *sample2* | *(dropdown)* |
| *sample2_fastq2* | *sample2* | *(dropdown)* |

</br>
**Sample Separation (Optional – Only necessary when using multiplexed reactions)** 

* ***Description**: 	
   * **sample_id**: The name to assign to an individual sample from a multiplexed reaction. Must be alphanumeric with hyphens and/or underscores, and less than 64 characters. Required for cell multiplexing libraries.  
   * **group**: Use this column to run Cell Ranger Multi separately for different groups of samples. Samples that share the same group value will be grouped together and processed together. Values in the group column should be consistent with the group column from the "Libraries" table above.
   * **cmo_ids**: Required for 3' cell multiplexing libraries. The cell multiplexing oligo IDs used to multiplex this sample. If multiple CMOs were used per sample, separate IDs with a pipe (e.g., CMO301|CMO302). When using this column for sample separation, you might provide a "cmo_set" file above to override cellranger's default cmo-set. https://cdn.10xgenomics.com/raw/upload/v1689285473/software-support/Cell%20Ranger/download-templates/default_cmo_ref.csv 
   * **hashtag_ids**: Required for 'Cell or sample hashing with Antibody Capture' libraries. The hashtag IDs used to multiplex this sample. If multiple antibody hashtags were used for the same sample, 	you can separate IDs with a pipe (e.g., ABHT-1| ABHT-2) When using this column for sample separation, you must provide a "feature_reference" file above.
   * **ocm_barcode_ids**: Required for 3' and 5' GEM-X On-chip multiplexing (OCM) libraries. The OCM barcode IDs used to multiplex this sample. Must be one of OB1, OB2, OB3, OB4. If multiple OCM Barcodes were used for the same sample, you can separate IDs with a pipe (e.g., OB1|OB2)
   * **probe_barcode_ids**: Required for FLEX libraries. The Fixed RNA Probe Barcode IDs, Antibody Multiplexing Barcode IDs, and CRISPR Multiplexing Barcode IDs used in the experiment.
   * **description**: Optional. A description of the sample.
   * Original documentation link: https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-3p-multi  
* ***Format**: 	TABLE
* ***Example for cmo_ids**:

| *sample_id* | *group* | *cmo_ids* | *hashtag_ids* | *ocm_barcode_ids* | *probe_barcode_ids* | *description* |
| :---: | :---: | :---: | :---: | ----- | :---: | :---: |
| *sample1_1* | *sample1*	 | *cmo_set_id1* |  |  |  |  |
| *sample1_2* | *sample1* | *cmo_set_id2* |  |  |  |  |
| *sample1_3* | *sample1* | *cmo_set_id3* |  |  |  |  |
| *sample2_1* | *sample2* | *cmo_set_id1* |  |  |  |  |

* ***Example for hashtag_ids**:

| *sample_id* | *group* | *cmo_ids* | *hashtag_ids* | *ocm_barcode_ids* | *probe_barcode_ids* | *description* |
| :---: | :---: | :---: | :---: | ----- | :---: | :---: |
| *sample1_1* | *sample1*	 |  | *Feature_reference_id1* |  |  |  |
| *sample1_2* | *sample1* |  | *Feature_reference_id2* |  |  |  |
| *sample1_3* | *sample1* |  | *Feature_reference_id3* |  |  |  |
| *sample2_1* | *sample2* |  | *Feature_reference_id1* |  |  |  |


**Other Cellranger Multi Settings** 

**create_bam**

* ***Description**: 	Setting the option to false reduces the overall computation time and 	decreases the size of the output directory, as it omits the 	generation of the BAM file* 	  
* ***Format**: 	DROPDOWN
* ***Options:** 	TRUE, FALSE
* ***Default**: 	TRUE

**cellranger_multi_chemistry**

* ***Description**: 	Assay configuration. You should only specify chemistry if there is 	an error in automatic detection. Select one of: auto for 	autodetection (default), threeprime for Single Cell 3', fiveprime 	for Single Cell 5', SC5P-PE for paired-end only, SC5P-R2 for 	R2-only, SC3Pv1 or SC3Pv2 or SC3Pv3 for Single Cell 3' v1/v2/v3, 	SC3Pv3HT for Single Cell 3' v3.1 HT or SC-FB for Single Cell 	Antibody-only 3' v2.
* ***Format**: 	DROPDOWN
* ***Options**: 	auto, threeprime, fiveprime, SC5P-PE, SC5P-R2, SC3Pv1, SC3Pv2, 	SC3Pv3, SC3Pv3HT, SC-FB*

**r1_length**

* ***Description**: 	Hard trim the input Read 1 of gene expression libraries to this 	length before analysis.
* ***Format**: 	INPUT
* ***Default**: 	0

**check_library_compatibility**

* ***Description**: 	This option allows users to disable the check that evaluates 10x 	Barcode overlap between libraries when multiple libraries are 	specified (e.g., Gene Expression + Antibody Capture). Setting this 	option to false will disable the check across all library 	combinations. We recommend running this check (default), however if 	the pipeline errors out, users can bypass the check to generate 	outputs for troubleshooting.
* ***Format**: 	DROPDOWN
* ***Options**: 	TRUE, FALSE
* ***Default:** 	TRUE

**include_introns**

* ***Description**: 	Include intronic reads in count.
* ***Format**: 	DROPDOWN
* ***Options**: 	TRUE, FALSE
* ***Default:** 	TRUE

**r1_length_feature**

* ***Description**: 	Hard trim the input Read 1 of Feature Barcode libraries to this length before analysis.
* ***Format**: 	INPUT
* ***Default**: 	0

**Run scRNA-seq Analysis Module

* ***Description**: 	Optional bclconvert parameters
* ***Format**: 	DROPDOWN
* ***Options**: 	"yes", "no"
* ***Default**: 	"yes"

**Run Velocity Analysis Module

* ***Description**: 	Run RNA velocity analysis? This opens options for Velocity module.
* ***Format**: 	DROPDOWN
* ***Options**: 	"yes", "no"
* ***Default**: 	"no"

**Run pySCENIC Module

* ***Description**: 	Run pySCENIC analysis? This opens options for pySCENIC module.
* ***Format**: 	DROPDOWN
* ***Options**: 	"yes", "no"
* ***Default**: 	"no"

**Run Slingshot Module

* ***Description**: 	Run Slingshot analysis? This opens options for Slingshot module.
* ***Format**: 	DROPDOWN
* ***Options**: 	"yes", "no"
* ***Default**: 	"no"

**Download Genomic Sources**

* ***Description**: 	You need to enable Download Genomic Sources option when pipeline needs to prepare custom fasta and gtf files. Here are the cases:
	* When custom genome_source or gtf_source is provided
	* When add_sequences_to_reference option is enabled
	* When custom sequences added to genome.  
* ***Format**: 	DROPDOWN
* ***Options:** 	"yes", "no"
* ***Default:** 	"no"

### **Optional Inputs**

**reads**

* ***Description**: 	Cloud/local dataset for input fastq.gz files. Needs to be given in case of Run BCL-Convert is selected as "no"
* ***Format**: 	FASTQ.GZ

**bclconvert_parameters**

* ***Description**: 	Optional bclconvert parameters
* ***Format**: 	INPUT
* ***Default**: 	 -	

**Demultiplexer_prep**

* ***Description**: 	The SampleSheet.csv. [Sample sheet documentation for BCL Convert](https://support-docs.illumina.com/SW/BCL_Convert/Content/SW/BCLConvert/SampleSheets_swBCL.htm) By default, this file should be located in the BCL directory. If it's missing from that directory, or you need to specify a different location, please provide the full path to the file. For example: s3://bio/BCLfolder/SampleSheet2.csv
* ***Format**: 	INPUT
* ***Default**: 	-
* ***Example**:

| *[Header]* |  |  |
| :---- | :---- | :---- |
| *FileFormatVersion* | *2* |  |
|  |  |  |
| *[BCLConvert_Settings]* |  |  |
| *CreateFastqForIndexReads* | *1* |  |
| *TrimUMI* | *0* |  |
|  |  |  |
| *[BCLConvert_Data]* |  |  |
| *Lane* | *Sample_ID* | *index* |
| *1* | *test_sample_atac* | *TTGTAAGA* |
| *1* | *test_sample_atac* | *GGCGTTTC* |
| *1* | *test_sample_atac* | *CCTACCAT* |
| *1* | *test_sample_atac* | *AAACGGCG* |




**feature_reference**

* ***Description**: 	Feature barcode reference for protein or CRISPR capture. **Required if selecting the 'Antibody Capture' option in the 'feature_type' table**. '**hashing ids**' section of the 'Sample Separation' table should be filled in using the 'id' column of this table. Details: https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-feature-ref-csv
* ***Format**: 	CSV
* ***Example:***

| *id* | *name* | *read* | *pattern* | *sequence* | *feature_type* |
| :---- | :---- | :---- | :---- | :---- | :---- |
| *CD3* | *CD3_UCHT1_TotalC* | *R2* | *^NNNNNNNNNN(BC)NNNNNNNNN* | *CTCATTGTAACTCCT* | *Antibody Capture* |
| *CD19* | *CD19_HIB19_TotalC* | *R2* | *^NNNNNNNNNN(BC)NNNNNNNNN* | *CTGGGCAATTACTCG* | *Antibody Capture* |
| *CD45RA* | *CD45RA_HI100_TotalC* | *R2* | *^NNNNNNNNNN(BC)NNNNNNNNN* | *TCAATCCTTCCGCTT* | *Antibody Capture* |

**cmo_set**

* ***Description**: 	You do not need to provide a Cell Multiplexing Oligos (CMO) reference file unless you want to use a custom one. By default, Cell Ranger will automatically use its built-in CMO set (https://cdn.10xgenomics.com/raw/upload/v1689285473/software-support/Cell%20Ranger/download-templates/default_cmo_ref.csv ) if no cmo_set file is specified. **The 'cmo ids' column in the Sample Separation table should be filled using values from the id column of the Cell Multiplexing Oligos (CMO) reference file.**
* ***Format**: 	CSV
* ***Example:***

***cmo_set:***

| id      | name   | read | pattern | sequence        | feature_type        |
|---------|--------|------|---------|------------------|---------------------|
| CMO301  | CMO301 | R2   | 5P(BC)  | ATGAGGAATTCCTGC | Multiplexing Capture |
| CMO302  | CMO302 | R2   | 5P(BC)  | CATGCCAATAGAGCG | Multiplexing Capture |
| CMO303  | CMO303 | R2   | 5P(BC)  | CCGTCGTCCAAGCAT | Multiplexing Capture |
| CMO304  | CMO304 | R2   | 5P(BC)  | AACGTTAATCACTCA | Multiplexing Capture |
| CMO305  | CMO305 | R2   | 5P(BC)  | CGCGATATGGTCGGA | Multiplexing Capture |
| CMO306  | CMO306 | R2   | 5P(BC)  | AAGATGAGGTCTGTG | Multiplexing Capture |
| CMO307  | CMO307 | R2   | 5P(BC)  | AAGCTCGTTGGAAGA | Multiplexing Capture |
| CMO308  | CMO308 | R2   | 5P(BC)  | CGGATTCCACATCAT | Multiplexing Capture |
| CMO309  | CMO309 | R2   | 5P(BC)  | GTTGATCTATAACAG | Multiplexing Capture |
| CMO310  | CMO310 | R2   | 5P(BC)  | GCAGGAGGTATCAAT | Multiplexing Capture |
| CMO311  | CMO311 | R2   | 5P(BC)  | GAATCGTGATTCTTC | Multiplexing Capture |
| CMO312  | CMO312 | R2   | 5P(BC)  | ACATGGTCAACGCTG | Multiplexing Capture |

**VDJ_reference**

* ***Description**: 	VDJ reference for immune profiling. Required if selecting VDJ-related options in the 'feature_type' of the 'libraries' 	table. A folder or tar.gz archive 
* ***Format**: 	INPUT

**expect_cells**

* ***Description**: 	Expected number of recovered cells.
* ***Format**: 	INPUT
* ***Default: 	-

**force_cells**

* ***Description**: 	Force pipeline to use this number of cells, bypassing cell 	detection.
* ***Format**: 	INPUT
* ***Default: 	-

**r1_length_vdj**

* ***Description**: 	Hard trim the input Read 1 of Vdj libraries to this length before analysis  
* ***Format**: 	INPUT
* ***Default: 	-

**Metadata**

* ***Description**: 	Tab-separated sample metadata file. Information entered here will be available in the metadata layer of the final seurat object.  
* ***Format**: 	TABLE
* ***Required** **Columns**: 	Sample ID, Condition, Replicate  
* ***Example**:

| *SampleID* | *Condition* | *Replicate* |
| :---- | :---- | :---- |
| *Sample01* | *Treated* | *1* |
| *Sample02* | *Control* | *2* |

**Run automated cell-type annotation**

* ***Description**: 	Enable automated cell type annotation. The cells will be annotated using sc-type tool. Link to the paper: https://www.nature.com/articles/s41467-022-28803-w
* ***Format**: 	DROPDOWN
* ***Options**: 	"yes", "no"
* ***Default**: 	"yes"

**Run TCR analysis**

* ***Description**: 	Run TCR analysis? Selecting this option will make TCR analysis shiny app available.
* ***Format**: 	DROPDOWN
* ***Options**: 	"yes", "no"
* ***Default**: 	"no"

**Tissue type**

* ***Description**: 	Primary tissue of origin for the input samples. This is used by the pipeline to select tissue-appropriate reference data and marker sets for automated cell-type annotation.
* ***Format**: 	DROPDOWN
* ***Options**: 	"Immune system", "Pancreas", "Liver", "Eye", "Kidney", "Brain", "Lung", "Adrenal", "Heart"," Intestine", "Muscle"," Placenta", "Spleen", "Stomach", "Thymus", "Hippocampus"
* ***Default:** 	"Immune system"

**Remove Mitochondrial Genes**

* ***Description**: 	When checked, mitochondrial genes will be completely removed prior 	to data filtering.
* ***Format**: 	checkbox
* ***Default:** 	unchecked

**Remove ribosomal RNA Genes**

* ***Description**: 	When checked, ribosomal RNA genes will be completely removed prior 	to data filtering
* ***Format**: 	checkbox
* ***Default:** 	unchecked

**Minimal Genes per Cell**

* ***Description**: 	Threshold quantile for minimum number of genes in a cell. Cells that 	don't meet this threshold will be removed. 
* ***Format**: 	INPUT  
* ***Default:** 	"0.01"

**Maximal Genes per Cell**

* ***Description**: 	Hard trim the input Read 1 of Vdj libraries to this length before	analysis
* ***Format**: 	INPUT
* ***Default:** 	"0.99"

**Minimal UMIs per Cell**

* ***Description**: 	Threshold quantile for minimum number of unique transcript molecules 	in a cell. Cells that don't meet this threshold will be removed.
* ***Format**: 	INPUT
* ***Default:** 	"0.01"

**Maximal UMIs per Cell**

* ***Description**: 	Cutoff quantile for maximum number of unique transcript molecules in 	a cell. Cells that exceed this cutoff will be removed.
* ***Format**: 	INPUT
* ***Default:** 	"0.99"

**Maximal Percent Mitochondrial Reads per Cell**

* ***Description**: 	Cutoff percentage for percentage of reads that come from the 	mitochondria. Cells that exceed this cutoff will be removed.
* ***Format**: 	INPUT
* ***Default:** 	"25"

**Maximal Percent Ribosomal RNA Reads per Cell**

* ***Description**: 	Cutoff percentage for percentage of reads that come from the 	ribosomal RNA. Cells that exceed this cutoff will be removed.
* ***Format**: 	INPUT
* ***Default:** 	"50"

**Remove Doublets**

* ***Description**: 	Whether doublet from sample should be removed. TRUE means that the 	doublet detection/removal will be run and FALSE means it will not. 	DEFAULT means the pipeline will try to detect whether the data is 	from a cellranger multi pipeline, if yes the doublet removal will 	not be run.
* ***Description**: 	When checked, doublets will be removed.
* ***Format**: 	checkbox
* ***Default:** 	unchecked

**Doublet Rate**

* ***Description**: 	Doublet rate per 1K cells to use during doublet removal. The default, 0.008 (e.g. 3.2% doublets among 4000 cells), is appropriate for standard 10X chips. For High Throughput (HT) 10X chips, use half, i.e. 0.004. (Some more recent chips might have this rate even lower)
* ***Format**: 	INPUT
* ***Default:** 	"0.008"

**Doublet Removal Tool**

* ***Description**: 	Tool for doublet removing.
* ***Format**: 	DROPDOWN
* ***Options:** 	"scDblFinder", "DoubletFinder"
* ***Default:** 	"scDblFinder"

**Normalization Method**

* ***Description**: 	Name of normalization method used.
* ***Format**: 	DROPDOWN
* ***Options:** 	"LogNormalize", "CLR", "RC", "SCT"
* ***Default:** 	"LogNormalize"*

**# of Variable Features - Normalization**

* ***Description**: 	Number of variable features to use after ranking by residual 	variance
* ***Format**: 	INPUT
* ***Default:** 	"3000"

**# of Variable Features - PCA_and_Batch_Effect_Correction**

* ***Description**: 	Use this many features as variable features after ranking by 	residual variance
* ***Format**: 	INPUT
* ***Default:** 	"3000"

**Selection Method**

* ***Description**: 	Doublet percentage to use during doublet removal.
* ***Format**: 	DROPDOWN
* ***Options:** 	"vst", "mean.var.plot", "dispersion"
* ***Default:** 	"vst"

**Correct Batch Effect**

* ***Description**: 	Choose whether to do batch effect correction.*  
* ***Format**: 	DROPDOWN
* ***Options**: 	"TRUE", "FALSE"
* ***Default:** 	"TRUE"

**Weighted Nearest Network assay**

* ***Description**: 	If the data is multi-modal/omics, it is possible to leverage more 	than one assay in the downstream analysis.
* ***Format**: 	INPUT
* ***Default:** 	-

**Minimum Resolution**

* ***Description**: 	Minimum resolution for clustering parameter selection.
* ***Format**: 	INPUT
* ***Default:** 	"0.1"

**Maximum Resolution**

* ***Description**: 	Maximum resolution for clustering parameter selection.
* ***Format**: 	INPUT
* ***Default:** 	"2.0"

**# of Principal Components**

* ***Description**: 	Number of principal components to build UMAP, tSNE and nearest 	neighbor graph. Enter 0 for automated prediction.
* ***Format**: 	INPUT
* ***Default:** 	"0"

**Find Markers for All Resolutions**

* ***Description**: 	Whether to find cluster markers for all resolutions. This can 	significantly increase computational time and cost.
* ***Format**: 	CHECKBOX
* ***Default:** 	unchecked

**Generate_loom_file**

* ***Description**: 	Whether to generate a loom file from the final analysis results.
* ***Format**: 	DROPDOWN
* ***Options**: 	"TRUE", "FALSE"
* ***Default:** 	"FALSE"

**replace_geneID_with_geneName**

* ***Description**: 	Hard trim the input Read 1 of Vdj libraries to this length before 	analysis
* ***Format**: 	DROPDOWN
* ***Options:** 	"yes", "no"
* ***Default:** 	"no"

**mask_gtf**

* ***Description**: 	.gtf file containing intervals to mask. It might be desirable to mask expressed repetitive elements, since those count could constitute a confounding factor in the downstream analysis.
* ***Format**: 	gtf

**Condition Column**

* ***Description**: 	Name of condition column specified in Metadata table
* ***Format**: 	INPUT
* ***Default:** 	"Condition"

**k_steps**

* ***Description**: 	K steps to compute transition probabilities in velocity analysis. Higher K = longer-term predictions. Only the step numbers written will be used.
* ***Format**: 	INPUT
* ***Default:** 	"1,2,3"

**Number of Macrostates**

* ***Description**: 	Number of course-grained cell states for fate probability analysis. High values = more detailed trajectory structure.
* ***Format**: 	INPUT
* ***Default:** 	"10"

**Target Clusters**

* ***Description**: 	Select specific clusters to include in analysis. Leave empty to analyse all clusters.
* ***Format**: 	INPUT
* ***Default:** 	""

**Method to use in GRN inference**

* ***Description**: 	Algorithm to use for gene regulatory network inference. GRNBoost2 is the original method used in pySCENIC workflow but it is the slowest one. Using `regdiffusiton_GPU` will speed up the analysis significantly but may cost more.
* ***Format**: 	DROPDOWN
* ***Options:** 	"GRNBoost2", "regdiffusion_GPU", "regdiffusion_CPU"
* ***Default:** 	"GRNBoost2"

**Number of Highly Variable Genes**

* ***Description**: 	Number of Highly Variable Genes to be included. This limits the analysis to use only a set of the most highly variable genes. This is required if `regdiffusion_GPU` GRN method is selected. A dataset consisting of 38K cells with top 15K highly variable genes requires 20GB of GPU ram, whereas using 3K highly variable genes requires 1.5GB of GPU ram. Please scale these parameters accordingly to avoid high costs.
* ***Format**: 	INPUT
* ***Default:** 	""

**Mask Dropouts**

* ***Description**: 	Using this option excludes zero-expression cells when calculating TF-target correlations, focusing only on cells where both genes are expressed. This affects how the correlation between TF and target genes is calculated
* ***Format**: 	CHECKBOX
* ***Default:** 	unchecked

**AUC Threshold**

* ***Description**: 	Fraction of the ranked gene list in each cell that is used when computing the AUC (Area Under the Curve) score for each regulon.
* ***Format**: 	INPUT
* ***Default:** 	"0.05"

**Reduction**

* ***Description**: 	Reduction to use in slingshot analysis
* ***Format**: 	DROPDOWN
* ***Options:** 	"umap", "tsne"
* ***Default:** 	"umap"

**Number of Genes**

* ***Description**: 	Number of top variable genes to use in fitgam analysis. If empty, all of top variable genes will be used. Set this value for big datasets to make the analysis faster. The lower the value, faster the analysis will be but the results will be less accurate.
* ***Format**: 	INPUT
* ***Default:** 	""

### **Outputs**

**geneexpressionmatrix.csv:**

* ***Description**: 	Normalized gene expression matrix
* ***Format**: 	 CSV
* ***Location:** 	results/geneexpressionmatrix.csv*

**summaryreport.html:**

* ***Description**: 	 Summary report of run
* ***Format**: 	 HTML
* ***Location:** 	results/report.html

**alignedreads.bam:**

* ***Description**: 	Aligned reads in BAM format
* ***Format**: 	BAM
* ***Location**: 	results/alignments/alignedreads.bam

**fastqcreport.html:**

* ***Description**: 	FastQC report if enabled
* ***Format**: 	HTML
* ***Location**: 	results/qc/fastqcreport.html

**scvelo_out.h5ad:**

* ***Description**: 	H5AD file including velocity analysis results. scVelo shiny app is used to explore the results.
* ***Format**: 	h5ad
* ***Location**: 	scVelo_out/scvelo_out.h5ad

**sce_fitgam.rds:**

* ***Description**: 	rds file including Seurat object with slingshot analysis results. Slingshot shiny app is used to explore the results.
* ***Format**: 	rds
* ***Location**: 	fitgam_rds/sce_fitgam.rds

**scenic_integrated.loom:**

* ***Description**: 	Updated loom file with the AUCell scores added
* ***Format**: 	loom
* ***Location**: 	pySCENIC_loom/scenic_integrated.loom

**pyscenic_out.zip:**

* ***Description**: 	Zip file including pySCENIC pipeline result files: adjacencies.csv, regulons.csv and aucell_matrix.csv.
* ***Format**: 	zip
* ***Location**: 	pySCENIC_out/pyscenic_out.zip


### Example Datasets:

#### **Example 1. PBMCs of a Healthy Donor (v1) - 5' Gene Expression + Cell Surface Protein Libraries (VDJ-B + VDJ-T)**
- **Source**: https://www.10xgenomics.com/datasets/pbm-cs-of-a-healthy-donor-v-1-1-1-standard-3-1-0
- **Dataset**: https://www.viafoundry.com/test_data/cellranger_multi/fastq_PBMC-VDJ-GEX-downsampled
- **feature_reference**:
https://www.viafoundry.com/test_data/cellranger_multi/fastq_PBMC-VDJ-GEX-downsampled/vdj_v1_hs_pbmc3_feature_ref.csv
| id                    | name                                  | read | pattern                  | sequence        | feature_type     |
| --------------------- | ------------------------------------- | ---- | ------------------------ | --------------- | ---------------- |
| CD3                   | CD3_UCHT1_TotalC                      | R2   | ^NNNNNNNNNN(BC)NNNNNNNNN | CTCATTGTAACTCCT | Antibody Capture |
| CD19                  | CD19_HIB19_TotalC                     | R2   | ^NNNNNNNNNN(BC)NNNNNNNNN | CTGGGCAATTACTCG | Antibody Capture |
| CD45RA                | CD45RA_HI100_TotalC                   | R2   | ^NNNNNNNNNN(BC)NNNNNNNNN | TCAATCCTTCCGCTT | Antibody Capture |
| CD4                   | CD4_RPA-T4_TotalC                     | R2   | ^NNNNNNNNNN(BC)NNNNNNNNN | TGTTCCCGCTCAACT | Antibody Capture |
| CD8a                  | CD8a_RPA-T8_TotalC                    | R2   | ^NNNNNNNNNN(BC)NNNNNNNNN | GCTGCGCTTTCCATT | Antibody Capture |
| CD14                  | CD14_M5E2_TotalC                      | R2   | ^NNNNNNNNNN(BC)NNNNNNNNN | TCTCAGACCTCCGTA | Antibody Capture |
| CD16                  | CD16_3G8_TotalC                       | R2   | ^NNNNNNNNNN(BC)NNNNNNNNN | AAGTTCACTCTTTGC | Antibody Capture |
| CD56                  | CD56_QA17A16_TotalC                   | R2   | ^NNNNNNNNNN(BC)NNNNNNNNN | TTCGCCGCATTGAGT | Antibody Capture |
| CD25                  | CD25_BC96_TotalC                      | R2   | ^NNNNNNNNNN(BC)NNNNNNNNN | TTTGTCCTGTACGCC | Antibody Capture |
| CD45RO                | CD45RO_UCHL1_TotalC                   | R2   | ^NNNNNNNNNN(BC)NNNNNNNNN | CTCCGAATCATGTTG | Antibody Capture |
| PD-1                  | PD-1_EH12.2H7_TotalC                  | R2   | ^NNNNNNNNNN(BC)NNNNNNNNN | ACAGCGCCGTATTTA | Antibody Capture |
| TIGIT                 | TIGIT_A15153G_TotalC                  | R2   | ^NNNNNNNNNN(BC)NNNNNNNNN | TTGCTTACCGCCAGA | Antibody Capture |
| isotype_control_IgG1  | isotype_control_IgG1_MOPC-21_TotalC   | R2   | ^NNNNNNNNNN(BC)NNNNNNNNN | GCCGGACGACATTAA | Antibody Capture |
| isotype_control_IgG2a | isotype_control_IgG2a_MOPC-173_TotalC | R2   | ^NNNNNNNNNN(BC)NNNNNNNNN | CTCCTACCTAAACTG | Antibody Capture |
| isotype_control_IgG2b | isotype_control_IgG2b_MPC-11_TotalC   | R2   | ^NNNNNNNNNN(BC)NNNNNNNNN | ATATGTATCACGCGA | Antibody Capture |
| CD127                 | CD127_A019D5_TotalC                   | R2   | ^NNNNNNNNNN(BC)NNNNNNNNN | GTGTGTTGTCCTATG | Antibody Capture |
| CD15                  | CD15_W6D3_TotalC                      | R2   | ^NNNNNNNNNN(BC)NNNNNNNNN | TCACCAGTACCTAGT | Antibody Capture |


- **Libraries Section**
| *fastq_id* | *group* | *feature_types* |
| -------- | -------- | -------- |
|5gex_protein_antibody	| |		Antibody Capture| 
|5gex_protein_gex		| |	Gene Expression| 
|vdj-b		| |	VDJ-B| 
|vdj-t	| |		VDJ-T| 

- **Sample Separation** 
Not required

#### **Example 2. PBMCs - Multiplexed - CMOs**

* **Source**: 10x Datasets
* **Dataset**: https://www.viafoundry.com/test_data/cellranger_multi/fastq_PBMCs_Multiplexed_2CMOs_subsampled
- **feature_reference**:
https://www.viafoundry.com/test_data/cellranger_multi/fastq_PBMCs_Multiplexed_2CMOs_subsampled/SC3_v3_NextGem_DI_CellPlex_Mouse_PBMC_10K_Multiplex_count_feature_reference.csv
| id     | name   | read | pattern | sequence        | feature_type         |
| ------ | ------ | ---- | ------- | --------------- | -------------------- |
| CMO301 | CMO301 | R2   | 5P(BC)  | ATGAGGAATTCCTGC | Multiplexing Capture |
| CMO302 | CMO302 | R2   | 5P(BC)  | CATGCCAATAGAGCG | Multiplexing Capture |
| CMO303 | CMO303 | R2   | 5P(BC)  | CCGTCGTCCAAGCAT | Multiplexing Capture |
| CMO304 | CMO304 | R2   | 5P(BC)  | AACGTTAATCACTCA | Multiplexing Capture |
| CMO305 | CMO305 | R2   | 5P(BC)  | CGCGATATGGTCGGA | Multiplexing Capture |
| CMO306 | CMO306 | R2   | 5P(BC)  | AAGATGAGGTCTGTG | Multiplexing Capture |
| CMO307 | CMO307 | R2   | 5P(BC)  | AAGCTCGTTGGAAGA | Multiplexing Capture |
| CMO308 | CMO308 | R2   | 5P(BC)  | CGGATTCCACATCAT | Multiplexing Capture |
| CMO309 | CMO309 | R2   | 5P(BC)  | GTTGATCTATAACAG | Multiplexing Capture |
| CMO310 | CMO310 | R2   | 5P(BC)  | GCAGGAGGTATCAAT | Multiplexing Capture |
| CMO311 | CMO311 | R2   | 5P(BC)  | GAATCGTGATTCTTC | Multiplexing Capture |
| CMO312 | CMO312 | R2   | 5P(BC)  | ACATGGTCAACGCTG | Multiplexing Capture |

- **Libraries section**
| *fastq_id* | *group* | *feature_types* |
| -------- | -------- | -------- |
| SC3_v3_NextGem_DI_CellPlex_Mouse_PBMC_10K_1_multiplexing_capture	| group1	| Multiplexing Capture| 
| SC3_v3_NextGem_DI_CellPlex_Mouse_PBMC_10K_2_multiplexing_capture| 	group1 | 	Multiplexing Capture| 
| SC3_v3_NextGem_DI_CellPlex_Mouse_PBMC_gex_sub	| group1	| Gene Expression| 


- **Sample Separation** 
| *sample_id* | *cmo_ids* | *description* |
| -------- | -------- | -------- |
| PBMCs_mouse_1	 |	CMO309	 |			PBMCs_mouse_1 |
| PBMCs_mouse_2	 |	CMO310		 |		PBMCs_mouse_2 |


#### **Example 3. Demultiplexing with BCLConvert:**
- Source: https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-direct-demultiplexing-bcl-convert
- **Dataset**:https://cf.10xgenomics.com/supp/spatial-exp/demultiplexing/iseq-DI.tar.gz
- **Sample Sheet:** https://cf.10xgenomics.com/supp/spatial-exp/demultiplexing/bcl_convert_samplesheet.csv

```
[Header]			
FileFormatVersion	2		
			
[BCLConvert_Settings]			
CreateFastqForIndexReads	0		
			
[BCLConvert_Data]			
Lane	Sample_ID	index	index2
1	iseq-DI	GTAACATGCG	AGGTAACACT
```

- **Introduction:**
	- First download BCL directory and untar it.(tar -xf /working-directory/iseq-DI.tar.gz)
	- Upload BCL directory and CSV file to location where Foundry can reach. (Cloud buckets or HPC clusters)
	- Enable **Run BCL-Convert** input of the pipeline
	- Under Demultiplexer_prep section, click **Add** button and enter necessary locations:
		- **bcl_directory**: enter location of bcl_directory (not tar or gzipped)
		- **SampleSheet**: Sample sheet documentation for BCL Convert. By default, SampleSheet.csv file should be located in the BCL directory. If it’s missing from that directory or you need to specify a different location, please provide the full path to the file. For example: s3://bio/BCLfolder/SampleSheet2.csv
- **Notes:**
	- To test BCL Convert you can use the iseq-DI example dataset. This dual-indexed iSeq dataset has been selected for its small size (541 MB). It should not be used to **run downstream pipelines**.


### References and Additional Documentation

**Additional CellRanger Multi information:**  
[https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/bcl2fastq-direct](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/bcl2fastq-direct)  
[https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/multi](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/multi)  
 **Related Papers/links:**  
 - https://www.10xgenomics.com/support/software/cell-ranger/latest  
 - https://www.10xgenomics.com/support/software/cell-ranger/downloads  
 **Pipeline Repository**: https://github.com/10XGenomics/cellranger  

 ***License**: https://github.com/10XGenomics/cellranger?tab=License-1-ov-file

