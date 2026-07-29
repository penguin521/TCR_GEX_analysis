Integrated TCR and Single-Cell Gene Expression Analysis

This repository contains an R workflow for integrating single-cell gene-expression (GEX) data with paired T-cell receptor (TCR) sequencing data. The script processes four matched sample sets, summarizes TCR clonotype structure, performs a standard Seurat workflow, matches TCR metadata to single-cell barcodes, and compares clonotype expansion across clusters and preliminary cell-type annotations.

The workflow was developed for the following sample pairs:

Set	Gene-expression sample	TCR sample
Set 1	Pool185_77	Pool185_81
Set 2	Pool185_78	Pool185_82
Set 3	Pool185_79	Pool185_83
Set 4	Pool185_80	Pool185_84
What the Script Does

For each paired GEX/TCR dataset, the script:

Reads Cell Ranger gene-expression and V(D)J output files
Ranks clonotypes by abundance
Assigns clonotypes to frequency categories
Calculates Shannon and Simpson diversity
Summarizes productive TCR chains and paired TRA/TRB recovery
Examines V-gene, J-gene, and CDR3 usage
Identifies expanded and low-frequency clonotypes
Creates and processes a Seurat object
Performs quality control, normalization, PCA, clustering, and UMAP
Identifies cluster marker genes
Matches TCR barcodes to gene-expression cell barcodes
Adds clonotype information to Seurat metadata
Displays TCR expansion categories on the UMAP
Summarizes clonotypes by cluster and preliminary cell type
Saves the final annotated Seurat objects
Generates combined tables and plots across all four sample sets
Compares clonotypes across samples using CDR3 and V/J/CDR3 signatures
Requirements

The script requires R and the following packages:

install.packages(c(
  "tidyverse",
  "janitor",
  "vegan",
  "scales"
))

install.packages("Seurat")

Load the packages with:

library(Seurat)
library(tidyverse)
library(janitor)
library(vegan)
library(scales)

A recent version of R and RStudio is recommended. The full analysis may require substantial memory because four single-cell datasets are processed sequentially.

Expected Input Files

The script expects standard Cell Ranger output folders.

Each GEX folder must contain either:

filtered_feature_bc_matrix.h5

or:

filtered_feature_bc_matrix/

Each TCR folder must contain:

clonotypes.csv
filtered_contig_annotations.csv

A metrics_summary.csv file is optional. When present, the script copies it into the corresponding analysis output folder.

Expected Folder Layout

Place the four GEX output folders under one parent directory:

GEX_outputs/
├── Pool185_77_outs/
├── Pool185_78_outs/
├── Pool185_79_outs/
└── Pool185_80_outs/

Place the four TCR output folders under a second parent directory:

TCR_outputs/
├── Pool185_81_outs/
├── Pool185_82_outs/
├── Pool185_83_outs/
└── Pool185_84_outs/

The output directory can be any separate folder where the generated tables, plots, and Seurat objects should be saved.

Configuration

Download Pool185_TCR_GEX_analysis_portable.R and open it in RStudio or another R editor.

Before running the analysis, edit the three paths in the User configuration section:

gex_root <- "path/to/GEX_outputs"
tcr_root <- "path/to/TCR_outputs"
output_root <- "path/to/analysis_outputs"

Example for Windows:

gex_root <- "C:/Users/YourName/Desktop/GEX_outputs"
tcr_root <- "C:/Users/YourName/Desktop/TCR_outputs"
output_root <- "C:/Users/YourName/Desktop/TCR_GEX_results"

Forward slashes are recommended in R paths, including on Windows.

The script checks the configured paths and stops with a readable error if a required input folder cannot be found.

Running the Analysis

In RStudio:

Download the R script.
Open the script in RStudio.
Update the three directory paths.
Install any missing packages.
Run the complete script by selecting Source.

The script can also be run from a terminal:

Rscript Pool185_TCR_GEX_analysis_portable.R

Progress messages are printed as each sample set is processed.

Main Analysis Settings

The current Seurat workflow uses:

At least 200 detected features when creating each object
Genes detected in at least 3 cells
More than 200 and fewer than 6,000 detected features per retained cell
Less than 15% mitochondrial expression
The first 20 principal components
A clustering resolution of 0.5
A fixed random seed of 185

Clonotypes are grouped into the following categories:

Category	Percentage of TCR barcodes
Dominant	At least 1%
Expanded	0.2% to less than 1%
Low-frequency	0.1% to less than 0.2%
Rare tail	Less than 0.1%

These values are analysis choices rather than universal biological cutoffs. They should be reviewed before applying the script to a different experiment.

Output

The script creates one analysis folder per paired sample:

analysis_outputs/
├── Pool185_81_MASTER_TCR_GEX_analysis/
├── Pool185_82_MASTER_TCR_GEX_analysis/
├── Pool185_83_MASTER_TCR_GEX_analysis/
├── Pool185_84_MASTER_TCR_GEX_analysis/
└── Pool185_81_to_84_MASTER_COMPARISON/

Each sample-level folder contains CSV summaries, PNG plots, and a final annotated Seurat object.

Major sample-level outputs include:

Ranked clonotype tables and frequency bins
TCR diversity summaries
Productive contig and chain summaries
Paired TRA/TRB recovery
V-gene and J-gene usage
CDR3 length distributions
Expanded clonotype CDR3 and V/J details
Seurat QC, PCA, cluster, and marker plots
Barcode-overlap summaries
TCR metadata matched to GEX cells
TCR expansion by cluster and cell type
A final .rds Seurat object for each sample pair

The comparison folder contains combined tables and figures for:

Diversity and clonal burden
Frequency-bin composition
GEX/TCR barcode overlap
TCR expansion by cluster and cell type
Cross-sample CDR3 signatures
Shared V/J/CDR3 clonotypes across sample sets
Preliminary Cell-Type Annotations

The script contains first-pass cluster labels based on marker genes previously reviewed for these four datasets.

These labels are specific to the current clustering results and are not automatically transferable to other datasets. If the analysis is rerun with different filtering, dimensions, clustering parameters, software versions, or samples, cluster numbering may change.

Before using the script with new data:

Review the marker genes for every cluster.
Update the celltype_labels list.
Confirm that each cluster number has an appropriate label.
Treat labels containing terms such as like, unclear, or unassigned as provisional.

The annotations should not be treated as independently validated cell identities.

Adapting the Script to Other Samples

The script is portable across computers, but it is currently organized around four specific Pool185 sample pairs.

To adapt it to another dataset, update:

The sample names and pairings in sample_map
The expected input folder names
The preliminary cluster labels
The marker-gene panel
QC thresholds
PCA dimensions and clustering resolution
Clonotype frequency categories, when appropriate
Any species-specific gene symbols

The barcode-matching function expects standard 10x Genomics barcodes containing a 16-base sequence followed by a numeric suffix, such as:

AAACCCAAGTACGCCC-1
Important Interpretation Notes

A clonotype_id is assigned within a sample and should not be assumed to represent the same receptor in another sample. For cross-sample comparisons, the script constructs receptor signatures from CDR3 sequences and V/J gene assignments.

The script creates both:

A CDR3-based signature using chain and amino-acid sequence
A more restrictive V/J/CDR3 signature using chain, V gene, J gene, and CDR3 sequence

Shared V/J/CDR3 signatures are reported when they occur in more than one sample set.

Barcode overlap depends on matching Cell Ranger barcodes between the GEX and TCR outputs. Low overlap may reflect biological, technical, filtering, or sample-pairing differences and should be investigated rather than interpreted by the script alone.

Limitations
The workflow was written for a specific four-sample analysis.
Preliminary cell-type labels are manually assigned and dataset-specific.
Cluster identities may change when analysis parameters or software versions change.
The script does not perform automated reference-based cell annotation.
The script does not correct batch effects or integrate the four GEX datasets into one joint Seurat object.
Clonotype frequency bins are user-defined analysis categories.
Shared receptor signatures do not by themselves establish shared antigen specificity.
CDR3 and gene assignments depend on the quality of the Cell Ranger V(D)J output.
The workflow should be reviewed and validated before being used for formal biological conclusions.
Reproducibility

The script uses:

set.seed(185)

This improves reproducibility for stochastic Seurat steps. Package versions and session information can also affect results. To document the R environment after running the analysis, use:

sessionInfo()

For stronger reproducibility, users may save the output of sessionInfo() with the analysis results or manage package versions with a tool such as renv.

Intended Use

This script was created as a reusable research workflow for exploring paired TCR repertoire and single-cell gene-expression data. It is intended to reduce repetitive analysis, preserve a consistent workflow across related samples, and generate organized tables and figures for downstream interpretation.

It should be treated as an analysis starting point rather than a fully automated biological interpretation pipeline.

Author

Developed by Steven Xie.

Questions, issue reports, and suggestions for improvement are welcome.
