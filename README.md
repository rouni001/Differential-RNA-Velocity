# Differentially Veloce Genes (DVG) 

This repository contains code to identify and visualize **Differentially Veloce Genes (DVGs)**—genes exhibiting a statistically significant differential distribution of RNA velocity values when compared to a reference distribution. The pipeline implements statistical testing, multiple hypothesis correction, and interactive cluster visualization using Python.

# Author & Date
Rachid Ounit, Ph.D.

01/17/2025

## Introduction

RNA velocity analysis provides insights into cellular dynamics by estimating the future state of individual cells. Building on this concept, **Differentially Veloce Genes (DVGs)** are defined as genes that display statistically significant differences in their velocity distributions across experimental conditions or regions. This pipeline is tailored to:
- Identify DVGs using differential testing (Welch’s t-test)
- Adjust for multiple comparisons (using the Benjamini-Hochberg procedure)
- Generate interactive clustermaps with dendrograms for in-depth exploration

## Definition of DVG

A **Differentially Veloce Gene (DVG)** is a gene that meets the following criteria:
- Exhibits a differential distribution of RNA velocity values between groups or conditions.
- Achieves statistical significance as determined by a hypothesis test (Welch’s t-test).
- Remains significant after correcting for multiple comparisons with an adjusted p-value threshold.

## Code Overview

The code is organized into functions that:
- Extract and process RNA velocity data.
- Perform differential analysis per region.
- Correct for multiple comparisons.
- Generate interactive visualizations of the resulting DVGs.

### Methods Description

#### `get_df_dv_per_region(name, ds, regions, adata, debug, p_lim)`

- **Purpose:** Iterates over given regions and compiles DVG results from each region into a single DataFrame.
- **Parameters:**
  - `name`: Dataset identifier.
  - `ds`: Index of the condition group.
  - `regions`: List of region names to analyze.
  - `adata`: Annotated data matrix containing RNA velocity layers.
  - `debug`: Boolean flag to enable/disable debug output.
  - `p_lim`: p-value threshold for significance.
- **Returns:** Combined DataFrame with region-specific DVGs.

#### `get_diff_velo_genes(name, ds, region, adata, debug=False, p_lim=1e-9, plot=True, save_plot_path='difference_distribution.png')`

- **Purpose:** Identifies DVGs within a specific region using the Welch's t-test.
- **Details:**
  - Extracts the RNA velocity layer from `adata`.
  - Filters data based on region and condition labels.
  - Compares two pre-defined condition groups.
  - Computes a test statistic and p-value for each gene.
  - Applies a vector normalization (using exponential means) to quantify differences.
  - Adjusts p-values using the Benjamini-Hochberg method.
  - Optionally plots the distribution of the difference metric.
- **Returns:** DataFrame of significant DVGs (genes passing the adjusted p-value cutoff).

#### `generate_interactive_clustermap_with_dendrogram(ds, df, sub_folder_name, value_col='adjusted_p_value', gene_col='gene', region_col='region', output_dir='output', cmap='Viridis', width=1200, height=800)`

- **Purpose:** Creates an interactive clustermap with dendrograms to visualize DVG significance across regions.
- **Procedure:**
  1. Pivots the input DataFrame to form a gene-by-region matrix.
  2. Replaces very small/NaN values and applies a -log10 transformation to the data.
  3. Generates row and column dendrograms using hierarchical clustering.
  4. Reorders the matrix based on the dendrogram ordering.
  5. Constructs a composite figure with interactive subplots (dendrograms and heatmap) using Plotly.
  6. Saves the interactive plot as an HTML file.
- **Returns:** No direct output; the function displays and saves the interactive figure.

#### `get_adata(idx)`

- **Purpose:** Loads annotated data matrices from provided file paths.
- **Details:**
  - Accesses a list of file paths and dataset names.
  - Reads the corresponding `.h5ad` file.
- **Returns:** A tuple containing the dataset name and the AnnData object.

#### `run_analysis(plim_lowest, plim_name)`

- **Purpose:** Automates the analysis across multiple datasets.
- **Procedure:**
  - Iterates through pre-defined datasets.
  - Extracts regions (excluding some that are not of interest).
  - Compiles DVG results using `get_df_dv_per_region`.
  - Generates interactive clustermaps with dendrograms.
- **Parameters:**
  - `plim_lowest`: The p-value threshold used for significance.
  - `plim_name`: Subdirectory name for the output.
- **Execution:** This function is invoked at the end of the script to run the complete pipeline.

## Execution
To run the analysis, simply execute the script. Ensure that:
- The AnnData `.h5ad` files are correctly referenced in the `FILE_ADATA` list.
- The required packages (`numpy`, `pandas`, `plotly`, `scipy`, `statsmodels`, etc.) are installed.

```bash
python <script_name>.py

