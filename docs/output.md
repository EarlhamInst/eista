# nf-core/eista: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Primary analysis](#primary-analysis)
  - [Vizgen post-processing analysis](#vizgen-post-processing-analysis) - Cell segmentation and quantification
- [Secondary analysis](#secondary-analysis)
  - [QC & cell filtering](#qc--cell-filtering) - Cell filtering and QC on raw data and filtered data
  - [Clustering analysis](#clustering-analysis) - Single-cell clustering analysis
  - [Spatial statistics analysis](#spatial-statistics-analysis) - Single-cell spatial statistics analysis
- Tertiary analysis
  - [Cell-type annotation analysis](#annotation-analysis) - Single-cell cell-type annotation analysis
  - [Cell-type annotation analysis with scvi-tools](#annotation-analysis) - Single-cell cell-type annotation analysis with scvi-tools
  - [Differential expression analysis](#dea-analysis) - Single-cell differential expression analysis
  - [Cell-cell communication analysis](#cellchat-analysis) - Single-cell cell-cell communication analysis
- [Pipeline reporting](#pipeline-reporting)
  - [Analysis report](#analysis-report) - Single-ell Analysis Report
  - [MultiQC](#multiqc) - Aggregate report describing results and QC for tools registered in nf-core
  - [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution


## Primary analysis

### <u>Vizgen post-processing analysis</u>

**Output directory: `results/vizgen`**
- `<sample_name>/`: Contains output files from Vizgen post-processing too VPT.
  - `counts_unfiltered/`
    - `cell_by_gene.csv`: The raw count of transcripts of each targeted gene in each cell.
    - `cellpose_micron_space.parquet`: The primary output of the segmentation. This table contains the EntityIDs for the cells and the geometries in units of microns.
    - `detected_transcripts.csv`: This file has the information about the detected genes and their locations.
    - `cell_metadata.csv`: The cell metadata file has annotation about the location, size, and shape of each cell that can be used to identify cell neighbors, sort cells into cell types, filter low quality cells, etc.
    - `sum_signals.csv`: The sum signals file has information about the brightness of each mosaic tiff image within each cell.
    - `updated_*.vzg`: An updated vzg file with new entity geometries and a new Cell-by-Entity matrix. The file is used for visualization in the MERSCOPE Desktop Vizualizer software.
- `mtx_conversions/`
  - `<sample_name>/`
    - `sample_name_st_matrix.h5ad`: AnnData object file for this sample.
  - `combined_st_matrix.h5ad`: AnnData object file for combined samples.


## Secondary analysis

### <u>QC & cell filtering</u>

**Output directory: `results/qc_cell_filter`**
- `sample_summary.csv`: overall summary of the single-cell count matrix
- `adata_filtered_normalized.h5ad`: AnnData object file after cell filtering and normalization
- `raw_counts/sample_*/`
  - `scatter_counts_genes_volume.png`: scatter plot shows the relationship between total read counts and the number of genes.
  - `knee_plot.png`: knee plot shows the relationship between cell segmentation rank and total transcript count per cell segmentation.
  - `violin*.png`: violin plots display the distribution of cells based on the number of genes, total counts, and cell volumes.
  - `histograms.png`: histogram plots display the distribution of cells based on the number of genes, total counts, and cell volumes. 
- `cell_filtering/`
  - `highly_variable_genes.png`: plot of mean expressions against dispersions of genes for highly variable genes.
  - `sample_summary_filtered.csv`: overall summary of the single-cell count matrix after cell filtering
  - `sample_*/`
    - `spatial_scatter_total_counts_genes.png`: spatial scatter plots for the number of genes, total counts.
    - `umap_total_counts_genes.png`: UMAP plots for the number of genes, total counts.
    - `violin*.png`: violin plots display the distribution of cells based on the number of genes, total counts, and cell volumes.
    - `histograms.png`: histogram plots display the distribution of cells based on the number of genes, total counts, and cell volumes.
- `parameters.json`: a JSON file containing the parameter settings in the analysis.
    

### <u>Clustering analysis</u>

**Output directory: `results/clustering`**
- `adata_clustering.h5ad`: AnnData object file after clustering analysis.
- `sample_*/` or `group_*/`
  - `umap_leiden_res_*.png`: UMAP plots showing clustering results with differnt resoultuion settings.
  - `spatial_scatter_leiden_res_*.png`: spatial scatter plots for the corresponding resolutions of clustering.
- `resolution_*/`
  - `prop_leiden_res_*.png`: plot showing a stacked bar chart that presents the proportions of clusters across samples/groups.
- `parameters.json`: a JSON file containing the parameter settings in the analysis.


### <u>Spatial statistics analysis</u>

**Output directory: `results/spatialstats`**
- `sample_*/` or `group_*/`
  - `centrality_scores_leiden_res_*.png`: a plot showing 3 centrality scores include: Closeness centrality, Degree centrality, and Clustering coefficient.
  - `spatial_scatter_closeness_high_*.png`: spatial scatter plots highlight clusters/groups with high closeness centrality
  - `spatial_scatter_closeness_low_*.png`: spatial scatter plots highlight clusters/groups with low closeness centrality
  - `spatial_scatter_degree_high_*.png`: spatial scatter plots highlight clusters/groups with high degree centrality
  - `spatial_scatter_degree_low_*.png`: spatial scatter plots highlight clusters/groups with low degree centrality
  - `spatial_scatter_clustering_high_*.png`: spatial scatter plots highlight clusters/groups with high clustering coefficient
  - `spatial_scatter_clustering_low_*.png`: spatial scatter plots highlight clusters/groups with low clustering coefficient
  - `neighbors_enrichment_leiden_res_*.png`: the heatmap shows the results of the neighborhood enrichment analysis on clusters. 
  - `autocorr_moranI.csv`: a csv file shows Moran’s I statistics for autocorrelation analysis of the top 50 and bottom 50 genes.
  - `spatial_scatter_top_*.png`: spatial scatter plots showing top 6 genes among the highest Moran's I scores.
  - `spatial_scatter_bot_*.png`: spatial scatter plots showing top 6 genes among the lowest Moran's I scores.
- `parameters.json`: a JSON file containing the parameter settings in the analysis.


### <u>Annotation analysis</u>

**Output directory: `results/annotation/celltypist`**
- `adata_annotation.h5ad`: AnnData object file after cell-type annotation analysis.
- `sample_*/` or `group_*/`
  - `umap_cell_type.png`: UMAP plots showing predicted cell-type clusters.
  - `umap_conf_score.png`: UMAP plots showing mapped confidence scores of the cells.
  - `spatial_scatter_*.png`: spatial scatter plot shows how cell-types are spatially mapped onto the tissue morphology.
- `prop_majority_voting.png`: plot showing a stacked bar chart that presents the proportions of cell-type clusters across samples/groups.
- `parameters.json`: a JSON file containing the parameter settings in the analysis.


**Output directory: `results/annotation/scvi`**
- `adata_annotation.h5ad`: AnnData object file after cell-type annotation analysis.
- `sample_*/` or `group_*/`
  - `umap_cell_type.png`: UMAP plots showing predicted cell-type clusters.
  - `umap_conf_score.png`: UMAP plots showing mapped confidence scores of the cells.
  - `spatial_scatter_all_labels.png`: spatial scatter plot shows how all cell-types are spatially mapped onto the tissue morphology.
  - `spatial_scatter_label_*.png`: spatial scatter plot highlight how a single cell-type are spatially mapped onto the tissue morphology.
- `prop_scanvi_label.png`: plot showing a stacked bar chart that presents the proportions of cell-type clusters across samples/groups.
- `parameters.json`: a JSON file containing the parameter settings in the analysis.


### <u>DEA analysis</u>

**Output directory: `results/dea`**
- `markers/`: subfolder for the identification of marker genes.
- `compare/`: subfolder for the identification of differentially expressed genes between groups across all cells.
- `compare_ct/`: subfolder for the identification of differentially expressed genes between groups across cell types.
- `sample_*/` or `group_*/` or `celltype_*/` (no subfolder for DEA betweeen groups)
  - `plot_genes_*.png`: plots showing top number of DE genes across groups.
  - `dotplot_genes_*.png`: dot plot showing top number of DE genes across groups.
  - `dea_*.csv`: a csv table file showing DEA results for all genes, e.g. log fold change, p-values.
  - `spatial_scatter_*.png`: spatial scatter plots show top marker genes onto the tissue morphology.
- `parameters.json`: a JSON file containing the parameter settings in the analysis.


### <u>CellChat analysis</u>

**Output directory: `results/cellchat`**
- `sample_*/` or `group_*/`
  - `aggregated_network_all.png`: circular network plot showing aggregated cell-cell communications.
  - `aggregated_network_all_weights.png`: circular network plot showing aggregated cell-cell communications with total interaction strength (weights) between any two cell groups.
  - `aggregated_network_groups.png`: circular network plot (with weights) showing aggregated cell-cell communications sent from each cell group.
  - `pathway_network_circle_*.png`: circular network plot showing aggregated cell-cell communications for a signaling pathway.
  - `pathway_network_chord_*.png`: chord diagram plot showing aggregated cell-cell communications for a signaling pathway.
  - `pathway_network_heatmap_*.png`: heatmap plot showing aggregated cell-cell communications for a signaling pathway.
  - `pathway_network_contribution_*.png`: bar plot showing the contribution of each ligand-receptor pair to a signaling pathway.
  - `pathway_network_LR_circle_*.png`: circular network plot showing aggregated cell-cell communications mediated by a single ligand-receptor pair for a signaling pathway.
  - `pathway_network_LR_chord_*.png`: chord diagram plot showing aggregated cell-cell communications mediated by a single ligand-receptor pair for a signaling pathway.
  - `cellcell_LR_bubble_*.png`: bubble plot showing all the significant interactions (L-R pairs) from each cell group to other cell groups.
  - `cellcell_LR_chord_*.png`: chord diagram plot showing all the significant interactions (L-R pairs) from each cell group to other cell groups.
  - `pathway_genes_violin_*.png`: violin plot showing the gene expression distribution of signaling genes related to the inferred significant communications for a signaling pathway.
  - `pathway_network_centrality_*.png`: heatmap plot showing dominant senders, receivers, mediators and influencers in the intercellular communication network by computing several network centrality measures for cell groups of a signaling pathway.
  - `heatmap_signaling_patterns.png`: heatmap plot showing which signals contributing most to outgoing or incoming signaling of certain cell groups. In this heatmap, colobar represents the relative signaling strength of a signaling pathway across cell groups (NB: values are row-scaled). The top colored bar plot shows the total signaling strength of a cell group by summarizing all signaling pathways displayed in the heatmap. The right grey bar plot shows the total signaling strength of a signaling pathway by summarizing all cell groups displayed in the heatmap.
  - `inferred_cellcell_comm.rds`: a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors.
  - `cellchat.rds`: the CellChat object created for the analysis.
- `parameters.json`: a JSON file containing the parameter settings in the analysis.



## Pipeline reporting

### <u>Analysis report</u>

**Output directory: `results/report`**
- `eisca_report.html`: this is HTML report file showing all major analysis results.


### <u>MultiQC</u>

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### <u>Pipeline information</u>

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
