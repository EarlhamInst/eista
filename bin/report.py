#!/usr/bin/env python

import os
os.environ[ 'MPLCONFIGDIR' ] = '/tmp/'

import dominate.tags as html
import ezcharts as ezc
from ezcharts.components.reports.labs import LabsReport, LabsNavigation, ILabsNavigationClasses
from ezcharts.layout.snippets import DataTable, Grid, Tabs
from ezcharts.components.theme import LAB_head_resources
from ezcharts.components.ezchart import EZChart

import argparse
from pathlib import Path
from report_util import *
import util
import json
import pandas as pd
import sys


logger = util.get_named_logger('REPORT')

report_title = 'EI Spatial Transcriptomics Analysis Report'
workflow_name = 'eista'


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""

    parser = argparse.ArgumentParser(
        description="Quality control before and after cell filtering",
        # epilog="python count_reads_from_bam.py --bam file.bam --bed file.bed --json output.json",
    )
    parser.add_argument(
        "--report_html", 
        help="Report output file"
    )
    parser.add_argument(
        "--results",
        metavar="RESULTS_DIR",
        type=Path,
        help="Results directory.",
        required=True,
    )
    # parser.add_argument(
    #     "--stats", nargs='+',
    #     help="Fastcat per-read stats, ordered as per entries in --metadata.")
    # parser.add_argument(
    #     "--images", nargs='+',
    #     help="Sample directories containing various images to put in report")
    # parser.add_argument(
    #     "--survival",
    #     help="Read survival data in TSV format")
    parser.add_argument(
        "--params",
        metavar="FILE_PARAMS",
        type=Path,
        help="Workflow params json file",
        required=True,
    )
    parser.add_argument(
        "--versions",
        metavar="FILE_VERSIONS",
        type=Path,
        help="Workflow versions file",
        required=True,
    )
    parser.add_argument(
        "--logo",
        metavar="FILE_LOGO",
        type=Path,
        help="Logo image file",
        required=True,
    )    
    # parser.add_argument(
    #     "--umap_dirs", nargs='+',
    #     help="Sample directories containing umap and gene expression files")
    # parser.add_argument(
    #     "--umap_genes", help="File containing list of genes to annnotate UMAPs")
    # parser.add_argument(
    #     "--metadata", default='metadata.json', required=True,
    #     help="sample metadata")
    parser.add_argument(
        "--wf-version", default='unknown',
        help="version of the executed workflow")
    parser.add_argument(
        "--samplesheet",
        metavar="SAMPLESHEET",
        type=Path,
        help="Input samplesheet file.",
        required=True,
    )
    parser.add_argument(
        "--tech",
        help="The spatial single-cell technology.",
        required=True,
    )           
    return parser.parse_args(argv)


def main(argv=None):
    logger.info('Building report')
    args = parse_args(argv)


    if not args.params.is_file():
        logger.error(f"The given input file {args.params} was not found!")
        sys.exit(2)

    if not args.versions.is_file():
        logger.error(f"The given input file {args.versions} was not found!")
        sys.exit(2)

    report = LabsReport(
        report_title, workflow_name,
        args.params, args.versions, args.wf_version,
        head_resources=[*LAB_head_resources],
    )
    
    report.nav.getElementsByTagName('a')[0].clear()
    navdiv = report.nav.getElementsByTagName('div')[0]
    link = html.a(href='https://www.earlham.ac.uk/', cls=ILabsNavigationClasses().logo)
    with link:
        # with open(Path('images/EI_logo.png'), 'rb') as f:
        with open(Path(args.logo), 'rb') as f:
            b64img = base64.b64encode(f.read()).decode()
            html.img(src=f'data:image/png;base64,{b64img}', width=120)
    navdiv[0] = link

    report.banner.clear()
    report.footer.clear()
    report.intro_content.add(EIBanner(report_title, workflow_name))
    report.footer.add(EILabsAddendum(workflow_name, args.wf_version))

    samplesheet = pd.read_csv(args.samplesheet)
    # batch = 'group' if 'group' in samplesheet.columns else 'sample'
    sample = 'plate' if 'plate' in samplesheet.columns else  'sample'
    batch = 'sample'
    if 'group' in samplesheet.columns:
        batch = 'group'
    elif 'plate' in samplesheet.columns:
        batch = 'plate'  
    # Nbatch = len(samplesheet[batch].unique())
    # samples = samplesheet['sample'].unique()

    cell_size_unit = 'volume'
    if args.tech == 'xenium': cell_size_unit = 'area'

    path_quant_qc = Path(args.results, 'qc_cell_filter')
    path_quant_qc_raw = Path(path_quant_qc, 'raw_counts')
    # path_quant_qc_scatter = Path(path_quant_qc, 'scatter')
    # path_quant_qc_violin = Path(path_quant_qc, 'violin')
    path_cell_filtering = Path(path_quant_qc, 'cell_filtering')
    # path_cell_filtering_dist = Path(path_cell_filtering, 'distribution')
    path_clustering = Path(args.results, 'clustering')
    path_spatialstats = Path(args.results, 'spatialstats')
    path_annotation = Path(args.results, 'annotation')
    path_dea = Path(args.results, 'dea')
    path_cellchat = Path(args.results, 'cellchat')

    # print(path_quant_qc) #tst

    if path_quant_qc.exists():
        summary = pd.read_csv(Path(path_quant_qc, 'sample_summary.csv')).set_index(f"{sample.capitalize()} ID")
        Nsample = summary.shape[0]
        with report.add_section('Single cell summary', 'Summary'):
            html.p(f"""This section gives an overall summary of the single-cell count matrix for 
                   each sample. The statistics include the total number of cells with at least 
                   one gene expressed, the total number of genes expressed in at least one cell, 
                   the median number of genes per cell, and the {cell_size_unit} of segmented cells.""")
            DataTable.from_pandas(summary)   

        with report.add_section('Quantification QC', 'QC'):
            html.p("""This section presents the QC plots of the raw count matrix generated during 
                   the quantification step. These plots provide insight into the quality of the 
                   experiments and guide the filtering of low-quality cells.""")
            html.p(f"""The following scatter plot shows the relationship between total 
                   read counts and the number of genes, with the {cell_size_unit} of segmented cells indicated by color.""")
            plots_from_image_files(path_quant_qc_raw, meta='sample', widths=['800'], suffix=['scatter*.png'])
            html.p("""The following knee plot shows the relationship between cell segmentation rank and 
                   total number of transcripts per barcode. The inflection (knee) point where the transcript 
                   frequency rapidly drops can be good cutoff to remove empty cell segmentation.""")
            plots_from_image_files(path_quant_qc_raw, meta='sample', widths=['500'], suffix=['knee_plot.png'])
            html.p(f"""The following violin plots display the distribution of cells based on the number of 
                   genes, total counts, and cell {cell_size_unit}s.""")
            plots_from_image_files(path_quant_qc_raw, meta='sample', ncol=3, suffix=['violin*.png'])
            html.p(f"""The following histogram plots display the distribution of cells based on the number of 
                   genes, total counts, and cell {cell_size_unit}s.""")
            plots_from_image_files(path_quant_qc_raw, meta='sample', suffix=['histograms.png'])
    else:
        logger.info('Skipping Quantification QC')

    if path_cell_filtering.exists():
        summary = pd.read_csv(Path(path_cell_filtering, 'sample_summary_filtered.csv')).set_index(f"{sample.capitalize()} ID")
        with report.add_section('Cell filtering', 'Cell filtering'):
            html.p(f"""This section presents the statistics and QC plots after cell filtering process. 
                   The filtering process includes hard thresholds for minimum genes, minimum counts, 
                   minimum cells and the {cell_size_unit} of segmented cells. Additionally, 
                   users can set quantile limits on the number of genes.""")
            html.p("""The following table shows summary statistics, with percentages in brackets 
                   indicating the comparison to the raw counts.""")
            DataTable.from_pandas(summary)
            html.p(f"""The following violin plots display the distribution of cells based on the number of 
                   genes, total counts, and cell {cell_size_unit}s after filtering.""")
            plots_from_image_files(path_cell_filtering, meta='sample', ncol=3, suffix=['violin*.png'])            
            html.p(f"""The following histogram plots display the distribution of cells based on the number of 
                   genes, total counts, and cell {cell_size_unit}s.""")
            plots_from_image_files(path_cell_filtering, meta='sample', suffix=['histograms.png'])
            # html.p("""The following plots show the UMAP plots for the number of genes, total counts.""")                        
            # plots_from_image_files(path_cell_filtering, meta='sample', suffix=['umap_total*.png'])
            html.p("""The following plots show the spatial scatter plots for the number of genes, total counts.""")                        
            plots_from_image_files(path_cell_filtering, meta='sample', suffix=['spatial_scatter*.png'])
            # if util.check_file(f"{path_cell_filtering}/sample_*", 'umap_doublet.png'):
            #     html.p("""The following plots show the UMAP plots for the predicted doublets and doublet scores.""")                        
            #     plots_from_image_files(path_cell_filtering, meta='sample', suffix=['umap_doublet.png'])
            show_analysis_parameters(f"{path_quant_qc}/parameters.json")          
    else:
        logger.info('Skipping Cell filtering')


    if path_clustering.exists():
        if util.check_file(f"{path_clustering}/sample_*", ''):
            batch = 'sample'
        elif util.check_file(f"{path_clustering}/group_*", ''):
            batch = 'group'
        Nbatch = len(samplesheet[batch].unique())
        if Nbatch == 1: Nbatch = 2
        with report.add_section('Clustering analysis', 'Clustering'):
            html.p(f"""This section shows clustering UMAP plots for each {batch}. The clustering 
                   was performed using Leiden graph-clustering method. The resolution parameter 
                   was set for different values to get different number of clusters which 
                   could match to biologically-meaningful cell types.""")
            plots_from_image_files(path_clustering, meta=batch, ncol=2, suffix=['umap*.png'])
            html.p("""The following plots show the spatial scatter plots for the corresponding resolutions.""")                        
            plots_from_image_files(path_clustering, meta=batch, ncol=2, suffix=['spatial_scatter*.png'])
            html.p(f"""The following plot shows a stacked bar chart that presents the proportions of clusters 
                   across {batch}s, calculated for each resolution value. The plot illustrates the distribution 
                   profiles of predicted clusters between {batch}s.""")                        
            plots_from_image_files(path_clustering, meta='resolution', widths=[str(min(Nbatch*280, 1200))])
            show_analysis_parameters(f"{path_clustering}/parameters.json")                   
    else:
        logger.info('Skipping clustering analysis')        


    if path_spatialstats.exists():
        if util.check_file(f"{path_spatialstats}/sample_*", ''):
            batch = 'sample'
        elif util.check_file(f"{path_spatialstats}/group_*", ''):
            batch = 'group'
        Nbatch = len(samplesheet[batch].unique())
        # pd_autocorr = pd.read_csv(Path(path_spatialstats, 'autocorr_moranI.csv'))
        with report.add_section('Spatial statistics analysis', 'Spatial-stats'):
            html.p(f"""This section presents the spatial statistics analysis results for each {batch}. The statistics, 
                   computed using the Squidpy package, include centrality scores, neighborhood enrichment scores, 
                   and Moran’s I score. These statistics indicate the relationship between 
                   expression patterns and biological morphology.""")
            
            html.p("""The following plots show the centrality analysis results of the 3 centrality scores include: 
                   1) Closeness centrality, which indicates how close a group is to other nodes. 2) Degree centrality, 
                   which represents the fraction of connected non-group members. 3) Clustering coefficient, which 
                   measures the degree to which nodes form clusters. The three centrality scores are visualized in the 
                   score plots. The spatial scatter plots highlight clusters with higher or lower values of these 
                   centrality scores, which can be explained as follows: clusters with high closeness centrality are 
                   close to other groups and tend to display a dispersed distribution throughout the tissue, while 
                   clusters with low closeness centrality tend to display an uneven and isolated distribution; clusters 
                   with high degree centrality scores have a high fraction of non-group member connections, while 
                   clusters with low degree centrality scores have fewer non-group member connections and tend to be 
                   lower-abundance clusters; a higher clustering coefficient indicates a stronger tendency for nodes 
                   to cluster together, while a lower clustering coefficient suggests a more even distribution 
                   throughout the tissue.""")                        
            plots_from_image_files(path_spatialstats, meta=batch, ncol=1, suffix=['centrality_scores*.png'])            
            plots_from_image_files(path_spatialstats, meta=batch, ncol=2, suffix=['spatial_scatter_closeness*.png'])            
            plots_from_image_files(path_spatialstats, meta=batch, ncol=2, suffix=['spatial_scatter_degree*.png'])            
            plots_from_image_files(path_spatialstats, meta=batch, ncol=2, suffix=['spatial_scatter_clustering*.png'])            

            html.p("""The following heatmap shows the results of the neighborhood enrichment analysis on clusters. This 
                   analysis indicates whether certain cell types or gene expression patterns are spatially enriched in 
                   neighboring regions of a tissue. It helps identify spatial dependencies between cells or tissue regions. 
                   From the heatmap, we can observe the strength of spatial co-localization between clusters/cell types 
                   in a tissue. A strong positive score indicates that they are frequently neighbors, while a strong 
                   negative score indicates that they avoid each other.""")                    
            plots_from_image_files(path_spatialstats, meta=batch, suffix=['neighbors_enrichment*.png'], widths=['700'])

            html.p("""The following table shows Moran’s I statistics for autocorrelation analysis of the top 50 and bottom 
                   50 genes. Here, 'I' represents Moran's score, 'pval_norm' is the p-value under the normality assumption, 
                   and 'var_norm' is the variance of the score under the normality assumption. Moran’s I global spatial 
                   autocorrelation statistic evaluates whether features (i.e., genes) exhibit spatial clustering (local 
                   patterns) or a random distribution across a spatial field. Genes with high spatial autocorrelation follow 
                   a clustered spatial pattern, while genes with low spatial autocorrelation display a more evenly 
                   distributed expression pattern.""")
            show_tab_table(path_spatialstats, meta=batch, suffix='autocorr_moranI.csv')
            # DataTable.from_pandas(pd_autocorr)
            html.p(f"""The following plots show the spatial scatter plots for top 6 genes among the highest Moran's I scores.""")                        
            plots_from_image_files(path_spatialstats, meta=batch, ncol=3, suffix=['spatial_scatter_top*.png'])
            html.p(f"""The following plots show the spatial scatter plots for bottom 6 genes among the lowest Moran's I scores.""")                        
            plots_from_image_files(path_spatialstats, meta=batch, ncol=3, suffix=['spatial_scatter_bot*.png'])

            show_analysis_parameters(f"{path_spatialstats}/parameters.json")                   
    else:
        logger.info('Skipping clustering analysis') 


    if path_annotation.exists():
        if util.check_file(f"{path_annotation}/sample_*", ''):
            batch = 'sample'
        elif util.check_file(f"{path_annotation}/group_*", ''):
            batch = 'group'
        Nbatch = len(samplesheet[batch].unique())
        if Nbatch == 1: Nbatch = 2
        with report.add_section('Cell-type annotation', 'Annotation'):
            html.p(f"""This section presents cell-type annotation results using CellTypist which is an 
            automated tool for cell type annotation based on pre-trained models, capable of accurately 
            classifying different cell types and subtypes.""")
            html.p("""The following UMAP plots show the predicted cell-type clusters and the mapped 
            confidence scores of the cells.""")            
            plots_from_image_files(path_annotation, suffix=['umap_cell_type.png'], meta=batch, widths=['1200'])
            plots_from_image_files(path_annotation, suffix=['umap_conf_score.png'], meta=batch, widths=['600'])
            html.p("""The following spatial scatter plot shows how cell-types are spatially mapped onto the tissue morphology.""") 
            plots_from_image_files(path_annotation, suffix=['spatial_scatter_*.png'], meta=batch, widths=['1000'])
            html.p(f"""The following plot shows a stacked bar chart that presents the proportions 
                   of cell-type clusters across {batch}s. The plot illustrates the distribution 
                   profiles of predicted cell-type clusters between {batch}s.""")                   
            plots_from_image_files(path_annotation, suffix=['prop_*.png'], widths=[str(min(Nbatch*500, 1200))])
            show_analysis_parameters(f"{path_annotation}/parameters.json")                 
    else:
        logger.info('Skipping cell-type annotation')   


    if path_dea.exists():
        path_dea_markers = Path(path_dea, 'markers')
        path_dea_compare = Path(path_dea, 'compare')
        path_dea_compare_ct = Path(path_dea, 'compare_ct')
        with report.add_section('Differential expression analysis', 'DEA'):
            if path_dea_markers.exists():
                if util.check_file(f"{path_dea_markers}/sample_*", ''):
                    batch = 'sample'
                elif util.check_file(f"{path_dea_markers}/group_*", ''):
                    batch = 'group'                
                html.p("""This section presents the results of the differential expression analysis performed using 
                        Scanpy’s rank_genes_groups function. These results enable the identification of marker genes 
                        by comparing genes that are highly ranked in one cluster against all other clusters.""")

                # showing plots for one cluster vs rest for each sample/group
                if util.check_file(f"{path_dea_markers}/{batch}_*", '*.png'):
                    html.p(f"""The following plots display the ranking of genes for one of the cell clusters against the rest of the clusters across {batch}s.""")                        
                    plots_from_image_files(path_dea_markers, meta=batch, suffix=['plot_genes_*.png'])
                    plots_from_image_files(path_dea_markers, meta=batch, suffix=['dotplot_genes_*.png'])
                    html.p("""The following spatial scatter plots show top marker genes onto the tissue morphology.""")
                    plots_from_image_files(path_dea_markers, meta=batch, ncol=4, suffix=['spatial_scatter_*.png'])

                show_analysis_parameters(f"{path_dea_markers}/parameters.json")
                if path_dea_compare.exists() or path_dea_compare_ct.exists(): html.hr(style="border: 1px solid grey;")

            if path_dea_compare.exists():       
                html.p("""This section presents the results of the differential expression analysis performed using 
                        Scanpy’s rank_genes_groups function. The analysis identifies differentially expressed genes at 
                        the cellular level by comparing one group of cells against another across all cells.""")

                # showing plots for DEA between conditions for all cells
                if util.check_file(f"{path_dea_compare}", '*.png'):
                    html.p("""The following plots show the ranking of differentially expressed genes of one group of 
                            cells against another.""")                        
                    plots_from_image_files(path_dea_compare, suffix=['plot_genes_*.png'])
                    html.div(style="height: 50px;")
                    plots_from_image_files(path_dea_compare, suffix=['dotplot_genes_*.png'])
                    html.div(style="height: 50px;")
                    html.p("""The following spatial scatter plots show top marker genes onto the tissue morphology.""")
                    plots_from_image_files(path_dea_compare, ncol=4, suffix=['spatial_scatter_*.png'])

                show_analysis_parameters(f"{path_dea_compare}/parameters.json")
                if path_dea_compare_ct.exists(): html.hr(style="border: 1px solid grey;")

            if path_dea_compare_ct.exists():
                if util.check_file(f"{path_dea_compare_ct}/sample_*", ''):
                    batch = 'sample'
                elif util.check_file(f"{path_dea_compare_ct}/group_*", ''):
                    batch = 'group'        
                html.p("""This section presents the results of the differential expression analysis performed using 
                        Scanpy’s rank_genes_groups function. The analysis identifies differentially expressed genes at 
                        the cellular level by comparing one group of cells against another across cell types/clusters.""")

                # showing plots for DEA between conditions for each celltype
                if util.check_file(f"{path_dea_compare_ct}/celltype_*", '*.png'):
                    html.p("""The following plots show differentially expressed genes between two conditions across cell types/clusters.""")                        
                    plots_from_image_files(path_dea_compare_ct, meta='celltype', suffix=['plot_genes_*.png'])
                    plots_from_image_files(path_dea_compare_ct, meta='celltype', suffix=['dotplot_genes_*.png'])
                    html.p("""The following spatial scatter plots show top marker genes onto the tissue morphology.""")
                    plots_from_image_files(path_dea_compare_ct, meta='celltype', ncol=4, suffix=['spatial_scatter_*.png'])

                show_analysis_parameters(f"{path_dea_compare_ct}/parameters.json")              
    else:
        logger.info('Skipping differential expression analysis')


    if path_cellchat.exists():
        if util.check_file(f"{path_cellchat}/sample_*", ''):
            batch = 'sample'
        elif util.check_file(f"{path_cellchat}/group_*", ''):
            batch = 'group'        
        with report.add_section('Cell-cell communication analysis', 'CellChat'):
            html.p("""This section presents the results of the cell–cell communication analysis using CellChat. These 
                   results enable users to identify, analyze, and visualize intercellular communication networks between 
                   cell groups from scRNA-seq data, based on the CellChat database of ligand–receptor interactions.""")

            # showing the aggregated cell-cell communication network
            if util.check_file(f"{path_cellchat}/{batch}_*", 'aggregated_network_all.png'):
                html.p(f"""The following plots show the aggregated cell-cell communication network, either by counting the 
                       number of links or by summarizing communication probabilities (weights). These circle plots provide 
                       an overview of cell-cell interactions: the left plot shows the number of interactions, while the right 
                       plot displays the total interaction strength between cell groups. Thicker edge lines indicate stronger 
                       signals, and the circle sizes are proportional to the number of cells in each group.""")                        
                plots_from_image_files(path_cellchat, meta=batch, ncol=2, suffix=['aggregated_network_all*.png'])

            # showing the aggregated cell-cell communication network for each signaling pathway
            if util.check_file(f"{path_cellchat}/{batch}_*", 'pathway_network_circle_*.png'):
                html.p(f"""The following plots show the aggregated (weighted) cell-cell communication networks for 
                       individual signalling pathways with significant ligand-receptor interactions.""")                        
                plots_from_image_files(path_cellchat, meta=batch, ncol=3, suffix=['pathway_network_circle_*.png'])

            # showing the aggregated cell-cell communication network for each signaling pathway (heatmap)
            # if util.check_file(f"{path_cellchat}/{batch}_*", 'pathway_network_heatmap_*.png'):
            #     html.p(f"""The following heatmaps show interaction strength in the cell-cell communication network between 
            #            pairs of cell groups, with the y-axis representing sender cells and the x-axis representing receiver 
            #            cells. The top bar plot shows the column-wise sums of the heatmap values, while the right bar plot 
            #            shows the row-wise sums.""")                        
            #     plots_from_image_files(path_cellchat, meta=batch, ncol=2, suffix=['pathway_network_heatmap_*.png'])

            # showing significant interactions (L-R pairs) from each cell group to other cell groups
            if util.check_file(f"{path_cellchat}/{batch}_*", 'cellcell_LR_bubble_*.png'):
                html.p(f"""The following bubble plots show all significant interactions (ligand-receptor pairs) from each 
                       cell group to others. The color and size of each dot represent the communication probability and 
                       p-value, respectively. This visualization shows how cell-cell communication is mediated by multiple 
                       ligand-receptor pairs.""")                        
                plots_from_image_files(path_cellchat, meta=batch, ncol=2, suffix=['cellcell_LR_bubble_*.png'])

            show_analysis_parameters(f"{path_cellchat}/parameters.json")                 
    else:
        logger.info('Skipping cell-cell communication analysis')


    report.write(args.report_html)
    logger.info('Report writing finished')
    

if __name__ == "__main__":
    sys.exit(main())
