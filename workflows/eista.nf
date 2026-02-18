/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_eista_pipeline'

include { VPT_QUANTIFICATION } from '../modules/local/vpt_quantification'
// include { VPT_SEGMENTATION } from '../modules/local/vpt_segmentation'
// include { VPT_PARTITION } from '../modules/local/vpt_partition'
// include { VPT_METADATA } from '../modules/local/vpt_metadata'
// include { VPT_SUM_SIGNALS } from '../modules/local/vpt_sum_signals'
// include { VPT_UPDATE_VZG } from '../modules/local/vpt_update_vzg'
// include { MTX_CONVERSION } from "../subworkflows/local/mtx_conversion"
include { SPATIAL_TO_H5AD } from '../modules/local/spatial_to_h5ad'
include { SPATIAL_TO_H5AD_XENIUM } from '../modules/local/spatial_to_h5ad_xenium'
include { CONCAT_H5AD } from '../modules/local/concat_h5ad'
include { QC_CELL_FILTER } from "../modules/local/qc_cell_filter"
include { CLUSTERING_ANALYSIS } from "../modules/local/clustering_analysis"
include { SPATIAL_STATISTICS } from "../modules/local/spatial_statistics"
include { ANNOTATE_CELLS } from '../modules/local/annotate_cells'
include { ANNOTATE_CELLS_SCVI } from '../modules/local/annotate_cells_scvi'
include { TRAIN_CT_MODEL } from '../modules/local/train_ct_model'
include { RANK_GENES } from '../modules/local/rank_genes'
include { H5AD_TO_MTX_META } from '../modules/local/h5ad_to_mtx_meta'
include { CELLCELL_COMMUNICATION } from '../modules/local/cellcell_communication'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow EISTA {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:
    samplesheet = file(params.input)
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_datacountsmeta = Channel.empty()
    ch_concat_h5ad_input = Channel.empty()

    //
    // MODULE: Run FastQC
    //
    // FASTQC (
    //     ch_samplesheet
    // )
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    // ch_versions = ch_versions.mix(FASTQC.out.versions.first())



    //===================================== Primary anaysis stage =====================================

    if (params.run_analyses.contains('primary')){
        if (params.technology == "vizgen") {
            if (!params.skip_analyses.contains('quantification')) {

                VPT_QUANTIFICATION (
                    ch_samplesheet,
                    ch_samplesheet.map{ meta, data -> file("${data[0]}/*.vzg") }
                )
                ch_versions = ch_versions.mix(VPT_QUANTIFICATION.out.versions)
                VPT_QUANTIFICATION.out.datacountsmeta.set{ch_datacountsmeta}
                // ch_counts = VPT_QUANTIFICATION.out.counts
                // ch_metadata = VPT_QUANTIFICATION.out.metadata
            }
            SPATIAL_TO_H5AD (
                Channel.value(file(params.outdir)),
                ch_datacountsmeta
            )
            // ch_versions = ch_versions.mix(SPATIAL_TO_H5AD.out.versions)
            ch_concat_h5ad_input = SPATIAL_TO_H5AD.out.h5ad.groupTuple()

        } else if (params.technology == "xenium") {

            SPATIAL_TO_H5AD_XENIUM (
                Channel.value(file(params.outdir)),
                ch_samplesheet
            )
            // ch_versions = ch_versions.mix(SPATIAL_TO_H5AD_XENIUM.out.versions)
            ch_concat_h5ad_input = SPATIAL_TO_H5AD_XENIUM.out.h5ad.groupTuple()
        }

        // Run mtx to h5ad conversion subworkflow
        // MTX_CONVERSION (
        //     ch_samplesheet,
        //     ch_counts,
        //     ch_metadata,
        //     samplesheet
        // )
        // ch_versions.mix(MTX_CONVERSION.out.ch_versions)

        // ch_samplesheet
        // .join(ch_counts, by: [0]).map{ meta, data, counts -> tuple(meta, data, counts) }
        // .join(ch_metadata, by: [0]).map{ meta, data, counts, metadata -> tuple(meta, data, counts, metadata) }
        // .set { ch_data_counts_meta }
        
        CONCAT_H5AD (
            ch_concat_h5ad_input,
            samplesheet
        )
        ch_versions = ch_versions.mix(CONCAT_H5AD.out.versions)  

    }

    //===================================== Secondary anaysis stage =====================================

    // if (params.run_analyses.contains('secondary')){
    if (params.run_analyses.any{ it in ['secondary', 'qccellfilter', 'clustering', 'spatialstats'] }){
    
        // MODULE: Run QC and cell filtering
        ch_h5ad = Channel.empty()
        if(params.run_analyses.contains('primary')){
            ch_h5ad = CONCAT_H5AD.out.h5ad
        }else if(params.h5ad){
            ch_h5ad = Channel.fromPath(params.h5ad)
        }else{
            def path1 = "${params.outdir}/vizgen/mtx_conversions/combined_st_matrix.h5ad"
            def path2 = "${params.outdir}/qc_cell_filter/adata_filtered_normalized.h5ad"
            def path3 = "${params.outdir}/clustering/adata_clustering.h5ad"
            if (params.run_analyses.any{it=='secondary' || it=='qccellfilter'} && 
                !params.skip_analyses.contains('qccellfilter') && new File(path1).exists()) {
                ch_h5ad = Channel.fromPath(path1)
            }else if (params.run_analyses.any{it=='secondary' || it=='clustering'} && 
                !params.skip_analyses.contains('clustering') && new File(path2).exists()) {
                ch_h5ad = Channel.fromPath(path2)
            }else if (params.run_analyses.any{it=='secondary' || it=='spatialstats'} && 
                !params.skip_analyses.contains('spatialstats') && new File(path3).exists()) {
                ch_h5ad = Channel.fromPath(path3)
            }else{
                log.warn("No specified h5ad file found in previous analysis, please specify an h5ad file by setting --h5ad.")
                return
            }
        }

        if (params.run_analyses.any{it=='secondary' || it=='qccellfilter'} && !params.skip_analyses.contains('qccellfilter')) {
            QC_CELL_FILTER (
                ch_h5ad,
                Channel.fromPath(params.input)
            )
            // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
            ch_versions = ch_versions.mix(QC_CELL_FILTER.out.versions)
            ch_h5ad = QC_CELL_FILTER.out.h5ad         
        }
        
        if (params.run_analyses.any{it=='secondary' || it=='clustering'} && !params.skip_analyses.contains('clustering')) {
            CLUSTERING_ANALYSIS (
                ch_h5ad
            )
            ch_versions = ch_versions.mix(CLUSTERING_ANALYSIS.out.versions)
            ch_h5ad = CLUSTERING_ANALYSIS.out.h5ad
        } 

        if (params.run_analyses.any{it=='secondary' || it=='spatialstats'} && !params.skip_analyses.contains('spatialstats')) {
            SPATIAL_STATISTICS (
                ch_h5ad
            )
            ch_versions = ch_versions.mix(SPATIAL_STATISTICS.out.versions)
        }          
        
    }


    //===================================== Tertiary anaysis stage =====================================

    if (params.run_analyses.any{ it in ['tertiary', 'annotation', 'annotation_scvi', 'dea', 'cellchat'] }){
    
        // Get input h5ad file
        ch_h5ad = Channel.empty()
        if(params.run_analyses.contains('secondary')){
            if (!params.skip_analyses.contains('clustering')) {
                ch_h5ad = CLUSTERING_ANALYSIS.out.h5ad
            }else {
                ch_h5ad = QC_CELL_FILTER.out.h5ad
            }
        }else if(params.h5ad){
            ch_h5ad = Channel.fromPath(params.h5ad)
        }else{
            def path1 = "${params.outdir}/clustering/adata_clustering.h5ad"
            def path2 = "${params.outdir}/qc_cell_filter/adata_filtered_normalized.h5ad"
            def path3 = "${params.outdir}/annotation/adata_annotation.h5ad"
            def path4 = "${params.outdir}/annotation_scvi/adata_annotation.h5ad"
            if(params.run_analyses.any{it in ['dea', 'cellchat']}){
                if(new File(path3).exists()){
                    ch_h5ad = Channel.fromPath(path3)
                }else if(new File(path4).exists()){
                    ch_h5ad = Channel.fromPath(path4)
                }else if(new File(path1).exists()){
                    ch_h5ad = Channel.fromPath(path1)
                }          
            }else if(new File(path1).exists()){
                ch_h5ad = Channel.fromPath(path1)
            }else if(new File(path2).exists()){
                ch_h5ad = Channel.fromPath(path2)
            }
        }
        ch_h5ad.ifEmpty {
            log.warn("No specified h5ad file found in secondary analysis, please specify an h5ad file by setting --h5ad.")
            return            
        }

        if (params.run_analyses.any{it=='tertiary' || it=='annotation'} and !params.skip_analyses.contains('annotation')) {
            // ch_ctmodel = params.ctmodel? Channel.fromPath(params.ctmodel) : []
            ANNOTATE_CELLS (
                ch_h5ad,
                // ch_ctmodel
                // MTX_CONVERSION.out.h5ad
            )
            // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
            ch_versions = ch_versions.mix(ANNOTATE_CELLS.out.versions)
            ch_h5ad = ANNOTATE_CELLS.out.h5ad      
        }
        
        if (params.run_analyses.any{it=='tertiary' || it=='annotation_scvi'} && 
            !params.skip_analyses.contains('annotation_scvi') && params.args_annotation_scvi) {
            // ch_h5ad_ref = params.h5ad_ref? Channel.fromPath(params.h5ad_ref) : []
            ANNOTATE_CELLS_SCVI (
                ch_h5ad,
                // ch_h5ad_ref,
                // ch_ctmodel
                // MTX_CONVERSION.out.h5ad
            )
            // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
            ch_versions = ch_versions.mix(ANNOTATE_CELLS_SCVI.out.versions)
            ch_h5ad = ANNOTATE_CELLS_SCVI.out.h5ad      
        }

        if (params.run_analyses.any{it=='tertiary' || it=='dea'} and !params.skip_analyses.contains('dea')) {
            RANK_GENES (
                ch_h5ad,
            )
            ch_versions = ch_versions.mix(RANK_GENES.out.versions)
        }

        if (params.run_analyses.any{it=='tertiary' || it=='cellchat'} and !params.skip_analyses.contains('cellchat')) {
            H5AD_TO_MTX_META (
                ch_h5ad,
            )
            ch_versions = ch_versions.mix(H5AD_TO_MTX_META.out.versions)
            CELLCELL_COMMUNICATION (
                H5AD_TO_MTX_META.out.h5ad_mtx,
                H5AD_TO_MTX_META.out.h5ad_meta,
                H5AD_TO_MTX_META.out.h5ad_genes,
                H5AD_TO_MTX_META.out.h5ad_cells,
            )
            ch_versions = ch_versions.mix(CELLCELL_COMMUNICATION.out.versions)
        }

        
    }

    // train cell-type models for CellTypist
    if (params.run_analyses.contains('ctmodel') and !params.skip_analyses.contains('ctmodel')) {
            ch_h5ad = Channel.fromPath(params.h5ad)
            TRAIN_CT_MODEL (
                ch_h5ad,
                // MTX_CONVERSION.out.h5ad
            )
            // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
            ch_versions = ch_versions.mix(TRAIN_CT_MODEL.out.versions)
            // ch_h5ad = ANNOTATE_CELLS.out.h5ad      
    }





    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
