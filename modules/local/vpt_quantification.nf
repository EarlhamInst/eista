process VPT_QUANTIFICATION {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vpt:1.3' :
        'docker.io/vzgdocker/vpt:1.3' }"

    input:
    tuple val(meta), path(data)
    path vzg

    output:
    tuple val(meta), path("cellpose_micron_space.parquet"), emit: parquet
    tuple val(meta), path("cell_by_gene.csv"), emit: counts
    tuple val(meta), path("detected_transcripts.csv"), emit: transcripts
    tuple val(meta), path("cell_metadata.csv"), emit: metadata
    tuple val(meta), path("sum_signals.csv"), emit: signals
    path "updated_${vzg ? vzg.getName() : '.vzg'}", optional: true
    path  "versions.yml", emit: versions
    path "micron_to_mosaic_pixel_transform.csv"
    tuple val(meta), path("cell_by_gene.csv"), path("cell_metadata.csv"), path("detected_transcripts.csv"), emit: datacountsmeta

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // def algorithm = (params.segmentation_algorithm)?: new File("${workflow.projectDir}/assets/algorithms/cellpose_default_1_ZLevel.json").text
    def algorithm = (params.segmentation_algorithm)?: "${workflow.projectDir}/assets/algorithms/cellpose_default_1_ZLevel.json"
    // def images = new File("${data}/mosaic_(?P<stain>[\w|-]+)_z(?P<z>[0-9]+).tif").text
    def images = """${data}/images/mosaic_(?P<stain>[\\w|-]+)_z(?P<z>[0-9]+).tif"""
    // def images = "${data}/images/mosaic_*_z*.tif"
    def transformation = "${data}/images/micron_to_mosaic_pixel_transform.csv"
    // if(!new File(transformation).exists()) {transformation = "${data}/micron_to_mosaic_pixel_transform.csv"}
    def transcripts = "${data}/detected_transcripts.csv"
    def vzg_update = "updated_${vzg ? vzg.getName() : '.vzg'}"

    // """
    // export HOME=$PWD
    // export CELLPOSE_LOCAL_MODELS_PATH=/ei/projects/0/05407428-a659-41d6-a7bd-4567bf45e494/data/vizgen/models
    // export CELLPOSE_LOCAL_MODELS_PATH=/home/wuhh/.cellpose/models

    """

    cp -n ${transformation} .

    vpt \\
        --verbose --processes ${task.cpus} run-segmentation \\
        --segmentation-algorithm ${algorithm} \\
        --input-images="${images}" \\
        --input-micron-to-mosaic ${transformation} \\
        --output-path "./" \\
        --tile-size ${params.tile_size} \\
        --tile-overlap ${params.tile_overlap} \\


    vpt \\
        --verbose partition-transcripts \\
        --input-boundaries cellpose_micron_space.parquet \\
        --input-transcripts ${transcripts} \\
        --output-entity-by-gene cell_by_gene.csv \\
        --output-transcripts detected_transcripts.csv \\


    vpt \\
        --verbose derive-entity-metadata \\
        --input-boundaries cellpose_micron_space.parquet \\
        --input-entity-by-gene cell_by_gene.csv \\
        --output-metadata cell_metadata.csv \\


    vpt \\
        --verbose sum-signals \\
        --input-images="${images}" \\
        --input-boundaries cellpose_micron_space.parquet \\
        --input-micron-to-mosaic ${transformation} \\
        --output-csv sum_signals.csv \\

           
    if [[ -f "${vzg}" ]]; then
        vpt \\
            --verbose update-vzg \\
            --input-vzg ${vzg} \\
            --input-boundaries cellpose_micron_space.parquet \\
            --input-entity-by-gene cell_by_gene.csv \\
            --input-metadata cell_metadata.csv \\
            --output-vzg ${vzg_update}
    else
        echo "File ${vzg} does not exist, skip update-vzg."
    fi  
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vpt: 1.3
    END_VERSIONS
    """

}