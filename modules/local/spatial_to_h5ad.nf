process SPATIAL_TO_H5AD {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::scanpy conda-forge::python-igraph conda-forge::leidenalg"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/villadsw/scanpy_squidpy:latest' :
        'docker.io/villadsw/scanpy_squidpy:latest' }"

    input:
    path outdir
    tuple val(meta), path(counts), path(metadata), path(transcripts)
    // tuple val(meta), path(counts)
    // tuple val(meta), path(metadata)

    output:
    tuple val("raw"), path("*.h5ad"),  emit: h5ad
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def transformation, datadir
    if (params.technology =='vizgen') {
        // transformation = "${data}/images/micron_to_mosaic_pixel_transform.csv"
        transformation = "micron_to_mosaic_pixel_transform.csv"
        datadir = "${outdir}/${params.technology}/${meta.id}"
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tech = params.technology
    // def transformation = "${data}/images/micron_to_mosaic_pixel_transform.csv"
    
    // run script
    """
    spatial_to_h5ad.py \\
        --tech ${tech} \\
        --datadir ${datadir} \\
        --counts ${counts} \\
        --metadata ${metadata} \\
        --transcripts ${transcripts} \\
        --transformation micron_to_mosaic_pixel_transform.csv \\
        --outfile "${meta.id}_st_matrix.h5ad" \\
        $args \\


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    mkdir ${meta.id}
    touch ${meta.id}/${meta.id}_matrix.h5ad
    touch versions.yml
    """
}
