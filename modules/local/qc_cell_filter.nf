process QC_CELL_FILTER {
    label 'process_high'

    conda "conda-forge::scanpy conda-forge::python-igraph conda-forge::leidenalg"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/villadsw/scanpy_squidpy:latest' :
        'docker.io/villadsw/scanpy_squidpy:latest' }"

    input:
    path h5ad_raw
    path samplesheet

    output:
    path "qc_cell_filter"
    path "qc_cell_filter/*.h5ad",  emit: h5ad
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    if (params.technology == "vizgen")
    """
    qc_cell_filter.py \\
        --h5ad ${h5ad_raw} \\
        --samplesheet ${samplesheet} \\
        --outdir qc_cell_filter \\
        $args \\


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS        
    """
    
    else if (params.technology == "xenium")
    """
    qc_cell_filter_xenium.py \\
        --h5ad ${h5ad_raw} \\
        --samplesheet ${samplesheet} \\
        --outdir qc_cell_filter \\
        $args \\


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS        
    """


    // stub:
    // """
    // touch combined_matrix.h5ad
    // """
}
