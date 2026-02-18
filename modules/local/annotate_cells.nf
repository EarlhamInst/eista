process ANNOTATE_CELLS {
    label 'process_medium'

    conda "conda-forge::scanpy conda-forge::python-igraph conda-forge::leidenalg"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/celltypist:1.6.3--pyhdfd78af_0' :
        'docker.io/myeihub/celltypist_squidpy:1.6.3_1.6.1' }"

    input:
    path h5ad_filtered
    // path model_file

    output:
    path "annotation/celltypist"
    path "annotation/celltypist/*.h5ad",  emit: h5ad
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    if(params.ctmodel)
    """
    export HOME=$PWD
    annotate_cells.py \\
        --h5ad ${h5ad_filtered} \\
        --outdir annotation/celltypist \\
        $args \\


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS        
    """
    
    else
    """
    export HOME=$PWD
    annotate_cells.py \\
        --h5ad ${h5ad_filtered} \\
        --outdir annotation/celltypist \\
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
