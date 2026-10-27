process HARMONIZE_AMRFINDER {
    tag "${meta.id}"
    label 'process_single'
    publishDir "${params.outdir}/harmonized", mode: 'copy', pattern: "${meta.id}.amr.tsv"

    input:
    tuple val(meta), path(raw), path(provenance)

    output:
    tuple val(meta), path("${meta.id}.amr.tsv"), emit: results
    path 'versions.yml', emit: versions

    script:
    // All workflow-generated sample IDs have already been validated.
    """
    amrscan.py harmonize --sample '${meta.id}' \
        --raw '${raw}' --provenance '${provenance}' --output '${meta.id}.amr.tsv'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
