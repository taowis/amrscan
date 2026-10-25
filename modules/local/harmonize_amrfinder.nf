process HARMONIZE_AMRFINDER {
    tag "${sample_id}"
    publishDir "${params.outdir}/harmonized", mode: 'copy'

    input:
    tuple val(sample_id), path(raw), path(provenance)

    output:
    tuple val(sample_id), path("${sample_id}.amr.tsv"), emit: results

    script:
    // All workflow-generated sample IDs have already been validated.
    """
    amrscan.py harmonize --sample '${sample_id}' \
        --raw '${raw}' --provenance '${provenance}' --output '${sample_id}.amr.tsv'
    """
}
