process AMRFINDERPLUS {
    tag "${sample_id}"
    publishDir "${params.outdir}/raw/amrfinderplus/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)
    path database, stageAs: 'amrfinder_db'

    output:
    tuple val(sample_id), path("${sample_id}.amrfinder.tsv"), path("${sample_id}.provenance.json"), emit: calls
    path "${sample_id}.amrfinder.log", emit: log

    script:
    def quote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
    def organism = params.amrfinder_organism ? "--organism ${quote(params.amrfinder_organism)}" : ''
    def plus = params.amrfinder_plus ? '--plus' : ''
    """
    amrscan.py call \
        --sample ${quote(sample_id)} --input ${quote(assembly)} \
        --database ${quote(database)} --database-source ${quote(params.amrfinder_db)} \
        --executable ${quote(params.amrfinder_bin)} --threads ${task.cpus} \
        --identity ${params.amrfinder_identity} --coverage ${params.amrfinder_coverage} \
        --raw ${quote(sample_id + '.amrfinder.tsv')} \
        --provenance ${quote(sample_id + '.provenance.json')} \
        --log ${quote(sample_id + '.amrfinder.log')} ${organism} ${plus}
    """
}
