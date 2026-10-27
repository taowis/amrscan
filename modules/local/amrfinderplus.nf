process AMRFINDERPLUS {
    tag "${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/raw/amrfinderplus/${meta.id}", mode: 'copy', pattern: "${meta.id}.*"

    input:
    tuple val(meta), path(assembly)
    path database, stageAs: 'amrfinder_db'

    output:
    tuple val(meta), path("${meta.id}.amrfinder.tsv"), path("${meta.id}.provenance.json"), emit: calls
    path "${meta.id}.amrfinder.log", emit: log
    path 'versions.yml', emit: versions

    script:
    def quote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
    def organism = params.amrfinder_organism ? "--organism ${quote(params.amrfinder_organism)}" : ''
    def plus = params.amrfinder_plus ? '--plus' : ''
    """
    amrscan.py call \
        --sample ${quote(meta.id)} --input ${quote(assembly)} \
        --database ${quote(database)} --database-source ${quote(params.amrfinder_db)} \
        --executable ${quote(params.amrfinder_bin)} --threads ${task.cpus} \
        --identity ${params.amrfinder_identity} --coverage ${params.amrfinder_coverage} \
        --raw ${quote(meta.id + '.amrfinder.tsv')} \
        --provenance ${quote(meta.id + '.provenance.json')} \
        --log ${quote(meta.id + '.amrfinder.log')} ${organism} ${plus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: \$(${quote(params.amrfinder_bin)} --version)
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
