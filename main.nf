#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
include { AMRSCAN } from './workflows/amrscan'

workflow {
    // Types, ranges and paths are checked against nextflow_schema.json; --help
    // is printed from the same file.
    validateParameters()
    log.info paramsSummaryLog(workflow)

    def identity = params.amrfinder_identity as double
    if (!(identity == -1 || (identity >= 0 && identity <= 1))) error '--amrfinder_identity must be -1 or between 0 and 1.'
    def database = file(params.amrfinder_db, checkIfExists: true)

    if (params.input.toString().endsWith('.csv')) {
        // Samplesheet: columns, paths, IDs and uniqueness are checked by assets/schema_input.json.
        samples = Channel.fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
    } else {
        samples = fastaChannel(params.input)
    }
    if (params.validate_only) {
        samples.view { meta, fasta -> "Validated ${meta.id}: ${fasta} (caller NOT executed)" }
    } else {
        AMRSCAN(samples, Channel.value(database))
        // One entry per process: every task writes its own, identical versions.yml.
        AMRSCAN.out.versions
            .unique { it.text }
            .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info", sort: true)
    }
}

// FASTA file or glob: the sample ID is the file name without its suffixes.
def fastaChannel(pattern) {
    // Collect before scheduling: reject duplicate IDs before any output can collide.
    Channel.fromPath(pattern, checkIfExists: true, type: 'file')
        .toList()
        .flatMap { files ->
            def seen = [] as Set
            files.sort().collect { fasta ->
                if (!(fasta.name ==~ /.+\.(fa|fna|fasta)(\.gz)?/)) {
                    error "Expected assembled FASTA (.fa/.fna/.fasta, optionally .gz): ${fasta.name}"
                }
                def sample = fasta.name.replaceFirst(/\.(fa|fna|fasta)(\.gz)?$/, '')
                if (!(sample ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) error "Unsafe sample ID: ${sample}"
                if (!seen.add(sample)) error "Duplicate sample ID: ${sample}"
                tuple([id: sample], fasta)
            }
        }
}
