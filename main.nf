#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { AMRSCAN } from './workflows/amrscan'

workflow {
    if (params.help) {
        log.info '''AMRScan v2: AMR determinant calling from assembled nucleotide FASTA
Usage: nextflow run main.nf --input 'assemblies/*.fna' --amrfinder_db /path/to/versioned/db
Options: --outdir results/v2 --threads 4 --amrfinder_bin amrfinder
         --amrfinder_organism Escherichia --amrfinder_plus
         --amrfinder_identity -1 --amrfinder_coverage 0.5
         --validate_only (validate configuration and input paths; do not run callers)
See README.md and docs/result-schema.md. Legacy: workflow/AMRScan.nf'''
    } else {
        if (!params.input) error 'Provide --input with assembled FASTA file(s).'
        if (!params.amrfinder_db) error 'Provide --amrfinder_db with a versioned AMRFinderPlus database directory.'
        if (!(params.threads.toString() ==~ /[1-9][0-9]*/)) error '--threads must be a positive integer.'
        def identity = params.amrfinder_identity.toString().toDouble()
        def coverage = params.amrfinder_coverage.toString().toDouble()
        if (!(identity == -1 || (identity >= 0 && identity <= 1))) error '--amrfinder_identity must be -1 or between 0 and 1.'
        if (!(coverage >= 0 && coverage <= 1)) error '--amrfinder_coverage must be between 0 and 1.'
        def database = file(params.amrfinder_db, checkIfExists: true)
        if (!database.isDirectory()) error '--amrfinder_db must be a directory.'

        // Collect before scheduling: reject duplicate IDs before any output can collide.
        samples = Channel.fromPath(params.input, checkIfExists: true, type: 'file')
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
                    tuple(sample, fasta)
                }
            }
        if (params.validate_only) {
            samples.view { sample, fasta -> "Validated ${sample}: ${fasta} (caller NOT executed)" }
        } else {
            AMRSCAN(samples, Channel.value(database))
        }
    }
}
