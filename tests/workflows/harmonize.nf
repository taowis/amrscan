nextflow.enable.dsl = 2
include { HARMONIZE_AMRFINDER } from '../../modules/local/harmonize_amrfinder'

workflow {
    HARMONIZE_AMRFINDER(Channel.of(tuple('upstream',
        file(params.raw, checkIfExists: true),
        file(params.provenance, checkIfExists: true))))
}
