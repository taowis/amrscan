include { AMRFINDERPLUS } from '../modules/local/amrfinderplus'
include { HARMONIZE_AMRFINDER } from '../modules/local/harmonize_amrfinder'

workflow AMRSCAN {
    take:
    assemblies  // [ meta, fasta ]
    database

    main:
    AMRFINDERPLUS(assemblies, database)
    HARMONIZE_AMRFINDER(AMRFINDERPLUS.out.calls)

    emit:
    raw = AMRFINDERPLUS.out.calls
    harmonized = HARMONIZE_AMRFINDER.out.results
    versions = AMRFINDERPLUS.out.versions.mix(HARMONIZE_AMRFINDER.out.versions)
}
