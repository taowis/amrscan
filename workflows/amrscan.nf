include { AMRFINDERPLUS } from '../modules/local/amrfinderplus'
include { HARMONIZE_AMRFINDER } from '../modules/local/harmonize_amrfinder'

workflow AMRSCAN {
    take:
    assemblies
    database

    main:
    AMRFINDERPLUS(assemblies, database)
    HARMONIZE_AMRFINDER(AMRFINDERPLUS.out.calls)

    emit:
    raw = AMRFINDERPLUS.out.calls
    harmonized = HARMONIZE_AMRFINDER.out.results
}
