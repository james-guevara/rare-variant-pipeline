include { FASTVEP_PICK; STANDALONE_LOFTEE } from '../modules/sites_annotation'

workflow SITES_ANNOTATION {
    take:
    sites
    runner
    benchmark
    main:
    FASTVEP_PICK(sites, runner, benchmark)
    STANDALONE_LOFTEE(FASTVEP_PICK.out.picked, runner)
    emit:
    picked = FASTVEP_PICK.out.picked
    loftee = STANDALONE_LOFTEE.out.annotated
}
