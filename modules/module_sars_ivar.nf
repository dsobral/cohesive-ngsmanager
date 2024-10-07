nextflow.enable.dsl=2

include { step_2AS_mapping__ivar } from '../steps/step_2AS_mapping__ivar'
include { step_4TY_lineage__pangolin } from '../steps/step_4TY_lineage__pangolin'
include { extractKey } from '../functions/common.nf'
include { getSingleInput } from '../functions/parameters.nf'

def referenceCode = 'NC_045512.2'
def referencePath = "${params.assets_dir}/module_covid_emergency/NC_045512.fasta"
def referenceRiscd = '220308-020220308005121273-2AS_import-external'

workflow module_sars_ivar {
    take: 
        trimmed
    main:
        trimmed.multiMap {
            trimmed: it
            reference: [ referenceRiscd, referenceCode, file(referencePath) ]
        }.set { trAndRef }
        consensus = step_2AS_mapping__ivar(trAndRef.trimmed, trAndRef.reference).consensus
        step_4TY_lineage__pangolin(consensus)
    }

workflow {
    module_sars_ivar(getSingleInput())
}