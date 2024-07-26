nextflow.enable.dsl=2

include { parseMetadataFromFileName;executionMetadata;taskMemory;taskTime } from '../functions/common.nf'
include { getInput;isIonTorrent;isIlluminaPaired;isCompatibleWithSeqType } from '../functions/parameters.nf'
include { stepInputs;getRisCd;extractKey } from '../functions/common.nf'

def ex = executionMetadata()

def STEP = '1PP_trimming'
def METHOD = 'fastp' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

process fastp {
    container "nanozoo/fastp:0.23.1--9f2e255"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 3.GB, task.attempt ) }
    time { taskTime( 10.m, task.attempt ) }
    when:
      isCompatibleWithSeqType(reads, ['ion','illumina_paired'], task.process)    
    input:
      tuple val(riscd_input), path(reads)
    output:
      path '*'
      tuple val(riscd), path('*.fastq.gz'), emit: trimmed
      path '{*.sh,*.log}', hidden: true
       afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, null)}' > ${base}_input.json"
       publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*.gz'
       publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '{*.json,*.html}'
       publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.log', saveAs: { "${base}.log" }
       publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}.cfg" }
    script:
      (r1,r2) = (reads instanceof java.util.Collection) ? reads : [reads, null]
      md = parseMetadataFromFileName(r1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}_${METHOD}"
      riscd = getRisCd(md, ex, STEP, METHOD)      
      if (isIlluminaPaired(reads)) { 
        """
          fastp --in1 ${r1} --out1 ${base}_R1.fastq.gz --in2 ${r2} --out2 ${base}_R2.fastq.gz \
          --unpaired1 ${base}_unpaired.fastq.gz --unpaired2 ${base}_unpaired.fastq.gz \
          --json ${base}_summary.json --html ${base}_summary.html --thread 8        """
      } else if (isIonTorrent(reads)) {
        """
          fastp --in1 ${r1} --out1 ${base}_R1.fastq.gz  \
          --json ${base}_summary.json --html ${base}_summary.html --thread 8        
        """        
      }
}

process sample_reads_check {
    container "quay.io/biocontainers/biopython:1.78"
    containerOptions = "-v ${workflow.projectDir}/scripts/${ENTRYPOINT}:/scripts:ro"
    memory { taskMemory( 14.GB, task.attempt ) }
    time { taskTime( 15.m, task.attempt ) }    
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    input:
      tuple val(_), path(reads)
      tuple val(_), path(trimmed)
    output:
      path '*'
      path '*.sh', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/qc/meta", pattern: '*SRC_raw.*[kg]'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/qc/result", pattern: '*SRC_raw.csv'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/qc/meta", pattern: '*SRC_treads.*[kg]'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/qc/result", pattern: '{*SRC_treads.csv,*_readsCheck.csv}'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/qc/meta", pattern: '.command.sh', saveAs: { "${base}_readscheck.cfg" }
    script:
      (r1,r2) = (reads instanceof java.util.Collection) ? reads : [reads, null]
      (t1,t2,u) = (trimmed instanceof java.util.Collection) ? trimmed : [trimmed, null, null]
      md = parseMetadataFromFileName(r1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      if (isIlluminaPaired(reads)) { 
        """
        /scripts/SampleReadsCheck.py -n $base -R1 $r1 -R2 $r2 -T1 $t1 -T2 $t2 -U $u > ${base}_SRC_raw.log
        cat ${base}_SRC_raw.log >	${base}_SRC_treads.log;
        """
      } else if (isIonTorrent(reads)) {
        """
        /scripts/SampleReadsCheck_ionTorrent.py -n $base -R1 $r1 -T1 $t1 > ${base}_SRC_raw.log
        """
      }      
}

process fastqc {
    container 'biocontainers/fastqc:v0.11.5_cv4'
    memory { taskMemory( 400.MB, task.attempt ) }
    time { taskTime( 5.m, task.attempt ) }    
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    input:
      tuple val(riscd_input), path(reads)
    output:
      path '*'
      path '*.sh', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/qc/result", pattern: '{*.zip,*.html}', saveAs: { filename -> filename.replaceFirst("-DT\\d+_", "-${ex.dt}_") }
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/qc/meta", pattern: '*.log'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/qc/meta", pattern: '.command.sh', saveAs: { "${base}_fastqc.cfg" }
    script:
      md = parseMetadataFromFileName(reads[0].getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      """
      fastqc $reads > "${base}_fastqc.log" 2>&1
      """
}

workflow step_1PP_trimming__fastp {
    take: 
      rawreads
    main:
      trimmed = fastp(rawreads).trimmed

      fastqc(trimmed)
      readsCheckInput = rawreads.cross(trimmed) { extractKey(it) }.multiMap { 
        rawreads: it[0]
        trimmed: it[1]
      }      
      sample_reads_check(readsCheckInput.rawreads, readsCheckInput.trimmed)

    emit:
      trimmed      
}

workflow {
    step_1PP_trimming__fastp(getInput())
}
