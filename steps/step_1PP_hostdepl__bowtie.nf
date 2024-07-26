nextflow.enable.dsl=2

include { parseMetadataFromFileName; executionMetadata;taskMemory } from '../functions/common'
include { param;getSingleInput;getHostUnkeyed;isIlluminaPaired;isCompatibleWithSeqType;isIonTorrent } from '../functions/parameters.nf'
include { stepInputs;getRisCd } from '../functions/common.nf'


def refDir = param('hosts_dir')

def ex = executionMetadata()

def STEP = '1PP_hostdepl'
def METHOD = 'bowtie' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

process bowtie2 {
    container "biocontainers/bowtie2:2.2.9"
    containerOptions = "-v ${refDir}:${refDir}:ro"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 5.GB, task.attempt ) }
    when:
      isCompatibleWithSeqType(reads, ['ion','illumina_paired'], task.process) 
    input:
      tuple val(riscd_input), path(reads), val(host)
    output:
      path '*'
      tuple path("${outname}.sam"), val(host), val(type), emit: sam
      path '*.sh', hidden: true
    afterScript "echo '${stepInputs(riscd_input, md, ex, STEP, METHOD, [reference:host])}' > ${base}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: "{*.log,*.json}"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${outname}.cfg" }
    script:
      (t1,t2) = (reads instanceof java.util.Collection) ? reads : [reads, null]
      md = parseMetadataFromFileName(t1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      outname = "${base}_vdhost_${host?.replace("_", "")}"
      type = isIlluminaPaired(reads) ? 'paired' : 'single' 
      if (isIlluminaPaired(reads)) {
        """
        bowtie2 -p 8 --very-fast -x ${refDir}/${host} -1 ${t1} -2 ${t2} -S ${outname}.sam 2>> ${outname}.log
        """
      } else if (isIonTorrent(reads)) {
        """
        bowtie2 -p 8 --very-fast -x ${refDir}/${host} -U ${t1} -S ${outname}.sam 2>> ${outname}.log
        """      
      }  
}

process samtools {
    container 'quay.io/biocontainers/samtools:0.1.19--hf89b575_7'
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 18.GB, task.attempt ) }
    input:
      tuple path(sam), val(host), val(type)
    output:
      tuple val(riscd), path('*.fastq.gz'), val(host), emit: depleted
      path '*', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: "${base}_R*.fastq.gz"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: "${base}*.log"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${outname}_samtools.cfg" }
    script:
      md = parseMetadataFromFileName(sam.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      tmpName = "${base}_tmp"
      hostWithoutUnderscore = host?.replace("_", "")
      outname = "${base}_vdhost_${hostWithoutUnderscore}"
      vdhostR1 = "${base}_R1_vdhost_${hostWithoutUnderscore}.fastq.gz"
      vdhostR2 = "${base}_R2_vdhost_${hostWithoutUnderscore}.fastq.gz"
      riscd = getRisCd(md, ex, STEP, METHOD)      
      if (type == 'paired') {
        """
		    samtools view -bS ${sam} -@ 16 > ${tmpName}.bam 2>> ${outname}_samtools.log
		    samtools view -b -f 12 -F 256 ${tmpName}.bam -@ 16 > ${tmpName}_um.bam 2>> ${outname}_samtools.log
		    samtools sort -n ${tmpName}_um.bam ${tmpName}.sorted -@ 16 2>> ${outname}_samtools.log
		    samtools bam2fq ${tmpName}.sorted.bam > ${tmpName}.fastq 2>> ${outname}_samtools.log
        grep '^@.*/1\$' -A 3 ${tmpName}.fastq | grep -v -- "^--\$" | gzip > ${vdhostR1}
		    grep '^@.*/2\$' -A 3 ${tmpName}.fastq | grep -v -- "^--\$" | gzip > ${vdhostR2}   
        rm ${tmpName}*.*
	    """
}     else if (type == 'single') {
        """
		    samtools view -bS ${sam} -@ 16 > ${tmpName}.bam 2>> ${outname}_samtools.log
		    samtools view -b -f 4 -F 256 ${tmpName}.bam -@ 16 > ${tmpName}_um.bam 2>> ${outname}_samtools.log
		    samtools sort -n ${tmpName}_um.bam ${tmpName}.sorted -@ 16 2>> ${outname}_samtools.log
		    samtools bam2fq ${tmpName}.sorted.bam > ${tmpName}.fastq 2>> ${outname}_samtools.log
        gzip -c ${tmpName}.fastq > ${vdhostR1} 
        rm ${tmpName}*.*
        """
      }   
}

workflow step_1PP_hostdepl__bowtie {
    take: 
      trimmedAndHost
    main:
      bowtie2(trimmedAndHost)
      samtools(bowtie2.out.sam)
    emit:
      samtools.out.depleted      
}

workflow {    
    input = getSingleInput()
      .combine(getHostUnkeyed())
    step_1PP_hostdepl__bowtie(input)    
}
