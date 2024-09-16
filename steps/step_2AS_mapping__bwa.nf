nextflow.enable.dsl=2

include { extractDsRef;parseMetadataFromFileName;executionMetadata;extractKey;taskMemory;getEmpty } from '../functions/common.nf'
include { getSingleInput;getReference;isIlluminaPaired;isCompatibleWithSeqType;isIonTorrent } from '../functions/parameters.nf'
include { stepInputs;getRisCd } from '../functions/common.nf'

def ex = executionMetadata()

def STEP = '2AS_mapping'
def METHOD = 'bwa' 
def ENTRYPOINT = "step_${STEP}__${METHOD}"

process bwa {
    container "biocontainers/bwa:v0.7.17_cv1"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 1.GB, task.attempt ) }
    when:
      referencePath && referencePath.exists() && !referencePath.empty() && (isCompatibleWithSeqType(reads, ['ion','illumina_paired'], task.process))
    input:
      tuple val(riscd_input), path(reads)
      tuple val(riscd_ref), val(reference), path(referencePath)
    output:
      path '*'
      tuple path("${base_ref}.sam"), val(reference), path(referencePath), emit: sam
      path '*.sh', hidden: true
    afterScript "echo '${stepInputs([riscd_input,riscd_ref], md, ex, STEP, METHOD, [reference:reference])}' > ${base_ref}_input.json"
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '{*.log,*.json}'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}.cfg" }
    script:
      (t1,t2) = (reads instanceof java.util.Collection) ? reads : [reads, null]
      md = parseMetadataFromFileName(t1.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      base_ref = "${base}_bwa_${reference}"
      if (isIlluminaPaired(reads)) {
        """
        bwa index -p ${reference} ${referencePath}
        bwa mem -t 8 ${reference} -1 ${t1} -2 ${t2} -o ${base_ref}.sam 2>> ${base_ref}.log
        """
      } else if (isIonTorrent(reads)) {
        """
        bwa index -p ${reference} ${referencePath} 
        bwa mem -t 8 ${reference} -U ${t1} -o ${base_ref}.sam 2>> ${base_ref}.log
        """      
      }      
}

process samtools {
    container "staphb/samtools:1.19"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 6.GB, task.attempt ) }
        input:
      tuple path(sam), val(reference), path(referencePath)
    output:
      tuple path("${base_ref}_sorted.bam"), val(reference), emit: bam
      tuple path("${base_ref}.fq"), val(reference), emit: fq
      path '*'
      path '*.sh', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '{*_sorted.bam*,*.vcf}'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}_samtools.cfg" }
    script:
      md = parseMetadataFromFileName(sam.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      base_ref = "${base}_bwa_${reference}"
      """
      samtools view -bS -o ${base_ref}.bam ${sam} 2>> ${base_ref}.log
      samtools sort -@ 8 ${base_ref}.bam ${base_ref}_sorted 2>> ${base_ref}.log
      samtools index ${base_ref}_sorted.bam 2>> ${base_ref}.log

    	samtools mpileup -uf ${referencePath} ${base_ref}_sorted.bam > ${base_ref}.bcf 2>> ${base_ref}.log
	    bcftools view -cg ${base_ref}.bcf > ${base_ref}.var.flt.vcf 2>> ${base_ref}.log
  	  vcfutils.pl vcf2fq ${base_ref}.var.flt.vcf > ${base_ref}.fq 2>> ${base_ref}.log
	    """
}

process seqio {
    container "python:3.10.1"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 250.MB, task.attempt ) }
    input:
      tuple path(fq), val(reference)
    output:
      tuple val(riscd), path("${base_ref}.fasta"), emit: consensus    
      path '*.sh', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*.fasta*'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}_seqio.cfg" }
    script:
      md = parseMetadataFromFileName(fq.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      base_ref = "${base}_bwa_${reference}"
      riscd = getRisCd(md, ex, STEP, METHOD)      
      """
      #!/usr/bin/env python3

      import os
      from Bio import SeqIO

      try:
          if os.path.getsize("${fq}") > 0:
              SeqIO.convert("${fq}", 'fastq', "${base_ref}.fasta", 'fasta')              
          else:
              print("WARNING no reads map on reference: ${reference}")
      except:
          print("Error: not found fastq file after mapping on: ${reference}")
	    """  
}

process samtools_depth {
    container 'quay.io/biocontainers/samtools:1.10--h9402c20_1'
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 200.MB, task.attempt ) }
        input:
      tuple path(bam), val(reference)
      val(method)
    output:
      tuple path("${base_ref}.coverage"), val(reference), emit: coverage
      path '*.sh', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}_samtools_depth.cfg" }
    script:
      md = parseMetadataFromFileName(bam.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      base_ref = "${base}_${method}_${reference}"
      """
      samtools depth -a ${bam} | awk '{ if (\$3!=0) c++;s+=\$3}{h++} END { if (c!=0) print s/c; else print 0;if (h!=0) print c/h; else print 0 }' > ${base_ref}.coverage
	    """
}

process coverage_minmax {
    container "staphb/samtools:1.19"
    containerOptions = "-v ${workflow.projectDir}/scripts/${ENTRYPOINT}:/scripts:ro"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 4.GB, task.attempt ) }
        input:
      tuple path(bam), val(reference)
      val(method)
    output:
      path '*.csv'
      path '*.sh', hidden: true
      tuple path("${base_ref}_samtools_depth.txt"), val(reference), emit: coverage_depth  
      tuple val(md.ds), val(reference), path("${base_ref}_import_coverage_minmax.csv"), emit: coverage_extra
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*.csv'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}_coverage_minmax.cfg" }
    script:
      md = parseMetadataFromFileName(bam.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      base_ref = "${base}_${method}_${reference}"
      """
      samtools view -F 4 -c ${bam} > samtools_view.txt
      samtools depth ${bam} > ${base_ref}_samtools_depth.txt
      /scripts/coverage_minmax.py ${md.cmp} ${md.ds} samtools_view.txt ${base_ref}_samtools_depth.txt ${base_ref}_import_coverage_minmax.csv
	    """
}

process coverage_check {
    container "python:3.10.1"
    containerOptions = "-v ${workflow.projectDir}/scripts/${ENTRYPOINT}:/scripts:ro"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 200.MB, task.attempt ) }
    input:
      tuple path(coverage), path(consensus), val(reference)
      val(context)
    output:
      tuple val(md.ds), val(reference), path("${base_ref}_import_coverage.csv"), emit: coverage_basic
      path '*'
      path '*.sh', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '{*.csv,*.check}'      
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}_coverage_check.cfg" }
    script:
      //TODO fix output folder
      md = parseMetadataFromFileName(consensus.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      base_ref = "${base}_${context}_${reference}"
      coverages= coverage.toRealPath().toFile().readLines()
      """
      /scripts/coverage.py ${coverage} ${consensus} ${reference} ${md.cmp} ${md.ds.substring(2)} ${base_ref}.check ${base_ref}_import_coverage.csv 
	    """
}

process coverage_check_group {
  container "ubuntu:20.04"
  tag "${md?.cmp}/${md?.ds}/${md?.dt}"
  memory { taskMemory( 200.MB, task.attempt ) }
  input:
    tuple val(key), path(files)
    val(method)
  output:
    path '*'
    path '*.sh', hidden: true
  publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*.csv'
  publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base}_${method}_coverage_check.cfg" }
  script:
    def coverage_file = (files instanceof java.util.Collection) ? files.flatten()[0] : files
    md = parseMetadataFromFileName(coverage_file.getName())
    base = "${md.ds}-${ex.dt}_${md.cmp}"
    """
    cat ${files} | sort -ur > ${base}_${method}_import_coverage_full.csv
    """  
}

process coverage_check_merge {
  container "ubuntu:20.04"
  tag "${md?.cmp}/${md?.ds}/${md?.dt}"
  memory { taskMemory( 200.MB, task.attempt ) }
  input:
    tuple val(key), val(reference), path(covMinMax), path(covBasic)
    val(method)
  output:
    tuple val(md.ds), path("${base_ref}_import_coverage_merged.csv"), emit: coverage_merged
    path '*'
    path '*.sh', hidden: true
  publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '*.csv'
  publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}_coverage_check_merge.cfg" }
  script:
    md = parseMetadataFromFileName(covMinMax.getName())
    base = "${md.ds}-${ex.dt}_${md.cmp}"
    base_ref = "${base}_${method}_${reference}"
    """
    paste -d, ${covMinMax} ${covBasic} | cut -d, -f1,2,3,4,5,8,9,10,11,12,13 > ${base_ref}_import_coverage_merged.csv
    """  
}

process coverage_plot {
    container "quay.io/biocontainers/matplotlib:3.1.2"
    containerOptions = "-v ${workflow.projectDir}/scripts/${ENTRYPOINT}:/scripts:ro"
    tag "${md?.cmp}/${md?.ds}/${md?.dt}"
    memory { taskMemory( 1.GB, task.attempt ) }
    input:
      tuple path(coverage_depth), val(reference)
    output:
      path '*'
      path '*.sh', hidden: true
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/result", pattern: '*.png'
    publishDir mode: 'rellink', "${params.outdir}/${md.anno}/${md.cmp}/${STEP}/${md.ds}-${ex.dt}_${METHOD}/meta", pattern: '.command.sh', saveAs: { "${base_ref}_coverage_plot.cfg" }
    script:
      md = parseMetadataFromFileName(coverage_depth.getName())
      base = "${md.ds}-${ex.dt}_${md.cmp}"
      base_ref = "${base}_bwa_${reference}"
      """
      /scripts/coverage_plot.py ${coverage_depth} ${base_ref}_coverage_plot.png
	    """
}

workflow step_2AS_mapping__bwa {
  take: 
    reads
    reference 
  main:
    bwa(reads, reference) //[ refId, refPath ]
    samtools(bwa.out.sam)
    consensus = seqio(samtools.out.fq).consensus

    coverage_minmax(samtools.out.bam, 'bwa')
    coverage_minmax.out.coverage_depth | coverage_plot
    
    coverage = samtools_depth(samtools.out.bam, 'bwa').coverage
    coverage.cross(consensus) { extractDsRef(it) }.map { 
        return [ it[0][0], it[1][1], it[0][1] ]
    }.set { coverageRefAndConsensus }
    coverageBasic = coverage_check(coverageRefAndConsensus, 'bwa').coverage_basic

    crossedChecks = coverage_minmax.out.coverage_extra.cross(coverageBasic) { it[0] + "-" + it[1] }
    .map { [ it[0][0], it[0][1], it[0][2], it[1][2] ] }
    coverage_check_group(coverage_check_merge(crossedChecks, 'bwa').coverage_merged | groupTuple, 'bwa')
  emit:
    consensus
  }

workflow {
    getSingleInput().cross(getReference('fa')) { extractKey(it) }
      .multiMap { 
          reads: it[0] // riscd, R[]
          refs:  it[1][1..3] // riscd, code, path
      }.set { input }  
    step_2AS_mapping__bwa(input.reads, input.refs)
}
