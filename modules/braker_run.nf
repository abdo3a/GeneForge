nextflow.enable.dsl=2

process BRAKER_RUN {
    tag "${meta.id}"
    label 'process_medium'

    publishDir "${params.outdir}/braker", mode: 'copy', pattern: '*.braker_error.log'

    container 'docker://teambraker/braker3:latest'

    input:
    tuple val(meta),
          path(genome_masked),
          path(protein_evidence),
          path(bam),
          path(plus_bam),
          path(minus_bam),
          val(gc_probability)

    output:
    tuple val(meta), path("braker.gtf"), emit: braker_gtf
    tuple val(meta), path("braker"), emit: braker_dir
    path "${meta.id}_braker_error.log", emit: error_log

    script:
    def prefix = meta.id
    def species = meta.species.replaceAll(/\s+/, '_')
    def busco_db = meta.busco_db
    
    // BAM placeholders to recognize as missing files
    def bam_placeholders = ['NO_BAM_FILE.bam', '', null]
    def plus_bam_placeholders = ['NO_PLUS_BAM_FILE.bam', '', null]
    def minus_bam_placeholders = ['NO_MINUS_BAM_FILE.bam', '', null]

    // Check BAM availability by excluding placeholders
    def bam_available = (bam != null && !bam_placeholders.contains(bam.getName()))
    def dual_bam_available = (plus_bam != null && minus_bam != null &&
                              !plus_bam_placeholders.contains(plus_bam.getName()) &&
                              !minus_bam_placeholders.contains(minus_bam.getName()))

   
    def bam_flag = ''
    def stranded_flag = ''

    if (meta.stranded in ['forward', 'reverse'] && dual_bam_available) {
        bam_flag = "--bam=${plus_bam},${minus_bam}"
        stranded_flag = "--stranded=+,-"
    } else if (bam_available) {
        bam_flag = "--bam=${bam}"
        stranded_flag = ''
    } else {
        // No BAM files, so don't pass bam-related flags
        bam_flag = ''
        stranded_flag = ''
    }

    def gc_flag = (gc_probability && gc_probability != '') ? "--gc_probability=${gc_probability}" : ""

    """
    #!/bin/bash
    set -euo pipefail

    echo '[INFO] BRAKER inputs:' > ${prefix}_braker_error.log
    echo '  meta.id=${meta.id}' >> ${prefix}_braker_error.log
    echo '  meta.use_dual_bams=${meta.use_dual_bams}' >> ${prefix}_braker_error.log
    echo '  meta.stranded=${meta.stranded}' >> ${prefix}_braker_error.log
    echo '  bam=${bam}' >> ${prefix}_braker_error.log
    echo '  plus_bam=${plus_bam}' >> ${prefix}_braker_error.log
    echo '  minus_bam=${minus_bam}' >> ${prefix}_braker_error.log
    echo '  protein_evidence=${protein_evidence}' >> ${prefix}_braker_error.log
    echo '  gc_probability=${gc_probability}' >> ${prefix}_braker_error.log
    echo '  bam_flag=${bam_flag}' >> ${prefix}_braker_error.log
    echo '  stranded_flag=${stranded_flag}' >> ${prefix}_braker_error.log
    echo '  gc_flag=${gc_flag}' >> ${prefix}_braker_error.log
     
    braker.pl \\
      --genome=${genome_masked} \\
      ${bam_flag} \\
      ${stranded_flag} \\
      --prot_seq=${protein_evidence} \\
      --species="${species}" \\
      ${gc_flag} \\
      --threads ${task.cpus} \\
      --busco_lineage=${busco_db} \\
      --useexisting 2>> ${prefix}_braker_error.log

    mv braker/braker.gtf .
  
    """
}

