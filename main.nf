#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads = "$baseDir/00_fastq_raw/*.fq"
params.outdir = "$baseDir/results"
params.gentrome = "$baseDir/ref_fasta/gentrome.fa.gz"
params.decoys = "$baseDir/ref_fasta/decoys.txt"

log.info """\
 R N A S E Q - N F   P I P E L I N E
 ===================================
 gentrome     : ${params.gentrome}
 reads        : ${params.reads}
 outdir       : ${params.outdir}
 """

process FASTQC {
    tag "FASTQC on $sample_id"
    conda 'bioconda::fastqc=0.12.1'
    publishDir "$params.outdir/01_fastq_raw_FastQC/", mode:'copy'

    cpus 2

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}" 

    script:
    """
    mkdir -p fastqc_${sample_id}
    fastqc \\
        -t ${task.cpus} \\
        -o fastqc_${sample_id} \\
        ${reads}
    """
}

process INDEX {
    tag "$gentrome.simpleName"
    conda 'bioconda::salmon=1.10.2'

    cpus 28

    input:
    path gentrome 
    path decoys

    output:
    path 'index' 

    script:
    """
    salmon index \\
        --threads $task.cpus \\
        --gencode \\
        -t $gentrome \\
        -d $decoys \\
        -i index
    """
}

process TRIMGALORE {
    tag "TRIMGALORE on $sample_id"
    conda "bioconda::trim-galore=0.6.10"
    publishDir "$params.outdir/02_fastq_trimmed/", mode:'copy'

    cpus 12

    input:
    tuple val(sample_id), path(reads) 

    output:
    tuple val(sample_id), path("trimgalore_${sample_id}/${sample_id}_trimmed.fq")

    script:
    """
    trim_galore \\
            --fastqc \\
            --cores 8 \\
            --output_dir trimgalore_${sample_id} \\
            ${reads}
    """
}

process QUANT {
    tag "QUANT on $sample_id"
    conda 'bioconda::salmon=1.10.2'
    publishDir "$params.outdir/03_salmon_quant/", mode:'copy'

    cpus 28

    input:
    path index 
    tuple val(sample_id), path(reads) 

    output:
    path "quant_${sample_id}"

    script:
    """
    salmon quant \\
        --threads $task.cpus \\
        --libType A \\
        -i $index \\
        --validateMappings \\
        --gcBias \\
        --numGibbsSamples 20 \\
        -r ${reads} \\
        -o quant_${sample_id}
    """
}

workflow {
    reads_ch = Channel
                .fromPath(params.reads, checkIfExists: true)
                .map{tuple(it.baseName, it)}

    FASTQC(reads_ch)
    INDEX(params.gentrome, params.decoys)
    TRIMGALORE(reads_ch)
    QUANT(INDEX.out, TRIMGALORE.out)
}

/* 
 * completion handler
 */
workflow.onComplete {
	log.info ( workflow.success ? "\nSuccessfully completed!\n" : "\nOops .. something went wrong\n" )
}
