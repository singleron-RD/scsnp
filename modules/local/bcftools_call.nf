process BCFTOOLS_CALL {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bcftools==1.20"
    container "biocontainers/bcftools:1.20--h8b25389_0"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*vcf.gz")     , emit: vcf
    path  "versions.yml"                 , emit: versions

    script:
    def args = '--max-depth 1000000 --indels-cns'
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    bcftools \\
        mpileup \\
        --fasta-ref $fasta \\
        $args \\
        $bam \\
        | bcftools call --multiallelic-caller --variants-only --output-type z -o ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}