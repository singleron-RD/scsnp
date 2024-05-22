process VCF_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "quay.io/singleron-rd/sccore:0.0.0"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.vcf.stats.json"), emit: stats_json
    tuple val(meta), path("*.GT.csv"), emit: gt_csv

    script:
    """
    vcf_stats.py \\
     --vcf $vcf \\
     --sample ${meta.id}
    """
}