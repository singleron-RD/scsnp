process VCF_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "quay.io/singleron-rd/sccore:0.0.0"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.count.json"), emit: count_json
    tuple val(meta), path("*.GT.csv"), emit: gt_csv
    tuple val(meta), path("*.meta.json"), emit: meta_json

    script:
    """
    vcf_stats.py \\
     --vcf $vcf \\
     --sample ${meta.id}
    """
}