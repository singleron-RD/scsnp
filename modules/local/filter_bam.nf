process FILTER_BAM {
    tag "$meta.id"
    label 'process_single'

    conda 'bioconda::pysam==0.22.1'
    container "biocontainers/pysam:0.22.1--py38h15b938a_0"

    input:
    tuple val(meta), path(bam), path(match_barcode)
    val genes
    val panel

    output:
    tuple val(meta), path("*.filtered.bam"), emit: bam

    script:
    def genes_args = genes ? "--genes $genes": ""
    def panel_args = panel ? "--panel $panel": ""
    """
    filter_bam.py \
     --bam $bam \
     --match_barcode_file $match_barcode \
     --sample ${meta.id} \
     $genes_args \
     $panel_args
    """
}