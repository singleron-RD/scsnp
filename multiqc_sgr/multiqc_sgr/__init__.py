def multiqc_sgr_config():
    from multiqc import config

    """ Set up MultiQC config defaults for this package """
    sgr_search_patterns = {
        "general_stats": {"fn": "*general_stats.json"},
        "scsnp/count": {"fn": "*count.json"},
        "scsnp/meta": {"fn": "*meta.json"},
    }
    config.update_dict(config.sp, sgr_search_patterns)
