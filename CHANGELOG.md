# singleron-RD/scsnp: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 1.1.0 - [2024-06-05]

### `Added`
- Make `--genes` and `--fasta` parameter required.
- Update freebayes_args and bcftools_filter_args to filter more variants.

## 1.2.0 - [2024-06-25]

### `Added`
- Add variant calling software `bcftools`

### `Changed`
- Change the default variant calling tools from `freebayes` to `bcftools` to improve speed.
