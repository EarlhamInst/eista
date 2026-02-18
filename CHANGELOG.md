# nf-core/eista: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.2 - [2026-02-18]
- add a new module annotation_scvi for cell-type prediction using scvi-tools

## v2.0 - [2025-09-02]
- add support for 10x Xenium in-situ data

## v1.8 - [2025-07-02]
- add tertiary analysis: Cell-cell communication analysis

## v1.5 - [2025-02-10]
- add tertiary analysis: Cell type annotation
- add tertiary analysis: Differential expression analysis

## v1.2 - [2025-03-10]
- add secondary analysis: Spatial statistics analysis

## v1.0 - [2025-02-10]

Initial release of nf-core/eista, created with the [nf-core](https://nf-co.re/) template.
The pipeline including following analyses:
1. Primary analysis:
    - Cell segmentation 
    - partition transcripts
    - calculate cell metadata
    - sum signals
    - update vzg
2. Secondary analysis
    - Single-cell quality control
    - Cell filtering
    - Clustering analysis

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`
