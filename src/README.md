# Source Code

This directory contains the source code for the automated bacterial genome curation and analysis pipeline.

Scripts are organized by functionality and are intended to be executed as part of an integrated workflow rather than as standalone tools.

## Source Files

- `sra_extractor.py`  
  Handles automated retrieval of sequencing data from public repositories.

- `sra_extractor_metadata.py`  
  Extracts and validates associated metadata, including isolate identifiers and phenotypic information.

- `genome_checker.py`  
  Performs genome quality checks and validation steps.

- `amr_transformer.py`  
  Processes and standardizes antimicrobial resistance gene outputs.

- `phenotype_checker.py`  
  Validates and filters isolates based on phenotypic susceptibility data.

- `DataIntegrator.py`  
  Integrates outputs from multiple analysis stages into a unified dataset.

- `ml_preprocessor.py`  
  Prepares curated genomic features for downstream data-driven analyses.

- `MLBuilder.py`  
  Constructs structured datasets suitable for machine learning applications.

## Notes

- Code emphasizes automation, scalability, and reproducibility.
- Detailed methodological descriptions are provided in the accompanying manuscript.
