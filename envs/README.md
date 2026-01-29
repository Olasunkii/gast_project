# Environment Definitions

This directory contains Conda environment definition files used to manage software dependencies across different stages of the pipeline.

Each environment is task-specific to reduce dependency conflicts and improve reproducibility.

## Environment Files

- `environment_amr.yaml`  
  Environment for antimicrobial resistance gene detection and annotation tools.

- `environment_bakta.yaml`  
  Environment for genome annotation using Bakta.

- `environment_checkm.yaml`  
  Environment for genome completeness and contamination assessment using CheckM.

- `environment_qc.yaml`  
  Environment for raw read quality control and preprocessing.

- `environment_quast.yaml`  
  Environment for genome assembly quality assessment using QUAST.

- `environment_wgs_depth.yaml`  
  Environment for whole-genome sequencing depth and coverage analysis.

- `environment_python.yaml`  
  General-purpose Python environment used for data processing, integration, and visualization.

## Notes

- Environments were created to ensure reproducibility and tool version consistency.
- Files are provided for reference and reuse.
