# Project Summary
Analysis source code for the paper titled "Autologous cell therapy with CD133+ bone marrow-derived stem cells for Asherman Syndrome: a phase 1/2 trial".

## Source code
- `AS_Patients` : Analysis of endometrium biopsies of AS patients

  - `data_preprocesing_and_integration/`
    - `00_AS_QC.Rmd` : This script processes single-cell RNA-seq count matrices and performs quality control and low quality cell filtering
    - `01_merge_Biop.Rmd`: Load cell ranger count matrices, adds metadata, and performs first SCT integration
    - `02_integration1_Biop.Rmd`: First integration clustering, main cell types annotation (cell_type_AS) and filtering low quality and erythrocyte cells 
    - `03_integration2_Biop.Rmd`: Second integration
    - `04_annotation_Biop.Rmd` : Second integration clustering and main cell types annotation (cell_type_AS_v2)

  - `Differential abundance/` 
    - `differential abundance_testing_Biop.Rmd`: Performs differential abundance testing between cell types

  - `Secretory phase scoring/` 
    - `secretory_phase_score_Biop.Rmd`: This script calculates enrichment scores of gene signature characteristic of the secretory phase and WOI in the menstrual cycle

  - `Differential gene expression/`
        - `differential_expression_testing_Biop.Rmd`: This script performs differential expression between cell types and conditions

  - `Cell-to-cell communication/`
      - `cell-to-cell-communications_Biop.Rmd`: This script preforms cell-to-cell communication using CellChat 
      - `cellchat_figures_Biop.Rmd`: This script makes all necesary figures 

  - ` Lineage tracing/`
    - `analyze_variant/`
      - `explore_variants_final.r`: Used to explore and extract informative variants obtained through `MAESTER_slurm.sh`. Also generates several plots.
    - `get_HQ_barcodes.R`: Obtains HQcb (High Quality cell-barcodes) from a scRNA-seq object.
    - `parse_metrics_slurm.py`: Calculates several quality statistics from a `MAESTER_slurm.sh` execution.

    - `processing_pipeline/`
      - `inc_assemble_fastq.R`: Assembles FASTQ files for high-quality cells from raw data (used in `MAESTER_slurm.sh`).
      - `inc_calculate_coverage_cut.R`: Calculates the "good quality coverage" threshold in HQcb (used in `MAESTER_slurm.sh`).
      - `inc_create_af_dp_table.R`: Creates an allele frequency and depth table for all variants in all cells (used in `MAESTER_slurm.sh`).
      - `inc_parse_vcf_lofreq.py`: Corrects `lofreq` VCF format (used in `MAESTER_slurm.sh`).
      - `inc_split_bam.py`: Splits the MAESTER BAM file into one BAM per HQcb (used in `MAESTER_slurm.sh`).
      - `MAESTER_slurm.conf`: Configuration file with paths, software, and parameters (used with `MAESTER_slurm.sh`).
      - `MAESTER_slurm.sh`: Main script for processing and analyzing MAESTER data.

- `AS_Patients_CD133` : Analysis of CD133+ cells derived from AS patients

  - `data_preprocesing_and_integration/`
    - `00_AS_QC_CD133.Rmd` : This script processes single-cell RNA-seq count matrices and performs quality control and low quality cell filtering
    - `01_merge_CD133.Rmd`: Load cell ranger count matrices and metadata, and first SCT integration
    - `02_integration_CD133.Rmd`: SCT integration and clustering
    - `03_annotation_CD133.Rmd` : Cell types annotation
    
- `Organoids` : Analysis of AS patients derived organoids

  - `data_preprocesing_and_integration/`
    - `00_QC_Bender_Org.Rmd`: This script processes single-cell RNA-seq count matrices and performs quality control and low quality cell filtering
    - `01_merge_Bender_Org.Rmd`: Merge all samples and metadata in one seurat object. 
    - `02_integration1_Org.Rmd`: Preliminary SCT integration to remove low quality clusters
    - `03_Integration2_annotation_Org.Rmd`: Second SCT integration and annotation from in vivo epithelium 
    - `04_Integration3_Org.Rmd`: New integration to improve overall integration quality removing LH+8 sample
    - `05_Integration4_annotation_Org.rmd`: Performs the final integration pipeline using SCT approach with more restrictive quality filters 

  - `celltype_prediction/`
    - `06_organoids-celltype_prediction_model.Rmd` : Calculates similarity scores beetween organoids and in vivo counterparts

- `README.md`

## Enviroment specifications
- `env/` 
  - `renv.lock` : R environment for everything except MAESTER
  - `renv_maester.lock` : R environment for the MAESTER analysis
  - `requirements_maester.txt` : Python environment for the MAESTER analysis
    	
## Code description

| Analysis                                                                                                                        | Step                             | Figures                   | Source Code                                               | Notes                                                                                     |
| ------------------------------------------------------------------------------------------------------------------------------- | -------------------------------- | ------------------------- | --------------------------------------------------------- | ----------------------------------------------------------------------------------------- |
| Single-cell RNA sequencing of the AS endometrium pre and post treatment                                                         | Data collection                  | -                         | -                                                         | -                                                                                         |
|                                                                                                                                 | QC/Processing                    | Supp Fig. 2C              | `src/AS_Patients/data_preprocesing_and_integration`       | -                                                                                         |
|                                                                                                                                 | Integration and annotation       | Figure 4A, Supp Fig. 2E   | `src/AS_Patients/data_preprocesing_and_integration`       | Extended and cell type markers adapted from Santamaria X et al. Nat Commun 14:5890 (2023) |
|                                                                                                                                 | Secretory phase enrichment score | Figure 4B                 | `src/AS_Patients/secretory_phase_score`                   | -                                                                                         |
|                                                                                                                                 | Differential abundance           | Figure 4C, Supp Fig. 2D   | `src/AS_Patients/differential_abundance`                  | Extended from Santamaria X et al. Nat Commun 14:5890 (2023)                               |
|                                                                                                                                 | Differential expression          | Supp. Fig. 3              | `src/AS_Patients/differential_expression`                 | Extended from Santamaria X et al. Nat Commun 14:5890 (2023)                               |
|                                                                                                                                 | Cell-to-cell communication       | Figure 4D, Supp Fig. 4A–E | `src/AS_Patients/cell-to-cell_communication`              | CellChat settings adapted from Santamaria X et al. Nat Commun 14:5890 (2023).             |
| Lineage tracing differentiation of engrafted CD133+ cells                                                                       | Data collection                  | -                         | -                                                         | -                                                                                         |
|                                                                                                                                 | QC/Processing                    | -                         | `src/AS_Patients/maester/processing_pipeline`             | -                                                                                         |
|                                                                                                                                 | Profiling of mtDNA variants      | Figure 4E                 | `src/AS_Patients/maester/analyze_variants`                | -                                                                                         |
| Single-cell RNA sequencing of endometrial epithelial organoids (EEOs) derived from endometrial biopsies pre- and post-treatment | Data collection                  | -                         | -                                                         | -                                                                                         |
|                                                                                                                                 | QC/Processing                    | Figure 5B                 | `src/Organoids/data_preprocesing_and_integration`         | -                                                                                         |
|                                                                                                                                 | Integration and annotation       | Figure 5C                 | `src/Organoids/data_preprocesing_and_integration`         | Extended and cell type markers adapted from Santamaria X et al. Nat Commun 14:5890        |
|                                                                                                                                 | Differential abundance           | -                         | `src/Organoids/differential_abundance`                    | Extended from Santamaria X et al. Nat Commun 14:5890 (2023).                              |
| scRNA-seq profiling of the CD133+ product                                                                                       | Data collection                  | -                         | -                                                         | -                                                                                         |
|                                                                                                                                 | QC/Processing                    | -                         | `src/AS_Patients_CD133/data_preprocesing_and_integration` | -                                                                                         |
|                                                                                                                                 | Integration and annotation       | Supp Fig. 2B              | `src/AS_Patients_CD133/data_preprocesing_and_integration` | -                                                                                         |

