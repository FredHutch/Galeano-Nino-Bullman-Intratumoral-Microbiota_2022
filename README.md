# Galeano-Nino-Bullman-Intratumoral-Microbiota-2022


[![DOI](https://zenodo.org/badge/530442339.svg)](https://zenodo.org/badge/latestdoi/530442339)


Analysis code used in Galeano Nino et al., Effect of the intratumoral microbiota on spatial and cellular heterogeneity in cancer. 2022

The code in this repository is organized to reflect the description in the Methods
section of Galeano Nino et al., Effect of the intratumoral microbiota on spatial and cellular heterogeneity in cancer. 2022.

## Environment and Reference Data

### Environment

All of the analysis code documented in this repository was run on the shared computing cluster
maintained at the Fred Hutchinson Cancer Research Center between May 2020 and August 2022.
The software dependencies used by these scripts are provided using the EasyBuild installation
maintained by the Fred Hutch Scientific Computing group.
Those software dependencies are loaded into the environment with the `ml` command (e.g. `ml CellRanger/6.1.1`).

### Reference Data

Prior to running the analysis scripts, reference databases were downloaded for PathSeq (December 2020)
and CellRanger (January 2022).
The location of those reference databases is provided to the analysis scripts using the environment variables `pathseqdb` and `cellrangerdb`.

# Overview of the Computational Pipeline for Bacteria-associated Spots/Cells Annotation

## Part 1: 10x Visium spatial transcriptomic data
## 10X Visium Scans have been uploaded to Zenodo: https://doi.org/10.5281/zenodo.7419806
   1. Identification of microbial reads within 10x Visium spatial transcriptomic data generated by 10x Space Ranger Count (`Visium_pipeline.sh`)
   2. Bioinformatic analysis of 10x Visium spatial transcriptomic data (`Visium.R`)
   3. summarize numbers of bacteria reads and UMIs in 10X Visium data (`validate_and_count.py`) The folder used as outputs from the previous steps should be provided as an argument to the `Pipeline_Visium.sh` script.
###   Output Data:
   - `CRC_16.visium.raw_matrix.genus.csv` and `OSCC_2.visium.raw_matrix.genus.csv` contain bacteria UMI counting matrix that can be used as metadata in visium data process
   - `CRC_16.visium.raw_matrix.validate.csv` and  `OSCC_2.visium.raw_matrix.validate.csv` contain validation data that can be used as the input of `validate_and_count.py`

## Part 2: 10x Single cell data (For cell culture samples and patient samples)
###   Input Data:
   - All of the input data for this analysis is provided in FASTQ format generated by the CellRanger `mkfastq` command
   - The folder containing those FASTQ files is set to the environment variable `raw_data_folder`
###   Preprocess:
   1. Identification of microbial reads within single cells GEX libraries (`patient_samples_GEX_pipeline.sh` and `cell_culture_samples_GEX_pipeline.sh`)
   2. INVADEseq bacterial 16S rRNA gene libraries (`patient_samples_16s_pipeline.sh` and `cell_culture_16s_pipeline.sh`). The variable `gex_bam_path` should be set to the output folder from the `patient_samples_GEX_pipeline.sh` and `cell_culture_samples_GEX_pipeline.sh` script.
   3. Combine and deduplication of microbial metadata from step 1 & 2 (`merge_metadata.py` and `metadata_dedup.py`). The folder used as outputs from the previous steps should be provided as an argument to the `merge_metadata.py` script.
###   Output Data:
   - `headneck_gex_16s_mix_dedup.csv` `HT_29_gex_16s_mix_dedup.csv` `HCT_116_csv_gex_16s_mix_dedup.csv` contain bacteria UMI counting matrix that can be used as Seurat object metadata in single cell process.

###   Processing of single cell data
   1. Seurat data processing, Harmony integration, SingleR annotation and copyKAT predication (`patient_samples_Seurat.r` and `cell_culture_Seurat.r`)
   2. Differentially expression analysis and GSEA (`DE.r`)
