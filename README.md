## Docker Image
You can pull the ready-to-use environment directly from Docker Hub:
`docker pull jasminezonzini/zonzini-jasmine-exam-jan26:latest`

## External Dependencies and Data Availability
GitHub Repository: Due to size limitations (>50MB), the Ensembl human genome annotation file (Homo_sapiens.GRCh38.115.gtf.gz) and the large sparse matrix files are not included in this GitHub repository.

Docker Image (Recommended): The ready-to-use Docker image already includes all necessary dependencies, including the GTF annotation file and the processed RDS/Matrix files. This ensures the pipeline is fully reproducible without requiring manual downloads.

To run the analysis without any manual setup, please use the Docker environment:
`docker pull jasminezonzini/zonzini-jasmine-exam-jan26:latest`

## Data Loading & Performance Optimization

To grant reproducibility and overcome RAM limits, scRNA-seq analysis does not starts from the original matrix-data `.mtx` or `.h5`, but utilizes a Seurat object pre-processed. This final object (`merged_5pazienti_RNA.rds`) has been generated following the workflow described in the RMarkdown of this repository:

1. **Loading and Subsampling:** Each PDX has been reda in whole and subsapled to retain only 500 random cells, using the function `prepare_pdx_rna` (defined in `Script_functions.R`). This passage has been necessary to overcome limited RAM availability.
   
2. **Intergation:** the samples have been merged in a single Seurat object (`merged_5pazienti_RNA`).

3. **Matrix export:** To make the analysis reproducible and stantard conformed, a matrix has been extrapolated using the command
   ```R
   esporta_matrice_10x(merged_5pazienti_RNA)
   ```
in the RMarkdown, which refers to the function `esporta_matrice_10x` defined in the `Script_functions.R`
