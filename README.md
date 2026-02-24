## External Dependencies: Gene Annotations
The analysis of scATAC-seq data requires the Ensembl human genome annotation file. Due to size limitations on GitHub (>50MB), the following file is not included in this repository:
- **File:** `Homo_sapiens.GRCh38.115.gtf.gz`
- **Source:** [Ensembl FTP Release 115](https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/)

**Instructions for the user:**
To run the full pipeline, please download the GTF file and place it in the project root folder (or inside the `/data` directory as specified in the scripts) before building the Docker image or running the RMarkdown.

## Data Loading & Performance Optimization

To grant reproducibility and overcome RAM limits, scRNA-seq analysis does not starts from the original matrix-data `.mtx` or `.h5`, but utilizes a Seurat object pre-processed. This final object (`merged_5pazienti_RNA.rds`) has been generated following the workflow described in the RMarkdown of this repository:

1. **Loading and Subsampling:** Each PDX has been reda in whole and subsapled to retain only 500 random cells, using the function `prepare_pdx_rna` (defined in `Script_functions.R`). This passage has been necessary to overcome limited RAM availability.
   
2. **Intergation:** the samples have been merged in a single Seurat object (`merged_5pazienti_RNA`).

3. **Matrix export:** To make the analysis reproducible and stantard conformed, a matrix has been extrapolated using
4. ```R
   esporta_matrice_10x(merged_5pazienti_RNA, path = "../data/RNA/")
   ```
   again defined in the same `prepare_pdx_rna` function in `Script_functions.R`
