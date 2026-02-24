## External Dependencies: Gene Annotations
The analysis of scATAC-seq data requires the Ensembl human genome annotation file. Due to size limitations on GitHub (>50MB), the following file is not included in this repository:
- **File:** `Homo_sapiens.GRCh38.115.gtf.gz`
- **Source:** [Ensembl FTP Release 115](https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/)

**Instructions for the user:**
To run the full pipeline, please download the GTF file and place it in the project root folder (or inside the `/data` directory as specified in the scripts) before building the Docker image or running the RMarkdown.
