FROM bioconductor/bioconductor_docker:RELEASE_3_18

RUN apt-get update && apt-get install -y --no-install-recommends \
    libhdf5-dev \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    tabix \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('Seurat', 'hdf5r', 'tidyverse', 'patchwork', 'ggplot2', 'devtools', 'roxygen2', 'knitr', 'rmarkdown', 'rtracklayer'))"

RUN R -e "BiocManager::install(c('Signac', 'rtracklayer', 'GenomicRanges', 'GenomeInfoDb'))"

WORKDIR /home/rstudio/app


COPY PACKAGE /home/rstudio/app/ZonziniJasmineExamJan26 

RUN R -e "devtools::document('/home/rstudio/app/ZonziniJasmineExamJan26')" 
RUN R -e "devtools::install('/home/rstudio/app/ZonziniJasmineExamJan26', build_vignettes = TRUE)"