# --- LIBRERIE ---
library(Seurat)
library(Signac)
library(rtracklayer)
library(GenomicRanges)
library(tidyverse)

#' Caricamento e Subsampling dei dati RNA
#'
#' Legge un file h5 di 10X e seleziona un numero prefissato di cellule casuali.
#' @param file_path Percorso del file .h5
#' @param project_name Nome del progetto per l'oggetto Seurat
#' @param n_cells Numero di cellule da mantenere (default 500)
#' @return Un oggetto Seurat con dati subsampled
#' @export
prepare_pdx_rna <- function(file_path, project_name, n_cells = 500) {
  input_data <- Seurat::Read10X_h5(file_path)
  if (is.list(input_data)) {
    temp_mat <- if ("Gene Expression" %in% names(input_data)) input_data[["Gene Expression"]] else input_data[[1]]
  } else {
    temp_mat <- input_data
  }
  set.seed(123)
  if (ncol(temp_mat) > n_cells) {
    cells_to_keep <- sample(colnames(temp_mat), n_cells)
    temp_mat <- temp_mat[, cells_to_keep]
  }
  obj <- Seurat::CreateSeuratObject(counts = temp_mat, project = project_name)
  return(obj)
}

#' Esportazione Matrice in formato 10X
#'
#' Salva la matrice RNA in una cartella compatibile con Read10X.
#' @param seurat_obj Oggetto Seurat
#' @param output_path Percorso della cartella di output
#' @export
esporta_matrice_10x <- function(seurat_obj, output_path = "../data/RNA/MATRICE_SPARSA") {
  if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  Seurat::Write10X(seurat_obj[["RNA"]], data.dir = output_path)
}

#' Rimozione Geni Murini (Stroma)
#'
#' Filtra l'oggetto mantenendo solo i geni espressi con nomenclatura umana (tutto maiuscolo).
#' @param seurat_obj Oggetto Seurat PDX
#' @return Oggetto Seurat senza geni di topo
#' @export
rimuovi_geni_topo <- function(seurat_obj) {
  geni_totali <- rownames(seurat_obj)
  geni_umani <- geni_totali[grep("^[A-Z0-9]+$", geni_totali)]
  obj_filtrato <- subset(seurat_obj, features = geni_umani)
  return(obj_filtrato)
}

#' Filtrazione Ribosomiale basata su GTF
#'
#' Filtra i geni mantenendo solo le proteine ribosomiali identificate via GTF.
#' @param seurat_obj Oggetto Seurat
#' @param gtf_path Percorso del file GTF Ensembl
#' @return Oggetto Seurat filtrato per geni RP
#' @export
filtra_ribosomiali <- function(seurat_obj, gtf_path) {
  gtf <- rtracklayer::import(gtf_path)
  gtf_rp <- gtf[grepl("^RP", gtf$gene_name), ]
  geni_comuni <- intersect(unique(gtf_rp$gene_name), rownames(seurat_obj))
  obj_filtrato <- subset(seurat_obj, features = geni_comuni)
  info_per_granges <- as.data.frame(gtf_rp[gtf_rp$gene_name %in% geni_comuni, ])
  info_per_granges <- info_per_granges[!duplicated(info_per_granges$gene_name), ]
  rownames(info_per_granges) <- info_per_granges$gene_name
  info_per_granges <- info_per_granges[rownames(obj_filtrato), ]
  obj_filtrato[["RNA"]] <- AddMetaData(obj_filtrato[["RNA"]], metadata = info_per_granges)
  return(obj_filtrato)
}

#' Analisi Standard RNA per geni RP
#'
#' Esegue normalizzazione, PCA e UMAP.
#' @param seurat_obj Oggetto Seurat filtrato
#' @return Oggetto Seurat processato
#' @export
analizza_rna_ribo <- function(seurat_obj) {
  obj <- Seurat::NormalizeData(seurat_obj, verbose = FALSE)
  obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", verbose = FALSE)
  obj <- Seurat::ScaleData(obj, verbose = FALSE)
  obj <- Seurat::RunPCA(obj, verbose = FALSE)
  obj <- Seurat::RunUMAP(obj, dims = 1:10, verbose = FALSE)
  return(obj)
}

#' Creazione GRanges da Metadati Seurat
#'
#' @param seurat_filtrato Oggetto Seurat filtrato con filtra_ribosomiali
#' @return Oggetto GRanges
#' @export
crea_granges_ribo <- function(seurat_filtrato) {
  gtf_filtrato <- seurat_filtrato[["RNA"]][[]]
  gr <- makeGRangesFromDataFrame(gtf_filtrato, keep.extra.columns = TRUE,
                                 seqnames.field = "seqnames", start.field = "start", end.field = "end")
  return(gr)
}

#' Setup File ATAC
#'
#' @param input_path Percorso file .tsv.gz
#' @return Percorso del file .bgz.gz indicizzato
#' @export
setup_atac_files <- function(input_path) {
  file_decompressed <- sub("\\.gz$", "", input_path)
  file_final <- paste0(file_decompressed, ".bgz.gz")
  file_index <- paste0(file_final, ".tbi")
  if (!file.exists(file_decompressed) && !file.exists(file_final)) system(paste("gunzip -k", input_path))
  if (!file.exists(file_final)) {
    system(paste("bgzip -c", file_decompressed, ">", file_final))
    file.remove(file_decompressed)
  }
  if (!file.exists(file_index)) system(paste("tabix -p bed", file_final))
  return(file_final)
}

#' Caricamento Oggetto ATAC
#'
#' @param path Percorso del file .rds
#' @return Oggetto Seurat ATAC
#' @export
carica_oggetto_atac <- function(path = "/app/data/ATAC/atac_merged_con_clusters.rds") {
  if (!file.exists(path)) stop(paste("Errore:", path))
  return(readRDS(path))
}

#' Analisi ATAC (TF-IDF, SVD, UMAP)
#'
#' @param atac_obj Oggetto ATAC
#' @param dims Dimensioni LSI (default 2:15)
#' @return Oggetto ATAC processato
#' @export
processa_pdx_atac <- function(atac_obj, dims = 2:15) {
  atac_obj <- Signac::RunTFIDF(atac_obj, verbose = FALSE)
  atac_obj <- Signac::FindTopFeatures(atac_obj, min.cutoff = 'q0', verbose = FALSE)
  atac_obj <- Signac::RunSVD(atac_obj, verbose = FALSE)
  atac_obj <- Seurat::RunUMAP(atac_obj, reduction = "lsi", dims = dims, verbose = FALSE)
  return(atac_obj)
}

#' Validazione Separazione Pazienti
#'
#' @param atac_obj Oggetto ATAC con clustering
#' @param res_col Colonna dei cluster da validare
#' @return Table di contingenza
#' @export
valida_separazione_pdx <- function(atac_obj, res_col = "ATAC_snn_res.0.5") {
  return(table(atac_obj[[res_col, drop = TRUE]], atac_obj$orig.ident))
}