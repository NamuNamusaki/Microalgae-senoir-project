library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(GSVA)

# =========================
# LOAD ANNOTATION
# =========================
load_annotation <- function(file) {
  
  anno <- fread(file, header=FALSE, data.table=FALSE)
  
  colnames(anno)[1]  <- "gene_id"
  colnames(anno)[12] <- "KEGG_ko"
  
  anno %>%
    mutate(KEGG_ko = str_replace_all(KEGG_ko, "ko:", "")) %>%
    separate_rows(KEGG_ko, sep=",")
}

# =========================
# LOAD TPM
# =========================
load_TPM <- function(path) {
  
  files <- list.files(path, pattern="rsem.genes.results", full.names=TRUE)
  
  if (length(files) == 0) {
    stop(" No TPM files found")
  }
  
  dat <- fread(files[1])
  TPM_matrix <- data.frame(Gene = dat$gene_id)
  
  for (f in files){
    tmp <- fread(f)
    sample_name <- sub("_rsem.genes.results","",basename(f))
    TPM_matrix[[sample_name]] <- tmp$TPM
  }
  
  return(TPM_matrix)
}

# =========================
# PREPROCESS
# =========================
prepare_expr <- function(TPM_matrix){
  
  TPM_log <- log2(TPM_matrix[,-1] + 1)
  rownames(TPM_log) <- TPM_matrix$Gene
  
  expr <- as.matrix(TPM_log)
  
  # filter low expression
  expr <- expr[rowMeans(expr) > 1, ]
  
  return(expr)
}

# =========================
# RUN MODEL (CORE)
# =========================
run_model <- function(expr, feature_genes, label_genes, model_name, output_path){
  
  feature_genes <- intersect(rownames(expr), feature_genes)
  label_genes   <- intersect(rownames(expr), label_genes)
  
  if (length(feature_genes) < 2) {
    stop(paste(" Not enough feature genes for", model_name))
  }
  
  if (length(label_genes) < 2) {
    stop(paste(" Not enough label genes for", model_name))
  }
  
  # ---------------------
  # GSVA (feature)
  # ---------------------
  gsva_feature <- gsva(expr,
                       list(feature = feature_genes),
                       method="gsva",
                       kcdf="Gaussian")
  
  feature_score <- gsva_feature[1,]
  
  # ---------------------
  # GSVA (label)
  # ---------------------
  gsva_label <- gsva(expr,
                     list(label = label_genes),
                     method="gsva",
                     kcdf="Gaussian")
  
  label_score <- gsva_label[1,]
  z <- as.numeric(scale(label_score))
  
  label <- ifelse(z >= 1, "High",
                  ifelse(z <= -1, "Low", "Medium"))
  
  # ---------------------
  # OUTPUT
  # ---------------------
  # ตั้งชื่อ column ตาม model
  score_name <- paste0(model_name, "_score")

  result <- data.frame(
    Sample = colnames(expr),
    #Feature_score = round(feature_score, 4),
    score = round(label_score, 4),
    Zscore = round(z, 4),
    Label = label
  )
  # เปลี่ยนชื่อ column Score → PUFAs_score หรือ Astaxanthin_score
  colnames(result)[2] <- score_name
  
  write.csv(result,
            file.path(output_path, paste0(model_name, "_model.csv")),
            row.names=FALSE)
  
  cat("", model_name, "DONE\n")
}

# =========================
# MAIN PIPELINE
# =========================
run_pipeline <- function(CONFIG){

  # ======================
  # 0.CREATE SPECIES FOLDER
  # ======================
  species_output <- file.path(CONFIG$output_gsva, CONFIG$species)

  dir.create(species_output, showWarnings = FALSE, recursive = TRUE)
  
  # ---------------------
  # 1. LOAD ANNOTATION
  # ---------------------
  anno <- load_annotation(CONFIG$annotation_file)
  
  # ---------------------
  # 2. MAP GENES
  # ---------------------
  fatty_genes <- anno %>%
    filter(KEGG_ko %in% CONFIG$fatty_ko) %>%
    pull(gene_id) %>%
    unique()
  
  asta_genes <- anno %>%
    filter(KEGG_ko %in% CONFIG$asta_ko_all) %>%
    pull(gene_id) %>%
    unique()
  
  label_genes <- anno %>%
    filter(KEGG_ko %in% CONFIG$asta_ko_label) %>%
    pull(gene_id) %>%
    unique()
  
  # ---------------------
  # 3. LOAD TPM
  # ---------------------
  TPM <- load_TPM(CONFIG$tpm_dir)
  
  # ---------------------
  # 4. PREPARE EXPRESSION
  # ---------------------
  expr <- prepare_expr(TPM)
  
  # ---------------------
  # 5. RUN MODELS
  # ---------------------
  
  # 🟢 PUFAs model → label = ทุก KO
  run_model(expr,
            feature_genes = fatty_genes,
            label_genes   = fatty_genes,
            model_name    = "PUFAs",
            output_path   = species_output)
  
  # 🟠 Astaxanthin model → label = crtZ + crtW
  run_model(expr,
            feature_genes = asta_genes,
            label_genes   = label_genes,
            model_name    = "Astaxanthin",
            output_path   = species_output)
  
  cat(" PIPELINE COMPLETED\n")
}