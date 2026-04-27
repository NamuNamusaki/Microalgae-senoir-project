# =====================================================
#   WGCNA + ANNOTATION + KEGG (FULL VERSION)
# =====================================================

# =====================================================
# 0) LOAD CONFIG
# =====================================================
source("config/config.R")

# =====================================================
# 1) LIBRARIES
# =====================================================
library(data.table)
library(WGCNA)
library(clusterProfiler)
library(dplyr)

options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# =====================================================
#  WGCNA PARAMETERS (แก้ได้)
# =====================================================
WGCNA_CONFIG <- list(
  powers = c(1:10, seq(12,20,2)),
  deepSplit = 2,
  minClusterSize = 30,
  mergeCutHeight = 0.25,
  kME_threshold = 0.8
)

# =====================================================
# 2) LOAD TPM
# =====================================================
files <- list.files(CONFIG$tpm_dir,
                    pattern="rsem.genes.results",
                    full.names=TRUE)

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

# =====================================================
# 3) PREPARE EXPRESSION
# =====================================================
TPM_log <- log2(TPM_matrix[,-1] + 1)
rownames(TPM_log) <- TPM_matrix$Gene

datExpr <- t(TPM_log)

# Filter low expression
datExpr <- datExpr[, colMeans(datExpr) > 1]

# Filter low variance
variance <- apply(datExpr, 2, var)
cutoff <- quantile(variance, 0.5)
datExpr <- datExpr[, variance > cutoff]

# Remove bad samples/genes
gsg <- goodSamplesGenes(datExpr)
datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

cat(" Expression ready\n")

# =====================================================
# 4) NETWORK CONSTRUCTION
# =====================================================
sft <- pickSoftThreshold(datExpr,
                         powerVector = WGCNA_CONFIG$powers)

chosen_power <- sft$powerEstimate

if (is.na(chosen_power)) {
  stop("Cannot determine soft threshold")
}

cat("Chosen power:", chosen_power, "\n")

adjacency <- adjacency(datExpr, power = chosen_power)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

geneTree <- hclust(as.dist(dissTOM), method="average")

dynamicMods <- cutreeDynamic(
  dendro = geneTree,
  distM = dissTOM,
  deepSplit = WGCNA_CONFIG$deepSplit,
  pamRespectsDendro = FALSE,
  minClusterSize = WGCNA_CONFIG$minClusterSize
)

dynamicColors <- labels2colors(dynamicMods)

merge <- mergeCloseModules(datExpr,
                           dynamicColors,
                           cutHeight = WGCNA_CONFIG$mergeCutHeight)

moduleColors <- merge$colors
MEs <- merge$newMEs

cat("Modules created\n")

# =====================================================
# 5) HUB GENES
# =====================================================
kME <- signedKME(datExpr, MEs)
kME_df <- as.data.frame(kME)

kME_df$gene_id <- rownames(kME_df)
kME_df$module  <- moduleColors

# =====================================================
# 6) LOAD ANNOTATION
# =====================================================
anno <- fread(CONFIG$annotation_file,
              sep="\t",
              header=TRUE,
              data.table=FALSE)

colnames(anno)[1] <- "gene_id"
anno$gene_id <- trimws(anno$gene_id)

# =====================================================
# 7) MERGE ANNOTATION
# =====================================================
geneInfo <- data.frame(
  gene_id = colnames(datExpr),
  module = moduleColors
)

geneInfo_annotated <- merge(
  geneInfo,
  anno,
  by="gene_id",
  all.x=TRUE
)

dir.create(CONFIG$output_wgcna, showWarnings = FALSE, recursive = TRUE)

write.csv(geneInfo_annotated,
          file.path(CONFIG$output_wgcna,
                    "All_Modules_with_Annotation.csv"),
          row.names=FALSE)

# =====================================================
# 8) KEGG ENRICHMENT (ALL MODULES)
# =====================================================
modules <- unique(geneInfo_annotated$module)

for (m in modules) {
  
  module_genes <- geneInfo_annotated %>%
    filter(module == m) %>%
    filter(!is.na(KEGG_ko))
  
  if (nrow(module_genes) < 10) next
  
  module_genes$KEGG_ko <- gsub("ko:", "", module_genes$KEGG_ko)
  
  kk <- enrichKEGG(
    gene = unique(module_genes$KEGG_ko),
    organism = "ko",
    pvalueCutoff = 0.05
  )
  
  if (!is.null(kk) && nrow(as.data.frame(kk)) > 0) {
    
    write.csv(as.data.frame(kk),
              file.path(CONFIG$output_wgcna,
                        paste0(m, "_KEGG.csv")),
              row.names=FALSE)
    
    cat(" KEGG saved:", m, "\n")
  }
}

# =====================================================
# 9) EXPORT HUB GENES
# =====================================================
hub_summary <- data.frame()

for (m in modules) {
  
  kME_column <- paste0("kME", m)
  
  if (kME_column %in% colnames(kME_df)) {
    
    hub_genes <- subset(
      kME_df,
      module == m &
        kME_df[[kME_column]] > WGCNA_CONFIG$kME_threshold
    )
    
    if (nrow(hub_genes) > 0) {
      
      write.csv(hub_genes,
                file.path(CONFIG$output_wgcna,
                          paste0("Hub_", m, ".csv")),
                row.names=FALSE)
      
      hub_summary <- rbind(
        hub_summary,
        data.frame(Module=m,
                   HubGeneCount=nrow(hub_genes))
      )
    }
  }
}

write.csv(hub_summary,
          file.path(CONFIG$output_wgcna,
                    "HubGene_Summary.csv"),
          row.names=FALSE)

cat(" WGCNA FULL ANALYSIS COMPLETE\n")