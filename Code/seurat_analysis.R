#!/usr/bin/env Rscript

library(optparse)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(plotly)
library(viridis)
library(htmlwidgets)
library(SingleR)
library(celldex)
library(future)
library(reshape2)

# Parse CLI arguments
option_list <- list(
  make_option("--input_dir", type = "character"),
  make_option("--output_dir", type = "character"),
  make_option("--prefix", type = "character"),
  make_option("--resolutions", type = "character", default = "0.1,0.3,0.5,0.8,1.0"),
  make_option("--seed", type = "integer", default = 0),
  make_option("--ncores", type = "integer", default = 4)
)
opt <- parse_args(OptionParser(option_list = option_list))

# Extract parameters
input_dir <- opt$input_dir
base_output_dir <- opt$output_dir
prefix <- opt$prefix
resolutions <- as.numeric(strsplit(opt$resolutions, ",")[[1]])
seed_param <- opt$seed
ncores <- opt$ncores

# Set up parallel processing
if (ncores > 1) {
  plan("multicore", workers = ncores)
  options(future.rng.onMisuse = "ignore")
} else {
  plan("sequential")
}

# Random seed handling
if (seed_param == 0) {
  actual_seed <- sample(1:10000, 1)
  set.seed(actual_seed)
  cat("Using random seed:", actual_seed, "\n")
} else {
  actual_seed <- seed_param
  set.seed(actual_seed)
  cat("Using specified seed:", actual_seed, "\n")
}

# Clean theme
theme_clean <- theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA)
  )

# ================================================================
# ENHANCED GENE SIGNATURES
# ================================================================

get_human_gene_sets <- function() {
  return(list(
    # Functional pathways (enhanced with multiple genes)
    inflammation = c("TNF", "IL1B", "IL6", "NFKB1", "CXCL8", "CCL2"),
    stress_response = c("HSP90AA1", "HSPA1A", "HSPB1", "DNAJB1", "HSP90AB1"),
    interferon_response = c("IFIT1", "IFIT2", "IFIT3", "MX1", "ISG15", "OAS1"),
    apoptosis = c("BAX", "BAK1", "CASP3", "CASP7", "PARP1"),
    
    # Quality markers
    mitochondrial_activity = c("MT-CO1", "MT-CO2", "MT-ND1", "MT-ND4", "MT-CYB"),
    ribosomal_activity = c("RPL3", "RPS3", "RPS18", "RPL4", "RPS4X", "RPL5"),
    
    # Metabolic pathways
    glycolysis = c("HK1", "HK2", "PFKP", "ALDOA", "GAPDH", "PKM"),
    oxidative_phosphorylation = c("COX4I1", "COX5A", "COX6A1", "ATP5F1A", "NDUFB1")
  ))
}

# Biology-based resolution selection
select_best_resolution_by_biology <- function(seurat_obj, resolutions) {
  cat("  Selecting resolution based on marker gene clarity...\n")
  
  best_res <- resolutions[1]
  best_score <- 0
  
  for (res in resolutions) {
    Idents(seurat_obj) <- seurat_obj[[paste0("RNA_snn_res.", res)]]
    n_clusters <- length(unique(Idents(seurat_obj)))
    
    if (n_clusters >= 3 && n_clusters <= 15) {
      markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, max.cells.per.ident = 200, 
                               min.pct = 0.15, logfc.threshold = 0.25, verbose = FALSE)
      
      if (nrow(markers) > 0) {
        biology_score <- mean(markers$avg_log2FC) * nrow(markers) / (n_clusters * 3)
        cat("    Resolution", res, ":", n_clusters, "clusters, score =", round(biology_score, 2), "\n")
        
        if (biology_score > best_score) {
          best_score <- biology_score
          best_res <- res
        }
      }
    } else {
      cat("    Resolution", res, ":", n_clusters, "clusters (too many/few)\n")
    }
  }
  
  cat("  Selected optimal resolution:", best_res, "\n")
  return(best_res)
}

# Enhanced cell type assignment (cell-first approach)
assign_celltypes_improved <- function(seurat_obj, markers) {
  cat("  Improved cell type assignment (cell-first approach)...\n")
  
  # Step 1: Assign individual cell types from SingleR
  n_cells <- ncol(seurat_obj)
  seurat_obj$individual_celltype <- seurat_obj$singler_fine
  seurat_obj$individual_confidence <- 1.0  # SingleR confidence
  
  # Step 2: For each cluster, find consensus cell type using top marker genes
  cluster_mappings <- data.frame()
  all_clusters <- unique(Idents(seurat_obj))
  
  for (cluster_id in all_clusters) {
    cluster_cells <- Idents(seurat_obj) == cluster_id
    
    if (sum(cluster_cells) == 0) next
    
    # Get individual cell type annotations for this cluster
    cluster_individual_types <- seurat_obj$individual_celltype[cluster_cells]
    cluster_individual_types <- cluster_individual_types[!is.na(cluster_individual_types)]
    
    # Find top 2-3 most common cell types in cluster
    if (length(cluster_individual_types) > 0) {
      type_table <- table(cluster_individual_types)
      # Get top 3 most common types
      top_types <- head(names(sort(type_table, decreasing = TRUE)), 3)
      dominant_type <- top_types[1]
      confidence <- max(type_table) / sum(type_table)
      
      # Get top marker genes for this cluster
      cluster_markers <- markers[markers$cluster == cluster_id, ]
      if (nrow(cluster_markers) > 0) {
        top_markers <- head(cluster_markers$gene, 3)
        avg_logfc <- round(mean(head(cluster_markers$avg_log2FC, 3)), 2)
      } else {
        top_markers <- "No markers found"
        avg_logfc <- 0
      }
      
      cluster_mappings <- rbind(cluster_mappings, data.frame(
        Cluster = cluster_id,
        N_cells = sum(cluster_cells),
        Dominant_Type = dominant_type,
        Confidence = round(confidence, 3),
        Alternative_Types = paste(top_types[-1], collapse = ", "),
        Top_Markers = paste(top_markers, collapse = ", "),
        Avg_LogFC = avg_logfc
      ))
    } else {
      cluster_mappings <- rbind(cluster_mappings, data.frame(
        Cluster = cluster_id,
        N_cells = sum(cluster_cells),
        Dominant_Type = "Unknown",
        Confidence = 0,
        Alternative_Types = "",
        Top_Markers = "",
        Avg_LogFC = 0
      ))
    }
  }
  
  # Step 3: Assign final cluster-based cell types
  seurat_obj$data_driven_celltype <- rep("Unknown", n_cells)
  seurat_obj$celltype_confidence <- rep(0, n_cells)
  
  names(seurat_obj$data_driven_celltype) <- colnames(seurat_obj)
  names(seurat_obj$celltype_confidence) <- colnames(seurat_obj)
  
  for (i in 1:nrow(cluster_mappings)) {
    cluster_id <- cluster_mappings$Cluster[i]
    cluster_cells <- Idents(seurat_obj) == cluster_id
    
    celltype <- cluster_mappings$Dominant_Type[i]
    confidence <- cluster_mappings$Confidence[i]
    
    cell_names_in_cluster <- colnames(seurat_obj)[cluster_cells]
    seurat_obj$data_driven_celltype[cell_names_in_cluster] <- celltype
    seurat_obj$celltype_confidence[cell_names_in_cluster] <- confidence
  }
  
  cat("  Cell type assignment completed\n")
  cat("  Individual annotations preserved in 'individual_celltype'\n")
  cat("  Cluster consensus in 'data_driven_celltype'\n")
  
  return(list(seurat_obj = seurat_obj, cluster_mappings = cluster_mappings))
}

# Enhanced G0 detection using CellCycleScoring
detect_g0_enhanced <- function(seurat_obj) {
  cat("  Enhanced G0 detection using Seurat's built-in genes...\n")
  
  s_genes <- cc.genes$s.genes
  g2m_genes <- cc.genes$g2m.genes
  
  if (length(s_genes) >= 3 && length(g2m_genes) >= 3) {
    # Use Seurat's CellCycleScoring function
    seurat_obj <- CellCycleScoring(seurat_obj, 
                                  s.features = s_genes, 
                                  g2m.features = g2m_genes, 
                                  verbose = FALSE)
    
    # Enhanced G0 detection with more stringent thresholds
    low_s_threshold <- quantile(seurat_obj$S.Score, 0.25)
    low_g2m_threshold <- quantile(seurat_obj$G2M.Score, 0.25)
    
    seurat_obj$Phase_enhanced <- as.character(seurat_obj$Phase)
    
    # G0 = cells with low S AND low G2M scores (bottom 25%)
    g0_cells <- seurat_obj$S.Score < low_s_threshold & seurat_obj$G2M.Score < low_g2m_threshold
    seurat_obj$Phase_enhanced[g0_cells] <- "G0"
    seurat_obj$G0_cells <- g0_cells
    
    # Clean up any problematic values before factor conversion
    seurat_obj$Phase_enhanced[is.na(seurat_obj$Phase_enhanced) | 
                             seurat_obj$Phase_enhanced == ""] <- "G1"
    
    # Convert to factor with proper order
    seurat_obj$Phase_enhanced <- factor(seurat_obj$Phase_enhanced, 
                                       levels = c("G0", "G1", "S", "G2M"))
    
    phase_counts <- table(seurat_obj$Phase_enhanced)
    cat("    Cell cycle phases detected:\n")
    for (phase in names(phase_counts)) {
      cat("      ", phase, ":", phase_counts[phase], "cells (", 
          round(100*phase_counts[phase]/ncol(seurat_obj), 1), "%)\n")
    }
    
  } else {
    cat("    Insufficient cycle genes found\n")
    seurat_obj$Phase_enhanced <- factor(rep("G1", ncol(seurat_obj)), 
                                       levels = c("G0", "G1", "S", "G2M"))
    seurat_obj$G0_cells <- rep(FALSE, ncol(seurat_obj))
  }
  
  return(seurat_obj)
}

# ================================================================
# MAIN ANALYSIS
# ================================================================

timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M")
run_dir <- file.path(base_output_dir, paste0(prefix, "_", timestamp))

dirs <- list(
  qc = file.path(run_dir, "01_QC"),
  clustering = file.path(run_dir, "02_Clustering"), 
  markers = file.path(run_dir, "03_Markers"),
  functional = file.path(run_dir, "04_Functional"),
  interactive = file.path(run_dir, "05_Interactive")
)

lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

cat("=== ENHANCED DATA-DRIVEN PIPELINE FOR HUMAN CELLS ===\n\n")

# ================================================================
# LOAD DATA AND QC
# ================================================================

cat("Loading data...\n")
data <- Read10X(input_dir, gene.column = 2, strip.suffix = TRUE)
seurat_obj <- CreateSeuratObject(counts = data, project = prefix, min.cells = 3, min.features = 200)
cat("  Loaded", ncol(seurat_obj), "cells and", nrow(seurat_obj), "genes\n")

# QC metrics
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Adaptive QC thresholds
nf_low <- quantile(seurat_obj$nFeature_RNA, 0.05)
nf_high <- quantile(seurat_obj$nFeature_RNA, 0.95)
mt_high <- quantile(seurat_obj$percent.mt, 0.95)

seurat_obj$qc_pass <- seurat_obj$nFeature_RNA > nf_low & 
                     seurat_obj$nFeature_RNA < nf_high & 
                     seurat_obj$percent.mt < mt_high

cells_before <- ncol(seurat_obj)

cat("  QC thresholds: nFeature >", round(nf_low), ", <", round(nf_high), ", MT% <", round(mt_high, 1), "%\n")

# QC plots
p_qc_raw <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) +
  plot_annotation(title = "Raw Data QC Metrics") & theme_clean
ggsave(file.path(dirs$qc, "QC_raw_data.png"), p_qc_raw, width = 12, height = 6, bg = "white")

p_boundaries <- FeatureScatter(seurat_obj, "nCount_RNA", "nFeature_RNA", group.by = "qc_pass", pt.size = 0.5) +
  geom_hline(yintercept = c(nf_low, nf_high), linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = "QC Filtering Boundaries", color = "Passes QC") + theme_clean
ggsave(file.path(dirs$qc, "QC_boundaries.png"), p_boundaries, width = 10, height = 8, bg = "white")

# Apply filtering
seurat_obj <- subset(seurat_obj, subset = qc_pass == TRUE)
cells_after <- ncol(seurat_obj)

cat("  Kept", cells_after, "of", cells_before, "cells (", round(100*cells_after/cells_before, 1), "%)\n")

p_qc_filtered <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) +
  plot_annotation(title = "Quality Control (After Filtering)") & theme_clean
ggsave(file.path(dirs$qc, "QC_filtered_data.png"), p_qc_filtered, width = 12, height = 6, bg = "white")

# ================================================================
# NORMALIZATION AND DIMENSIONALITY REDUCTION
# ================================================================

cat("Normalization and dimensionality reduction...\n")
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 3000, verbose = FALSE)
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)

seurat_obj <- RunPCA(seurat_obj, npcs = 50, verbose = FALSE)
optimal_pcs <- min(40, ncol(seurat_obj) - 1)

seurat_obj <- RunUMAP(seurat_obj, dims = 1:optimal_pcs, 
                     n.neighbors = 15,
                     min.dist = 0.1,
                     seed.use = actual_seed, verbose = FALSE)
perplexity <- min(30, floor(ncol(seurat_obj) / 5))
seurat_obj <- RunTSNE(seurat_obj, dims = 1:optimal_pcs, perplexity = perplexity, 
                      seed.use = actual_seed, verbose = FALSE)

cat("  Using", optimal_pcs, "PCs\n")

# ================================================================
# CLUSTERING
# ================================================================

cat("Biology-driven clustering...\n")
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:optimal_pcs, verbose = FALSE)

for (res in resolutions) {
  seurat_obj <- FindClusters(seurat_obj, resolution = res, random.seed = actual_seed, verbose = FALSE)
}

optimal_res <- select_best_resolution_by_biology(seurat_obj, resolutions)
Idents(seurat_obj) <- seurat_obj[[paste0("RNA_snn_res.", optimal_res)]]
n_clusters <- length(unique(Idents(seurat_obj)))

cat("  Final clustering: resolution", optimal_res, "with", n_clusters, "clusters\n")

# Clustering plots
p_umap_clusters <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, label.size = 6) +
  labs(title = paste("UMAP Clusters (res =", optimal_res, ")")) + theme_clean
ggsave(file.path(dirs$clustering, "UMAP_clusters.png"), p_umap_clusters, width = 10, height = 8, bg = "white")

p_tsne_clusters <- DimPlot(seurat_obj, reduction = "tsne", label = TRUE, label.size = 6) +
  labs(title = paste("t-SNE Clusters (res =", optimal_res, ")")) + theme_clean
ggsave(file.path(dirs$clustering, "tSNE_clusters.png"), p_tsne_clusters, width = 10, height = 8, bg = "white")

# ================================================================
# MODULE 3: IMPROVED CELL TYPE DEFINITION
# ================================================================

cat("Module 3: Improved cell type definition...\n")

# Find cluster markers
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, 
                         logfc.threshold = 0.5, verbose = FALSE)

# SingleR annotations (cell-by-cell first)
cat("  Getting SingleR annotations...\n")
ref <- celldex::HumanPrimaryCellAtlasData()
seurat_sce <- as.SingleCellExperiment(seurat_obj)
singler_results <- SingleR(test = seurat_sce, ref = ref, labels = ref$label.fine)

cell_names <- colnames(seurat_obj)
if (length(singler_results$labels) == length(cell_names)) {
  names(singler_results$labels) <- cell_names
  seurat_obj$singler_fine <- singler_results$labels
} else {
  stop("SingleR results length mismatch")
}

singler_broad <- SingleR(test = seurat_sce, ref = ref, labels = ref$label.main)
if (length(singler_broad$labels) == length(cell_names)) {
  names(singler_broad$labels) <- cell_names
  seurat_obj$singler_broad <- singler_broad$labels
} else {
  stop("SingleR broad results length mismatch")
}

# Improved cell type assignment (cell-first approach)
celltype_analysis <- assign_celltypes_improved(seurat_obj, markers)
seurat_obj <- celltype_analysis$seurat_obj
cluster_mappings <- celltype_analysis$cluster_mappings

write.csv(cluster_mappings, file.path(dirs$markers, "cluster_to_celltype_mapping.csv"), row.names = FALSE)

seurat_obj$cell_type <- seurat_obj$data_driven_celltype

cat("  Identified", length(unique(seurat_obj$data_driven_celltype)), "cell types\n")

# Cell type visualization
p_individual_celltypes <- DimPlot(seurat_obj, group.by = "individual_celltype", reduction = "umap", 
                                 label = TRUE, repel = TRUE) +
  labs(title = "Individual Cell Annotations (SingleR)") + theme_clean + theme(legend.position = "none")
ggsave(file.path(dirs$clustering, "UMAP_individual_celltypes.png"), p_individual_celltypes, 
       width = 12, height = 10, bg = "white")

p_consensus_celltypes <- DimPlot(seurat_obj, group.by = "data_driven_celltype", reduction = "umap", 
                               label = TRUE, repel = TRUE) +
  labs(title = "Consensus Cell Types (Cluster-based)") + theme_clean + theme(legend.position = "right")
ggsave(file.path(dirs$clustering, "UMAP_data_driven_celltypes.png"), p_consensus_celltypes, 
       width = 14, height = 10, bg = "white")

p_confidence <- FeaturePlot(seurat_obj, features = "celltype_confidence", reduction = "umap") +
  scale_color_viridis_c(name = "Confidence") +
  labs(title = "Cell Type Assignment Confidence") + theme_clean
ggsave(file.path(dirs$clustering, "UMAP_celltype_confidence.png"), p_confidence, 
       width = 10, height = 8, bg = "white")

# Marker heatmaps
top_markers <- markers %>% 
  group_by(cluster) %>% 
  slice_max(n = 3, order_by = avg_log2FC)

if (nrow(top_markers) > 0) {
  p_heatmap_clusters <- DoHeatmap(seurat_obj, features = top_markers$gene, size = 4, angle = 90) +
    labs(title = "Top Marker Genes by Cluster") + theme_clean
  ggsave(file.path(dirs$markers, "marker_heatmap_by_cluster.png"), p_heatmap_clusters, 
         width = 12, height = 8, bg = "white")
}

# ================================================================
# ENHANCED FUNCTIONAL ANALYSIS
# ================================================================

cat("Enhanced functional analysis...\n")

gene_sets <- get_human_gene_sets()

# Calculate signature scores
cat("  Calculating cellular signatures...\n")

signature_results <- list()
for (sig_name in names(gene_sets)) {
  available_genes <- intersect(gene_sets[[sig_name]], rownames(seurat_obj))
  
  if (length(available_genes) >= 3) {
    seurat_obj <- AddModuleScore(seurat_obj, 
                                features = list(available_genes), 
                                name = paste0(sig_name, "_Score"),
                                nbin = 24, ctrl = 100)
    
    signature_results[[sig_name]] <- paste0(sig_name, "_Score1")
    cat("    ", sig_name, ": ", length(available_genes), " genes\n")
  }
}

# Enhanced pathway plots using signatures instead of single genes
cat("  Creating pathway signature plots...\n")

pathway_signatures <- list(
  "inflammation_Score1" = "Inflammation Response",
  "stress_response_Score1" = "Stress Response",
  "interferon_response_Score1" = "Interferon Response",
  "apoptosis_Score1" = "Apoptosis Pathway"
)

available_pathway_sigs <- intersect(names(pathway_signatures), colnames(seurat_obj[[]]))

if (length(available_pathway_sigs) >= 4) {
  p_pathways <- FeaturePlot(seurat_obj, features = available_pathway_sigs, ncol = 2, reduction = "umap") & 
    theme_clean & scale_color_viridis_c(option = "viridis")
  
  for (i in seq_along(available_pathway_sigs)) {
    p_pathways[[i]] <- p_pathways[[i]] + 
      labs(title = pathway_signatures[[available_pathway_sigs[i]]])
  }
  
  p_pathways_final <- p_pathways + 
    plot_annotation(title = "Functional Pathway Signatures",
                   subtitle = "Multi-gene signature scores for cellular pathways")
  
  ggsave(file.path(dirs$functional, "pathway_signatures.png"), p_pathways_final, 
         width = 14, height = 10, bg = "white")
}

# Universal signatures plot
if (length(signature_results) >= 4) {
  core_signatures <- c()
  core_labels <- c()
  
  priority_sigs <- c("cell_cycle_g2m", "stress_response", "inflammation", "mitochondrial_activity")
  priority_labels <- c("Cell Cycle (G2M)", "Stress Response", "Inflammation", "Mitochondrial Activity")
  
  for (i in seq_along(priority_sigs)) {
    if (priority_sigs[i] %in% names(signature_results)) {
      core_signatures <- c(core_signatures, signature_results[[priority_sigs[i]]])
      core_labels <- c(core_labels, priority_labels[i])
    }
  }
  
  if (length(core_signatures) >= 3) {
    p_universal <- FeaturePlot(seurat_obj, features = head(core_signatures, 4), 
                              ncol = 2, reduction = "umap") & 
      theme_clean & scale_color_viridis_c(option = "plasma")
    
    for (i in seq_along(head(core_signatures, 4))) {
      if (i <= length(core_labels)) {
        p_universal[[i]] <- p_universal[[i]] + labs(title = core_labels[i])
      }
    }
    
    p_universal_final <- p_universal + 
      plot_annotation(title = "Universal Cellular Signatures",
                     subtitle = "Core states generalizable across human datasets")
    
    ggsave(file.path(dirs$functional, "universal_cellular_signatures.png"), 
           p_universal_final, width = 14, height = 10, bg = "white")
  }
}

# Metabolic state analysis
if (all(c("glycolysis_Score1", "oxidative_phosphorylation_Score1") %in% colnames(seurat_obj[[]]))) {
  p_metabolism <- FeatureScatter(seurat_obj, feature1 = "glycolysis_Score1", 
                                feature2 = "oxidative_phosphorylation_Score1",
                                group.by = "data_driven_celltype", pt.size = 1) +
    labs(title = "Metabolic State Analysis", 
         subtitle = "Glycolysis vs Oxidative Phosphorylation by Cell Type",
         x = "Glycolysis Score", y = "Oxidative Phosphorylation Score") + 
    theme_clean
  ggsave(file.path(dirs$functional, "metabolic_states_by_celltype.png"), p_metabolism, 
         width = 12, height = 10, bg = "white")
}

# Enhanced cell cycle analysis
seurat_obj <- detect_g0_enhanced(seurat_obj)

if ("Phase_enhanced" %in% colnames(seurat_obj[[]]) && !all(is.na(seurat_obj$Phase_enhanced)) && 
    length(unique(seurat_obj$Phase_enhanced)) > 1) {
  
  p_phases_enhanced <- DimPlot(seurat_obj, group.by = "Phase_enhanced", reduction = "umap") +
    labs(title = "Cell Cycle Phases (Including G0)",
         subtitle = "G0 = quiescent cells with low proliferation gene expression") + theme_clean
  ggsave(file.path(dirs$functional, "cell_cycle_phases_enhanced.png"), p_phases_enhanced, 
         width = 10, height = 8, bg = "white")
  
  if ("G0_cells" %in% colnames(seurat_obj[[]])) {
    p_g0 <- DimPlot(seurat_obj, group.by = "G0_cells", reduction = "umap") +
      labs(title = "G0 (Quiescent) Cells", 
           subtitle = "Based on low S and G2M gene expression") + theme_clean
    ggsave(file.path(dirs$functional, "G0_cells.png"), p_g0, width = 10, height = 8, bg = "white")
  }
  
  # Phase distribution by cell type
  if (length(unique(seurat_obj$cell_type)) > 1) {
    phase_by_type <- table(seurat_obj$cell_type, seurat_obj$Phase_enhanced)
    
    # Create phase distribution plot
    phase_df <- as.data.frame(prop.table(phase_by_type, margin = 1))
    colnames(phase_df) <- c("Cell_Type", "Phase", "Proportion")
    
    p_phase_dist <- ggplot(phase_df, aes(x = Cell_Type, y = Proportion, fill = Phase)) +
      geom_bar(stat = "identity", position = "stack") +
      labs(title = "Cell Cycle Phase Distribution by Cell Type") +
      theme_clean + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_fill_viridis_d()
    ggsave(file.path(dirs$functional, "phase_distribution_by_celltype.png"), p_phase_dist, 
           width = 12, height = 8, bg = "white")
  }
}

# ================================================================
# POPULATION PLOTS - DUAL CELL TYPE DISTRIBUTIONS
# ================================================================

cat("Creating dual cell type distribution plots...\n")

# 1. Cluster-based cell type distribution (rename existing)
celltype_counts_cluster <- as.data.frame(table(seurat_obj$cell_type))
p_population_cluster <- ggplot(celltype_counts_cluster, aes(x = reorder(Var1, Freq), y = Freq)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
  coord_flip() +
  labs(title = "Cell Type Distribution (Cluster-based)", 
       subtitle = "Based on cluster consensus annotations",
       x = "Cell Type", y = "Number of Cells") +
  theme_clean
ggsave(file.path(dirs$functional, "cell_type_distribution_cluster_based.png"), p_population_cluster, 
       width = 10, height = 8, bg = "white")

# 2. Individual cell type distribution (SingleR-based, top 20)
individual_celltype_counts <- as.data.frame(table(seurat_obj$individual_celltype))
individual_celltype_counts <- individual_celltype_counts[!is.na(individual_celltype_counts$Var1), ]
individual_celltype_counts <- individual_celltype_counts[order(individual_celltype_counts$Freq, decreasing = TRUE), ]

# Take top 20 most abundant individual cell types
top_n_individual <- min(20, nrow(individual_celltype_counts))
top_individual_types <- head(individual_celltype_counts, top_n_individual)

p_population_individual <- ggplot(top_individual_types, aes(x = reorder(Var1, Freq), y = Freq)) +
  geom_bar(stat = "identity", fill = "darkorange", alpha = 0.7) +
  coord_flip() +
  labs(title = "Cell Type Distribution (Individual-based)", 
       subtitle = paste("Top", top_n_individual, "most abundant cell types from SingleR annotations"),
       x = "Cell Type", y = "Number of Cells") +
  theme_clean +
  theme(axis.text.y = element_text(size = 8))
ggsave(file.path(dirs$functional, "cell_type_distribution.png"), p_population_individual, 
       width = 12, height = 10, bg = "white")

cat("  Created cluster-based distribution (", nrow(celltype_counts_cluster), " types)\n")
cat("  Created individual-based distribution (top ", top_n_individual, " of ", nrow(individual_celltype_counts), " types)\n")

# ================================================================
# COMPREHENSIVE INTERACTIVE DASHBOARD
# ================================================================

cat("Module 4: Creating comprehensive interactive dashboard with ALL plots...\n")

# Prepare data for interactive plots
umap_coords <- Embeddings(seurat_obj, "umap")
tsne_coords <- Embeddings(seurat_obj, "tsne")

plot_data <- data.frame(
  UMAP1 = umap_coords[, 1],
  UMAP2 = umap_coords[, 2],
  TSNE1 = tsne_coords[, 1],
  TSNE2 = tsne_coords[, 2],
  Cluster = as.factor(Idents(seurat_obj)),
  Cell_Type = seurat_obj$data_driven_celltype,
  Individual_Type = seurat_obj$individual_celltype,
  Cell_ID = colnames(seurat_obj),
  nFeature_RNA = seurat_obj$nFeature_RNA,
  nCount_RNA = seurat_obj$nCount_RNA,
  percent_mt = seurat_obj$percent.mt,
  Confidence = seurat_obj$celltype_confidence,
  stringsAsFactors = FALSE
)

# Add cell cycle information
if ("Phase_enhanced" %in% colnames(seurat_obj[[]])) {
  plot_data$Cell_Cycle = as.character(seurat_obj$Phase_enhanced)
}
if ("S.Score" %in% colnames(seurat_obj[[]])) {
  plot_data$S_Score = seurat_obj$S.Score
}
if ("G2M.Score" %in% colnames(seurat_obj[[]])) {
  plot_data$G2M_Score = seurat_obj$G2M.Score
}

# Add available signatures
signature_cols <- grep("_Score1$", colnames(seurat_obj[[]]), value = TRUE)
for (sig_col in signature_cols) {
  clean_name <- gsub("_Score1$", "_Score", sig_col)
  plot_data[[clean_name]] <- seurat_obj[[sig_col]][, 1]
}

# Function to create plotly scatter plots
create_interactive_plot <- function(data, x_col, y_col, color_col, title_text, plot_type = "categorical") {
  
  # Create hover text
  hover_text <- paste0(
    "Cell ID: ", data$Cell_ID,
    "<br>Cluster: ", data$Cluster,
    "<br>Cell Type: ", data$Cell_Type,
    "<br>Individual Type: ", data$Individual_Type,
    "<br>nFeatures: ", data$nFeature_RNA,
    "<br>nCounts: ", data$nCount_RNA,
    "<br>MT%: ", round(data$percent_mt, 2)
  )
  
  # Add specific information based on color variable
  if (color_col %in% colnames(data)) {
    if (is.numeric(data[[color_col]])) {
      hover_text <- paste0(hover_text, "<br>", gsub("_", " ", color_col), ": ", round(data[[color_col]], 3))
    } else {
      hover_text <- paste0(hover_text, "<br>", gsub("_", " ", color_col), ": ", data[[color_col]])
    }
  }
  
  p <- plot_ly(
    data = data,
    x = ~get(x_col),
    y = ~get(y_col),
    color = ~get(color_col),
    type = "scatter",
    mode = "markers",
    marker = list(size = 4, opacity = 0.7),
    text = hover_text,
    hovertemplate = "%{text}<extra></extra>",
    showlegend = TRUE
  )
  
  # Apply appropriate color scale
  if (plot_type == "continuous") {
    p <- p %>% layout(coloraxis = list(colorscale = "Viridis"))
  }
  
  # Apply layout
  p <- p %>% layout(
    title = list(
      text = paste0("<b>", title_text, "</b><br><sup>", prefix, " - ", ncol(seurat_obj), " cells, ", n_clusters, " clusters</sup>"),
      font = list(size = 16, family = "Arial"),
      x = 0.5,
      xanchor = "center"
    ),
    xaxis = list(
      title = gsub("([12])$", " \\1", x_col),
      showgrid = TRUE,
      gridcolor = "rgba(200,200,200,0.3)",
      zeroline = FALSE,
      tickfont = list(size = 12)
    ),
    yaxis = list(
      title = gsub("([12])$", " \\1", y_col),
      showgrid = TRUE,
      gridcolor = "rgba(200,200,200,0.3)",
      zeroline = FALSE,
      tickfont = list(size = 12)
    ),
    plot_bgcolor = "white",
    paper_bgcolor = "white",
    font = list(family = "Arial, sans-serif", size = 12),
    legend = list(
      orientation = "v",
      x = 1.02,
      y = 0.5,
      xanchor = "left",
      yanchor = "middle",
      font = list(size = 10)
    ),
    autosize = FALSE,
    width = 900,
    height = 600,
    margin = list(l = 60, r = 120, t = 80, b = 60)
  ) %>%
  config(
    displayModeBar = TRUE,
    modeBarButtonsToRemove = c("pan2d", "select2d", "lasso2d", "autoScale2d"),
    displaylogo = FALSE,
    responsive = TRUE
  )
  
  return(p)
}

# Create all interactive plots
cat("  Creating interactive plots...\n")

# 1. Main dimension reduction plots
plots_list <- list()

# UMAP plots
plots_list[["umap_clusters"]] <- create_interactive_plot(plot_data, "UMAP1", "UMAP2", "Cluster", "UMAP - Clusters")
plots_list[["umap_celltypes"]] <- create_interactive_plot(plot_data, "UMAP1", "UMAP2", "Cell_Type", "UMAP - Cell Types")
plots_list[["umap_individual"]] <- create_interactive_plot(plot_data, "UMAP1", "UMAP2", "Individual_Type", "UMAP - Individual Annotations")

# t-SNE plots
plots_list[["tsne_clusters"]] <- create_interactive_plot(plot_data, "TSNE1", "TSNE2", "Cluster", "t-SNE - Clusters")
plots_list[["tsne_celltypes"]] <- create_interactive_plot(plot_data, "TSNE1", "TSNE2", "Cell_Type", "t-SNE - Cell Types")
plots_list[["tsne_individual"]] <- create_interactive_plot(plot_data, "TSNE1", "TSNE2", "Individual_Type", "t-SNE - Individual Annotations")

# 2. QC metrics plots
plots_list[["umap_nfeatures"]] <- create_interactive_plot(plot_data, "UMAP1", "UMAP2", "nFeature_RNA", "UMAP - Gene Count", "continuous")
plots_list[["umap_ncounts"]] <- create_interactive_plot(plot_data, "UMAP1", "UMAP2", "nCount_RNA", "UMAP - UMI Count", "continuous")
plots_list[["umap_mtpercent"]] <- create_interactive_plot(plot_data, "UMAP1", "UMAP2", "percent_mt", "UMAP - Mitochondrial %", "continuous")
plots_list[["umap_confidence"]] <- create_interactive_plot(plot_data, "UMAP1", "UMAP2", "Confidence", "UMAP - Cell Type Confidence", "continuous")

# 3. Cell cycle plots
if ("Cell_Cycle" %in% colnames(plot_data)) {
  plots_list[["umap_cellcycle"]] <- create_interactive_plot(plot_data, "UMAP1", "UMAP2", "Cell_Cycle", "UMAP - Cell Cycle Phases")
}
if ("S_Score" %in% colnames(plot_data)) {
  plots_list[["umap_sscore"]] <- create_interactive_plot(plot_data, "UMAP1", "UMAP2", "S_Score", "UMAP - S Phase Score", "continuous")
}
if ("G2M_Score" %in% colnames(plot_data)) {
  plots_list[["umap_g2mscore"]] <- create_interactive_plot(plot_data, "UMAP1", "UMAP2", "G2M_Score", "UMAP - G2M Phase Score", "continuous")
}

# 4. Signature plots
signature_plot_names <- grep("_Score$", colnames(plot_data), value = TRUE)
for (sig_name in signature_plot_names) {
  clean_title <- gsub("_Score$", "", sig_name)
  clean_title <- gsub("_", " ", clean_title)
  clean_title <- tools::toTitleCase(clean_title)
  
  plot_key <- paste0("umap_", tolower(gsub("_Score$", "", sig_name)))
  plots_list[[plot_key]] <- create_interactive_plot(plot_data, "UMAP1", "UMAP2", sig_name, 
                                                   paste("UMAP -", clean_title, "Signature"), "continuous")
}

# 5. Create special comparison plots
if (all(c("glycolysis_Score", "oxidative_phosphorylation_Score") %in% colnames(plot_data))) {
  # Metabolic comparison plot
  metabolic_plot <- plot_ly(
    data = plot_data,
    x = ~glycolysis_Score,
    y = ~oxidative_phosphorylation_Score,
    color = ~Cell_Type,
    type = "scatter",
    mode = "markers",
    marker = list(size = 5, opacity = 0.7),
    text = ~paste0("Cell: ", Cell_ID, "<br>Cell Type: ", Cell_Type, 
                   "<br>Glycolysis: ", round(glycolysis_Score, 3),
                   "<br>OxPhos: ", round(oxidative_phosphorylation_Score, 3)),
    hovertemplate = "%{text}<extra></extra>"
  ) %>%
  layout(
    title = list(text = "<b>Metabolic State Analysis</b><br><sup>Glycolysis vs Oxidative Phosphorylation</sup>"),
    xaxis = list(title = "Glycolysis Score"),
    yaxis = list(title = "Oxidative Phosphorylation Score"),
    plot_bgcolor = "white",
    paper_bgcolor = "white"
  )
  
  plots_list[["metabolic_comparison"]] <- metabolic_plot
}

# Save individual interactive plots
cat("  Saving individual interactive plots...\n")
for (plot_name in names(plots_list)) {
  htmlwidgets::saveWidget(plots_list[[plot_name]], 
                         file.path(dirs$interactive, paste0(plot_name, ".html")), 
                         selfcontained = TRUE)
}

# ================================================================
# CREATE COMPREHENSIVE DASHBOARD HTML
# ================================================================

create_comprehensive_dashboard <- function(plots_list, plot_data, seurat_obj) {
  
  # Organize plots into categories
  plot_categories <- list(
    "Dimension Reduction" = list(
      "UMAP - Clusters" = "umap_clusters.html",
      "UMAP - Cell Types" = "umap_celltypes.html", 
      "UMAP - Individual Types" = "umap_individual.html",
      "t-SNE - Clusters" = "tsne_clusters.html",
      "t-SNE - Cell Types" = "tsne_celltypes.html",
      "t-SNE - Individual Types" = "tsne_individual.html"
    ),
    "Quality Control" = list(
      "Gene Count" = "umap_nfeatures.html",
      "UMI Count" = "umap_ncounts.html",
      "Mitochondrial %" = "umap_mtpercent.html",
      "Cell Type Confidence" = "umap_confidence.html"
    ),
    "Cell Cycle" = list(),
    "Functional Signatures" = list(),
    "Special Analyses" = list()
  )
  
  # Add cell cycle plots if they exist
  if (file.exists(file.path(dirs$interactive, "umap_cellcycle.html"))) {
    plot_categories[["Cell Cycle"]][["Cell Cycle Phases"]] <- "umap_cellcycle.html"
  }
  if (file.exists(file.path(dirs$interactive, "umap_sscore.html"))) {
    plot_categories[["Cell Cycle"]][["S Phase Score"]] <- "umap_sscore.html"
  }
  if (file.exists(file.path(dirs$interactive, "umap_g2mscore.html"))) {
    plot_categories[["Cell Cycle"]][["G2M Phase Score"]] <- "umap_g2mscore.html"
  }
  
  # Add signature plots
  signature_files <- list.files(dirs$interactive, pattern = "umap_(inflammation|stress|interferon|apoptosis|mitochondrial|ribosomal|glycolysis|oxidative).html")
  for (sig_file in signature_files) {
    sig_name <- gsub("umap_", "", gsub(".html", "", sig_file))
    display_name <- tools::toTitleCase(gsub("_", " ", sig_name))
    plot_categories[["Functional Signatures"]][[display_name]] <- sig_file
  }
  
  # Add special analyses
  if (file.exists(file.path(dirs$interactive, "metabolic_comparison.html"))) {
    plot_categories[["Special Analyses"]][["Metabolic States"]] <- "metabolic_comparison.html"
  }
  
  # Remove empty categories
  plot_categories <- plot_categories[sapply(plot_categories, length) > 0]
  
  # Create navigation HTML
  create_navigation <- function(categories) {
    nav_html <- ""
    
    for (category_name in names(categories)) {
      nav_html <- paste0(nav_html, '
        <div class="category-group">
          <div class="category-header">', category_name, '</div>
          <div class="category-buttons">')
      
      plots_in_category <- categories[[category_name]]
      for (plot_name in names(plots_in_category)) {
        file_name <- plots_in_category[[plot_name]]
        nav_html <- paste0(nav_html, '
            <button class="nav-btn" onclick="showPlot(\'', file_name, '\', this)">', plot_name, '</button>')
      }
      
      nav_html <- paste0(nav_html, '
          </div>
        </div>')
    }
    
    return(nav_html)
  }
  
  # Generate summary statistics
  create_summary_stats <- function() {
    paste0('
      <div class="stats-grid">
        <div class="stat-item">
          <div class="stat-number">', ncol(seurat_obj), '</div>
          <div class="stat-label">Cells</div>
        </div>
        <div class="stat-item">
          <div class="stat-number">', nrow(seurat_obj), '</div>
          <div class="stat-label">Genes</div>
        </div>
        <div class="stat-item">
          <div class="stat-number">', n_clusters, '</div>
          <div class="stat-label">Clusters</div>
        </div>
        <div class="stat-item">
          <div class="stat-number">', length(unique(seurat_obj$data_driven_celltype)), '</div>
          <div class="stat-label">Cell Types</div>
        </div>
        <div class="stat-item">
          <div class="stat-number">', round(optimal_res, 2), '</div>
          <div class="stat-label">Resolution</div>
        </div>
        <div class="stat-item">
          <div class="stat-number">', round(mean(seurat_obj$nFeature_RNA)), '</div>
          <div class="stat-label">Avg Genes/Cell</div>
        </div>
      </div>')
  }
  
  # Main HTML content
  html_content <- paste0('
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Comprehensive Single-Cell Analysis Dashboard - ', prefix, '</title>
    <style>
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }
        
        body {
            font-family: "Segoe UI", "Inter", Arial, sans-serif;
            background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
            height: 100vh;
            overflow: hidden;
        }
        
        .dashboard-container {
            display: flex;
            flex-direction: column;
            height: 100vh;
        }
        
        .header {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            box-shadow: 0 4px 20px rgba(0,0,0,0.1);
        }
        
        .header-content {
            display: flex;
            justify-content: space-between;
            align-items: center;
            max-width: 1200px;
            margin: 0 auto;
        }
        
        .header-left h1 {
            font-size: 2rem;
            margin: 0;
            font-weight: 300;
        }
        
        .header-left .subtitle {
            font-size: 1rem;
            opacity: 0.9;
            margin-top: 5px;
        }
        
        .header-right {
            text-align: right;
        }
        
        .stats-grid {
            display: grid;
            grid-template-columns: repeat(3, 1fr);
            gap: 10px;
            margin-top: 10px;
        }
        
        .stat-item {
            background: rgba(255,255,255,0.2);
            padding: 8px;
            border-radius: 8px;
            text-align: center;
        }
        
        .stat-number {
            font-size: 1.2rem;
            font-weight: 600;
        }
        
        .stat-label {
            font-size: 0.7rem;
            opacity: 0.9;
        }
        
        .main-content {
            display: flex;
            flex: 1;
            overflow: hidden;
        }
        
        .sidebar {
            width: 300px;
            background: white;
            border-right: 1px solid #e1e5e9;
            overflow-y: auto;
            box-shadow: 2px 0 10px rgba(0,0,0,0.1);
        }
        
        .category-group {
            margin-bottom: 15px;
        }
        
        .category-header {
            background: #f8f9fa;
            padding: 12px 15px;
            font-weight: 600;
            font-size: 0.9rem;
            color: #495057;
            border-bottom: 1px solid #dee2e6;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }
        
        .category-buttons {
            padding: 10px;
        }
        
        .nav-btn {
            width: 100%;
            background: #f8f9fa;
            border: 1px solid #dee2e6;
            color: #495057;
            padding: 10px 15px;
            margin-bottom: 5px;
            border-radius: 6px;
            cursor: pointer;
            font-size: 0.85rem;
            transition: all 0.2s ease;
            text-align: left;
        }
        
        .nav-btn:hover {
            background: #e9ecef;
            border-color: #adb5bd;
            transform: translateX(2px);
        }
        
        .nav-btn.active {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            border-color: #667eea;
            font-weight: 500;
        }
        
        .plot-container {
            flex: 1;
            background: white;
            margin: 15px;
            border-radius: 12px;
            box-shadow: 0 4px 20px rgba(0,0,0,0.1);
            overflow: hidden;
            position: relative;
        }
        
        .plot-frame {
            width: 100%;
            height: 100%;
            border: none;
            display: none;
        }
        
        .plot-frame.active {
            display: block;
        }
        
        .loading {
            position: absolute;
            top: 50%;
            left: 50%;
            transform: translate(-50%, -50%);
            text-align: center;
            color: #6c757d;
        }
        
        .loading-spinner {
            width: 40px;
            height: 40px;
            border: 4px solid #f3f3f3;
            border-top: 4px solid #667eea;
            border-radius: 50%;
            animation: spin 1s linear infinite;
            margin: 0 auto 15px;
        }
        
        @keyframes spin {
            0% { transform: rotate(0deg); }
            100% { transform: rotate(360deg); }
        }
        
        .welcome-screen {
            display: flex;
            flex-direction: column;
            justify-content: center;
            align-items: center;
            height: 100%;
            padding: 40px;
            text-align: center;
            background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
        }
        
        .welcome-icon {
            font-size: 4rem;
            margin-bottom: 20px;
            color: #667eea;
        }
        
        .welcome-title {
            font-size: 2rem;
            color: #2c3e50;
            margin-bottom: 10px;
            font-weight: 300;
        }
        
        .welcome-subtitle {
            font-size: 1.1rem;
            color: #7f8c8d;
            margin-bottom: 30px;
        }
        
        .welcome-instructions {
            background: white;
            padding: 20px;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            max-width: 500px;
        }
        
        @media (max-width: 768px) {
            .header-content {
                flex-direction: column;
                text-align: center;
            }
            
            .stats-grid {
                grid-template-columns: repeat(2, 1fr);
                margin-top: 15px;
            }
            
            .main-content {
                flex-direction: column;
            }
            
            .sidebar {
                width: 100%;
                height: auto;
                max-height: 200px;
            }
            
            .category-buttons {
                display: flex;
                flex-wrap: wrap;
                gap: 5px;
            }
            
            .nav-btn {
                flex: 1;
                min-width: 120px;
                margin-bottom: 0;
            }
        }
    </style>
</head>
<body>
    <div class="dashboard-container">
        <div class="header">
            <div class="header-content">
                <div class="header-left">
                    <h1>🔬 Single-Cell Analysis Dashboard</h1>
                    <div class="subtitle">Dataset: ', prefix, ' | Analyzed: ', format(Sys.time(), "%Y-%m-%d %H:%M"), '</div>
                </div>
                <div class="header-right">
                    ', create_summary_stats(), '
                </div>
            </div>
        </div>
        
        <div class="main-content">
            <div class="sidebar">
                ', create_navigation(plot_categories), '
            </div>
            
            <div class="plot-container">
                <div class="welcome-screen" id="welcome-screen">
                    <div class="welcome-icon">📊</div>
                    <div class="welcome-title">Welcome to Your Analysis</div>
                    <div class="welcome-subtitle">Explore your single-cell RNA-seq data interactively</div>
                    <div class="welcome-instructions">
                        <p><strong>Getting Started:</strong></p>
                        <p>• Select any plot from the sidebar to begin</p>
                        <p>• Hover over points for detailed cell information</p>
                        <p>• Use plot controls to zoom and pan</p>
                        <p>• Switch between different views and analyses</p>
                    </div>
                </div>
                
                <div class="loading" id="loading" style="display: none;">
                    <div class="loading-spinner"></div>
                    <div>Loading visualization...</div>
                </div>')
  
  # Add iframes for all plots
  all_files <- unlist(plot_categories, use.names = FALSE)
  for (i in seq_along(all_files)) {
    html_content <- paste0(html_content, '
                <iframe class="plot-frame" src="', all_files[i], '" id="frame-', i, '"></iframe>')
  }
  
  html_content <- paste0(html_content, '
            </div>
        </div>
    </div>
    
    <script>
        function showPlot(fileName, buttonElement) {
            // Hide welcome screen
            document.getElementById("welcome-screen").style.display = "none";
            
            // Show loading
            document.getElementById("loading").style.display = "block";
            
            // Hide all frames
            const frames = document.querySelectorAll(".plot-frame");
            frames.forEach(frame => {
                frame.classList.remove("active");
            });
            
            // Remove active class from all buttons
            const buttons = document.querySelectorAll(".nav-btn");
            buttons.forEach(btn => {
                btn.classList.remove("active");
            });
            
            // Show selected frame
            const targetFrame = document.querySelector(`iframe[src="${fileName}"]`);
            if (targetFrame) {
                targetFrame.classList.add("active");
                
                // Hide loading after frame loads
                targetFrame.onload = function() {
                    document.getElementById("loading").style.display = "none";
                };
            }
            
            // Activate clicked button
            buttonElement.classList.add("active");
        }
        
        // Auto-select first plot
        window.addEventListener("load", function() {
            const firstButton = document.querySelector(".nav-btn");
            if (firstButton) {
                setTimeout(() => {
                    firstButton.click();
                }, 500);
            }
        });
    </script>
</body>
</html>')
  
  return(html_content)
}

# Create and save the comprehensive dashboard
dashboard_result <- create_comprehensive_dashboard(plots_list, plot_data, seurat_obj)
writeLines(dashboard_result, file.path(dirs$interactive, "comprehensive_dashboard.html"))

# Create a simple redirect index
index_html <- paste0('
<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <meta http-equiv="refresh" content="0; url=comprehensive_dashboard.html">
    <title>Redirecting to Comprehensive Dashboard...</title>
</head>
<body>
    <p>Redirecting to <a href="comprehensive_dashboard.html">Comprehensive Interactive Dashboard</a>...</p>
</body>
</html>')

writeLines(index_html, file.path(dirs$interactive, "index.html"))

# Count plot categories for reporting
plot_categories_count <- 5  # We know we have 5 categories: Dimension Reduction, Quality Control, Cell Cycle, Functional Signatures, Special Analyses

cat("  ✅ Comprehensive interactive dashboard created successfully!\n")
cat("  📊 Main dashboard: comprehensive_dashboard.html\n")
cat("  📈 Individual plots: ", length(plots_list), " interactive visualizations\n")
cat("  🎯 Categories: ", plot_categories_count, " plot categories organized\n")
cat("  🔗 Index redirect: index.html\n")

# ================================================================
# JSON EXPORT - KEY PLOTS AND CSV OUTPUTS
# ================================================================

cat("Creating JSON export with key plots and data files...\n")

# Load required library
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite")
}
library(jsonlite)

# Create JSON export directory
json_dir <- file.path(run_dir, "06_JSON_Export")
dir.create(json_dir, recursive = TRUE, showWarnings = FALSE)

# Function to safely convert data for JSON
prepare_for_json <- function(x) {
  if (is.factor(x)) {
    return(as.character(x))
  } else if (is.matrix(x)) {
    return(as.data.frame(x))
  } else if (any(is.infinite(x)) || any(is.nan(x))) {
    x[is.infinite(x) | is.nan(x)] <- NA
    return(x)
  }
  return(x)
}

# Function to get relative paths from run_dir
get_relative_path <- function(file_path) {
  return(gsub(paste0(run_dir, "/"), "", file_path, fixed = TRUE))
}

# ================================================================
# 1. KEY PLOTS CATALOG
# ================================================================

cat("  Cataloging key plots...\n")

plots_catalog <- list(
  qc_plots = list(
    raw_data_qc = list(
      file = "01_QC/QC_raw_data.png",
      title = "Raw Data QC Metrics",
      description = "Violin plots of nFeature_RNA, nCount_RNA, and percent.mt before filtering",
      type = "violin_plot",
      exists = file.exists(file.path(dirs$qc, "QC_raw_data.png"))
    ),
    filtering_boundaries = list(
      file = "01_QC/QC_boundaries.png", 
      title = "QC Filtering Boundaries",
      description = "Scatter plot showing filtering thresholds for feature count vs UMI count",
      type = "scatter_plot",
      exists = file.exists(file.path(dirs$qc, "QC_boundaries.png"))
    ),
    filtered_data_qc = list(
      file = "01_QC/QC_filtered_data.png",
      title = "Filtered Data QC Metrics",
      description = "Violin plots after quality control filtering",
      type = "violin_plot", 
      exists = file.exists(file.path(dirs$qc, "QC_filtered_data.png"))
    )
  ),
  
  clustering_plots = list(
    umap_clusters = list(
      file = "02_Clustering/UMAP_clusters.png",
      title = paste("UMAP Clustering (Resolution", optimal_res, ")"),
      description = paste("UMAP visualization showing", n_clusters, "clusters"),
      type = "dimension_reduction",
      method = "UMAP",
      resolution = optimal_res,
      exists = file.exists(file.path(dirs$clustering, "UMAP_clusters.png"))
    ),
    tsne_clusters = list(
      file = "02_Clustering/tSNE_clusters.png",
      title = paste("t-SNE Clustering (Resolution", optimal_res, ")"),
      description = paste("t-SNE visualization showing", n_clusters, "clusters"),
      type = "dimension_reduction",
      method = "t-SNE", 
      resolution = optimal_res,
      exists = file.exists(file.path(dirs$clustering, "tSNE_clusters.png"))
    ),
    umap_individual_celltypes = list(
      file = "02_Clustering/UMAP_individual_celltypes.png",
      title = "UMAP Individual Cell Types",
      description = "UMAP colored by individual cell type annotations from SingleR",
      type = "dimension_reduction",
      annotation_method = "SingleR",
      exists = file.exists(file.path(dirs$clustering, "UMAP_individual_celltypes.png"))
    ),
    umap_consensus_celltypes = list(
      file = "02_Clustering/UMAP_data_driven_celltypes.png", 
      title = "UMAP Consensus Cell Types",
      description = "UMAP colored by cluster consensus cell type annotations",
      type = "dimension_reduction",
      annotation_method = "cluster_consensus",
      exists = file.exists(file.path(dirs$clustering, "UMAP_data_driven_celltypes.png"))
    ),
    celltype_confidence = list(
      file = "02_Clustering/UMAP_celltype_confidence.png",
      title = "Cell Type Assignment Confidence",
      description = "UMAP colored by confidence scores for cell type assignments",
      type = "feature_plot",
      exists = file.exists(file.path(dirs$clustering, "UMAP_celltype_confidence.png"))
    )
  ),
  
  marker_plots = list(
    cluster_heatmap = list(
      file = "03_Markers/marker_heatmap_by_cluster.png",
      title = "Marker Genes Heatmap",
      description = "Heatmap showing top marker genes for each cluster",
      type = "heatmap",
      exists = file.exists(file.path(dirs$markers, "marker_heatmap_by_cluster.png"))
    )
  ),
  
  functional_plots = list(),
  
  interactive_plots = list(
    comprehensive_dashboard = list(
      file = "05_Interactive/comprehensive_dashboard.html",
      title = "Comprehensive Interactive Dashboard",
      description = "Main interactive dashboard with all plot categories organized",
      type = "interactive_html",
      exists = file.exists(file.path(dirs$interactive, "comprehensive_dashboard.html"))
    ),
    index_redirect = list(
      file = "05_Interactive/index.html",
      title = "Dashboard Index",
      description = "Index page redirecting to main dashboard",
      type = "html_redirect",
      exists = file.exists(file.path(dirs$interactive, "index.html"))
    )
  )
)

# Add all individual interactive plots
interactive_files <- list.files(dirs$interactive, pattern = "\\.html$", full.names = FALSE)
interactive_files <- interactive_files[!interactive_files %in% c("comprehensive_dashboard.html", "index.html")]

for (file_name in interactive_files) {
  plot_key <- gsub("\\.html$", "", file_name)
  clean_title <- gsub("_", " ", plot_key)
  clean_title <- tools::toTitleCase(clean_title)
  
  plots_catalog$interactive_plots[[plot_key]] <- list(
    file = paste0("05_Interactive/", file_name),
    title = paste("Interactive", clean_title),
    description = paste("Interactive plotly visualization:", clean_title),
    type = "interactive_plotly",
    exists = TRUE
  )
}

# Add functional plots that may exist
functional_plot_candidates <- list(
  pathway_signatures = list(
    file = "04_Functional/pathway_signatures.png",
    title = "Pathway Signatures",
    description = "Feature plots showing functional pathway signature scores",
    type = "feature_plot_grid"
  ),
  universal_signatures = list(
    file = "04_Functional/universal_cellular_signatures.png", 
    title = "Universal Cellular Signatures",
    description = "Core cellular signatures for cross-dataset validation",
    type = "feature_plot_grid"
  ),
  metabolic_analysis = list(
    file = "04_Functional/metabolic_states_by_celltype.png",
    title = "Metabolic States by Cell Type",
    description = "Scatter plot of glycolysis vs oxidative phosphorylation by cell type",
    type = "scatter_plot"
  ),
  cell_cycle_enhanced = list(
    file = "04_Functional/cell_cycle_phases_enhanced.png",
    title = "Enhanced Cell Cycle Phases",
    description = "UMAP showing cell cycle phases including G0 detection", 
    type = "dimension_reduction"
  ),
  g0_cells = list(
    file = "04_Functional/G0_cells.png",
    title = "G0 (Quiescent) Cells",
    description = "UMAP highlighting quiescent cells",
    type = "dimension_reduction"
  ),
  phase_distribution = list(
    file = "04_Functional/phase_distribution_by_celltype.png",
    title = "Cell Cycle Phase Distribution",
    description = "Stacked bar plot of cell cycle phases by cell type",
    type = "stacked_bar_plot"
  ),
  celltype_distribution = list(
    file = "04_Functional/cell_type_distribution.png", 
    title = "Cell Type Distribution (Individual-based)",
    description = paste("Bar plot showing top", top_n_individual, "most abundant individual cell types from SingleR"),
    type = "bar_plot"
  ),
  celltype_distribution_cluster_based = list(
    file = "04_Functional/cell_type_distribution_cluster_based.png", 
    title = "Cell Type Distribution (Cluster-based)",
    description = "Bar plot showing cell type distribution based on cluster consensus annotations",
    type = "bar_plot"
  )
)

# Check which functional plots exist and add them
for (plot_name in names(functional_plot_candidates)) {
  plot_info <- functional_plot_candidates[[plot_name]]
  full_path <- file.path(run_dir, plot_info$file)
  plot_info$exists <- file.exists(full_path)
  if (plot_info$exists) {
    plots_catalog$functional_plots[[plot_name]] <- plot_info
  }
}

# ================================================================
# 2. KEY DATA FILES CATALOG  
# ================================================================

cat("  Cataloging key data files...\n")

data_files_catalog <- list(
  main_data = list(
    seurat_object = list(
      file = "seurat_object.rds",
      title = "Seurat Object",
      description = "Complete Seurat object with all analysis results",
      type = "rds",
      size_mb = if(file.exists(file.path(run_dir, "seurat_object.rds"))) round(file.size(file.path(run_dir, "seurat_object.rds")) / 1024^2, 1) else 0,
      exists = file.exists(file.path(run_dir, "seurat_object.rds"))
    ),
    analysis_summary = list(
      file = "ANALYSIS_SUMMARY.txt",
      title = "Analysis Summary",
      description = "Text summary of analysis parameters and results",
      type = "txt",
      exists = file.exists(file.path(run_dir, "ANALYSIS_SUMMARY.txt"))
    )
  ),
  
  csv_data = list(
    cluster_mappings = list(
      file = "03_Markers/cluster_to_celltype_mapping.csv",
      title = "Cluster to Cell Type Mappings", 
      description = "CSV file mapping clusters to cell types with confidence scores and marker genes",
      type = "csv",
      columns = c("Cluster", "N_cells", "Dominant_Type", "Confidence", "Alternative_Types", "Top_Markers", "Avg_LogFC"),
      exists = file.exists(file.path(dirs$markers, "cluster_to_celltype_mapping.csv"))
    )
  ),
  
  json_exports = list(
    comprehensive_analysis = list(
      file = "06_JSON_Export/comprehensive_analysis.json",
      title = "Comprehensive Analysis JSON",
      description = "Complete analysis results in JSON format",
      type = "json"
    ),
    cell_data = list(
      file = "06_JSON_Export/cell_data.json",
      title = "Cell-level Data",
      description = "Per-cell coordinates, annotations, and metrics",
      type = "json"
    ),
    cluster_data = list(
      file = "06_JSON_Export/cluster_data.json", 
      title = "Cluster-level Data",
      description = "Per-cluster statistics and information",
      type = "json"
    ),
    plots_data = list(
      file = "06_JSON_Export/plots_data.json",
      title = "Plot Data",
      description = "Plot-ready data for web visualization",
      type = "json" 
    )
  )
)

# ================================================================
# 3. CELL-LEVEL DATA EXPORT
# ================================================================

cat("  Exporting cell-level data...\n")

# Get coordinates
umap_coords <- Embeddings(seurat_obj, "umap")
tsne_coords <- Embeddings(seurat_obj, "tsne")
pca_coords <- Embeddings(seurat_obj, "pca")[, 1:10]  # First 10 PCs

# Compile cell metadata
cell_data <- data.frame(
  cell_id = colnames(seurat_obj),
  cluster = as.character(Idents(seurat_obj)),
  cell_type = seurat_obj$data_driven_celltype,
  individual_celltype = seurat_obj$individual_celltype,
  celltype_confidence = seurat_obj$celltype_confidence,
  
  # QC metrics
  nFeature_RNA = seurat_obj$nFeature_RNA,
  nCount_RNA = seurat_obj$nCount_RNA,
  percent_mt = seurat_obj$percent.mt,
  
  # Coordinates
  UMAP_1 = umap_coords[, 1],
  UMAP_2 = umap_coords[, 2],
  tSNE_1 = tsne_coords[, 1],
  tSNE_2 = tsne_coords[, 2],
  
  stringsAsFactors = FALSE
)

# Add PCA coordinates
for (i in 1:10) {
  cell_data[[paste0("PC_", i)]] <- pca_coords[, i]
}

# Add cell cycle information if available
if ("Phase_enhanced" %in% colnames(seurat_obj[[]])) {
  cell_data$cell_cycle_phase <- as.character(seurat_obj$Phase_enhanced)
  if ("S.Score" %in% colnames(seurat_obj[[]])) {
    cell_data$S_score <- seurat_obj$S.Score
  }
  if ("G2M.Score" %in% colnames(seurat_obj[[]])) {
    cell_data$G2M_score <- seurat_obj$G2M.Score
  }
}

# Add signature scores if available
signature_cols <- grep("_Score1$", colnames(seurat_obj[[]]), value = TRUE)
for (sig_col in signature_cols) {
  clean_name <- gsub("_Score1$", "_signature", sig_col)
  cell_data[[clean_name]] <- seurat_obj[[sig_col]][, 1]
}

# Clean and prepare cell data
cell_data <- lapply(cell_data, prepare_for_json)
cell_data <- as.data.frame(cell_data)

# ================================================================
# 4. CLUSTER-LEVEL DATA EXPORT
# ================================================================

cat("  Exporting cluster-level data...\n")

# Get cluster statistics
cluster_stats <- data.frame()
for (cluster_id in sort(unique(Idents(seurat_obj)))) {
  cluster_cells <- Idents(seurat_obj) == cluster_id
  
  cluster_info <- data.frame(
    cluster_id = as.character(cluster_id),
    n_cells = sum(cluster_cells),
    dominant_cell_type = names(sort(table(seurat_obj$data_driven_celltype[cluster_cells]), decreasing = TRUE))[1],
    mean_features = round(mean(seurat_obj$nFeature_RNA[cluster_cells]), 1),
    mean_counts = round(mean(seurat_obj$nCount_RNA[cluster_cells]), 1),
    mean_mt_percent = round(mean(seurat_obj$percent.mt[cluster_cells]), 2),
    
    # Centroid coordinates
    umap_centroid_x = round(mean(umap_coords[cluster_cells, 1]), 3),
    umap_centroid_y = round(mean(umap_coords[cluster_cells, 2]), 3),
    tsne_centroid_x = round(mean(tsne_coords[cluster_cells, 1]), 3),
    tsne_centroid_y = round(mean(tsne_coords[cluster_cells, 2]), 3),
    
    stringsAsFactors = FALSE
  )
  
  cluster_stats <- rbind(cluster_stats, cluster_info)
}

# ================================================================
# 5. MARKER GENES EXPORT
# ================================================================

cat("  Exporting marker genes...\n")

# Prepare marker genes data
if (exists("markers") && nrow(markers) > 0) {
  markers_clean <- data.frame(
    gene = markers$gene,
    cluster = as.character(markers$cluster),
    avg_log2FC = round(markers$avg_log2FC, 4),
    pct_1 = round(markers$pct.1, 3),
    pct_2 = round(markers$pct.2, 3),
    p_val_adj = markers$p_val_adj,
    stringsAsFactors = FALSE
  )
  
  # Get top markers per cluster
  top_markers_per_cluster <- markers_clean %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC) %>%
    ungroup() %>%
    as.data.frame()
} else {
  markers_clean <- data.frame()
  top_markers_per_cluster <- data.frame()
}

# ================================================================
# 6. PLOTS DATA FOR WEB INTEGRATION
# ================================================================

cat("  Exporting plot data for web integration...\n")

# Create plot-ready data structures
plots_data <- list(
  umap_clusters = list(
    x = umap_coords[, 1],
    y = umap_coords[, 2],
    color = as.character(Idents(seurat_obj)),
    text = paste("Cell:", colnames(seurat_obj), "| Cluster:", Idents(seurat_obj), "| Type:", seurat_obj$data_driven_celltype),
    type = "categorical"
  ),
  
  umap_celltypes = list(
    x = umap_coords[, 1],
    y = umap_coords[, 2], 
    color = seurat_obj$data_driven_celltype,
    text = paste("Cell:", colnames(seurat_obj), "| Type:", seurat_obj$data_driven_celltype),
    type = "categorical"
  ),
  
  tsne_clusters = list(
    x = tsne_coords[, 1],
    y = tsne_coords[, 2],
    color = as.character(Idents(seurat_obj)),
    text = paste("Cell:", colnames(seurat_obj), "| Cluster:", Idents(seurat_obj), "| Type:", seurat_obj$data_driven_celltype),
    type = "categorical"
  ),
  
  tsne_celltypes = list(
    x = tsne_coords[, 1],
    y = tsne_coords[, 2],
    color = seurat_obj$data_driven_celltype,
    text = paste("Cell:", colnames(seurat_obj), "| Type:", seurat_obj$data_driven_celltype),
    type = "categorical"
  )
)

# Add signature plots if available
for (sig_col in signature_cols) {
  clean_name <- gsub("_Score1$", "", sig_col)
  plots_data[[paste0("umap_", clean_name)]] <- list(
    x = umap_coords[, 1],
    y = umap_coords[, 2],
    color = seurat_obj[[sig_col]][, 1],
    text = paste("Cell:", colnames(seurat_obj), "| Score:", round(seurat_obj[[sig_col]][, 1], 3)),
    type = "continuous"
  )
}

# ================================================================
# 7. COMPREHENSIVE JSON COMPILATION
# ================================================================

cat("  Compiling comprehensive JSON export...\n")

# Main comprehensive export
comprehensive_export <- list(
  # Analysis metadata
  metadata = list(
    dataset_info = list(
      name = prefix,
      n_cells_original = cells_before,
      n_cells_filtered = ncol(seurat_obj),
      n_genes = nrow(seurat_obj),
      filtering_rate = round(100 * ncol(seurat_obj) / cells_before, 1)
    ),
    processing_info = list(
      analysis_date = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      seed_used = actual_seed,
      optimal_resolution = optimal_res,
      n_clusters = n_clusters,
      n_cell_types_cluster_based = length(unique(seurat_obj$data_driven_celltype)),
      n_cell_types_individual = length(unique(seurat_obj$individual_celltype)),
      pcs_used = optimal_pcs,
      n_interactive_plots = length(plots_list)
    ),
    qc_thresholds = list(
      min_features = round(nf_low),
      max_features = round(nf_high),
      max_mt_percent = round(mt_high, 1)
    )
  ),
  
  # Key outputs catalog
  outputs = list(
    plots = plots_catalog,
    data_files = data_files_catalog
  ),
  
  # Analysis data
  data = list(
    cells = cell_data,
    clusters = cluster_stats,
    markers = list(
      all_markers = if(nrow(markers_clean) > 0) markers_clean else NULL,
      top_markers = if(nrow(top_markers_per_cluster) > 0) top_markers_per_cluster else NULL
    ),
    plots = plots_data
  ),
  
  # Summary statistics
  summary = list(
    cell_type_distribution_cluster_based = as.list(table(seurat_obj$data_driven_celltype)),
    cell_type_distribution_individual = as.list(table(seurat_obj$individual_celltype)),
    cluster_distribution = as.list(table(Idents(seurat_obj))),
    qc_statistics = list(
      mean_features_per_cell = round(mean(seurat_obj$nFeature_RNA), 1),
      median_features_per_cell = round(median(seurat_obj$nFeature_RNA), 1),
      mean_counts_per_cell = round(mean(seurat_obj$nCount_RNA), 1),
      median_counts_per_cell = round(median(seurat_obj$nCount_RNA), 1),
      mean_mt_percent = round(mean(seurat_obj$percent.mt), 2)
    )
  ),
  
  # Export information
  export_info = list(
    format_version = "1.0",
    generated_by = "Enhanced Single-Cell Pipeline with Comprehensive Interactive Dashboard",
    r_version = R.version.string,
    seurat_version = as.character(packageVersion("Seurat")),
    plotly_version = as.character(packageVersion("plotly")),
    export_timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
)

# ================================================================
# 8. WRITE JSON FILES
# ================================================================

cat("  Writing JSON files...\n")

# Write comprehensive export
jsonlite::write_json(
  comprehensive_export, 
  file.path(json_dir, "comprehensive_analysis.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  na = "null"
)

# Write individual components for modular access
jsonlite::write_json(
  cell_data,
  file.path(json_dir, "cell_data.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  na = "null"
)

jsonlite::write_json(
  cluster_stats,
  file.path(json_dir, "cluster_data.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  na = "null"
)

jsonlite::write_json(
  plots_data,
  file.path(json_dir, "plots_data.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  na = "null"
)

if (nrow(markers_clean) > 0) {
  jsonlite::write_json(
    markers_clean,
    file.path(json_dir, "marker_genes.json"),
    pretty = TRUE,
    auto_unbox = TRUE,
    na = "null"
  )
}

# Write plots and data files catalog
jsonlite::write_json(
  list(plots = plots_catalog, data_files = data_files_catalog),
  file.path(json_dir, "outputs_catalog.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  na = "null"
)

# ================================================================
# 9. CREATE DOCUMENTATION
# ================================================================

json_documentation <- paste0('
# JSON Export - Key Plots and Data Files

## Overview
This directory contains JSON exports of the key analytical outputs from the single-cell RNA-seq analysis.

## Key Files:

### comprehensive_analysis.json
Complete analysis results including:
- Analysis metadata and parameters
- **Plots catalog** with file paths and descriptions
- **Data files catalog** (CSV, RDS files) 
- Cell-level data and coordinates
- Cluster statistics and marker genes
- Plot-ready data for web visualization

### Individual Files:
- **cell_data.json** - Per-cell coordinates, annotations, QC metrics (', nrow(cell_data), ' cells)
- **cluster_data.json** - Cluster statistics and centroids (', nrow(cluster_stats), ' clusters)
- **plots_data.json** - Plot-ready coordinates and colors for web visualization
- **marker_genes.json** - Differentially expressed genes by cluster
- **outputs_catalog.json** - Catalog of all plots and data files with metadata

## Interactive Dashboard:
- **Comprehensive Dashboard**: ', length(plots_list), ' interactive plots organized in categories
- **Plot Categories**: Dimension Reduction, Quality Control, Cell Cycle, Functional Signatures, Special Analyses
- **Features**: Hover tooltips, zoom/pan controls, organized sidebar navigation

## Key Plots Included:
', paste0("- ", names(plots_catalog), " (", sapply(plots_catalog, length), " plots each)", collapse = "\n"), '

## Cell Type Distributions:
- **Individual-based**: Top ', top_n_individual, ' most abundant SingleR annotations
- **Cluster-based**: ', length(unique(seurat_obj$data_driven_celltype)), ' consensus cell types from cluster analysis

## Data Files Included:
- Seurat object (.rds)
- Cluster mappings (.csv) 
- Analysis summary (.txt)
- JSON exports

## Usage Example:
```javascript
// Load complete analysis
fetch("comprehensive_analysis.json")
  .then(response => response.json())
  .then(data => {
    // Access plot files
    console.log("QC plots:", data.outputs.plots.qc_plots);
    console.log("Interactive dashboard:", data.outputs.plots.interactive_plots.comprehensive_dashboard.file);
    
    // Compare cell type distributions
    console.log("Individual cell types:", data.summary.cell_type_distribution_individual);
    console.log("Cluster-based cell types:", data.summary.cell_type_distribution_cluster_based);
    
    // Use plot data
    const umap_data = data.data.plots.umap_clusters;
    createPlot(umap_data.x, umap_data.y, umap_data.color);
  });
```

Generated: ', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '
Dataset: ', prefix, ' (', ncol(seurat_obj), ' cells, ', n_clusters, ' clusters)
Cell Types: ', length(unique(seurat_obj$data_driven_celltype)), ' cluster-based, ', length(unique(seurat_obj$individual_celltype)), ' individual-based
Interactive Plots: ', length(plots_list), ' visualizations
')

writeLines(json_documentation, file.path(json_dir, "README.md"))

cat("  JSON export completed successfully\n")
cat("  Export directory:", json_dir, "\n")
cat("  Comprehensive file: comprehensive_analysis.json\n")
cat("  Plots catalog: ", length(unlist(plots_catalog, recursive = FALSE)), " plots catalogued\n")
cat("  Data files: ", length(unlist(data_files_catalog, recursive = FALSE)), " data files catalogued\n")

# ================================================================
# SAVE RESULTS
# ================================================================

cat("Saving results...\n")

saveRDS(seurat_obj, file.path(run_dir, "seurat_object.rds"))

summary_text <- c(
  "=== ENHANCED DATA-DRIVEN PIPELINE FOR HUMAN CELLS ===",
  paste("Completed:", Sys.time()),
  paste("Random seed:", actual_seed),
  "",
  "=== RESULTS ===", 
  paste("Cells analyzed:", ncol(seurat_obj), "of", cells_before),
  paste("Optimal resolution:", optimal_res),
  paste("Clusters found:", n_clusters),
  paste("Cell types identified (cluster-based):", length(unique(seurat_obj$data_driven_celltype))),
  paste("Cell types identified (individual-based):", length(unique(seurat_obj$individual_celltype))),
  paste("Universal signatures:", length(signature_results)),
  paste("Interactive plots created:", length(plots_list)),
  ""
)

writeLines(summary_text, file.path(run_dir, "ANALYSIS_SUMMARY.txt"))

cat("\n================================================================\n")
cat("    COMPREHENSIVE SINGLE-CELL PIPELINE COMPLETE\n") 
cat("================================================================\n")
cat(paste("Results:", run_dir), "\n")
cat(paste("Analyzed:", ncol(seurat_obj), "cells,", n_clusters, "clusters"), "\n")
cat(paste("Cell types (cluster-based):", length(unique(seurat_obj$data_driven_celltype))), "\n")
cat(paste("Cell types (individual-based):", length(unique(seurat_obj$individual_celltype))), "\n")
cat(paste("Interactive plots:", length(plots_list)), "\n")
cat("Comprehensive dashboard: comprehensive_dashboard.html\n")
cat("================================================================\n")