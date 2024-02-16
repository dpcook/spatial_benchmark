library(Seurat)
library(dplyr)
library(ggplot2)
library(shadowtext)
library(scales)
library(cowplot)
library(data.table)
library(Matrix)
library(matrixStats)
library(SingleCellExperiment)
library(SpatialExperiment)
library(SpatialFeatureExperiment)
library(bluster)
library(BiocParallel)

#######
# I/O
#######

# Functions assume data is stored in specific formats:
# Xenium: Typical Xenium bundle; CosMx: Flat CSV (TAP style); MERSCOPE: Not implemented yet

#####
# readSpatial() reads in data from either Xenium, CosMx, of MERSCOPE.
# It outputs a seurat object with some common metadata e for downstream comparison.
# Regardless of platform, data is stored in an assay named "RNA" for convenient
# function
readSpatial <- function(sample_id, path, platform){
  print(paste0("Reading: ", sample_id))
  
  if(platform == "Xenium"){
    print("Loading Xenium data")
    seu_obj <- LoadXenium(path, assay = "RNA")
    seu_obj@meta.data$sample_id <- sample_id
    seu_obj@meta.data$platform <- platform
    seu_obj@meta.data$path <- path #used in some functions to pull tx data
    
    ##### Cell Metadata
    print("Getting additional cell metadata")
    cell_meta <- data.table::fread(
      file.path(path, "cells.csv.gz")
    )
    #Set up a few defined metadata columns
    seu_obj@meta.data$cell_area <- cell_meta$cell_area
    seu_obj@meta.data$nucleus_area <- cell_meta$nucleus_area
    seu_obj@meta.data$transcript_counts <- cell_meta$transcript_counts
    seu_obj@meta.data$negprobe_counts <- cell_meta$control_probe_counts
    
    ### Add tissue coordinates as an embedding for custom plotting
    print("Adding tissue coordinates as embedding")
    coords <- GetTissueCoordinates(seu_obj)
    coords <- as.matrix(coords[,1:2])
    colnames(coords) <- c("Tissue_1", "Tissue_2")
    rownames(coords) <- colnames(seu_obj)
    seu_obj[["tissue"]] <- CreateDimReducObject(coords, key="Tissue_", assay="RNA")
    
  } else if(platform == "CosMx"){
    print("Loading CosMx data")
    seu_obj <- LoadNanostring(path, fov="fov")
    seu_obj$sample_id <- sample_id
    seu_obj@meta.data$platform <- platform
    seu_obj@meta.data$path <- path #used in some functions to pull tx data
    
    #Fix assays to separate targeting and non-targeting probes
    ## Add negative control probe assay
    sys_probes <- grep("SystemControl", rownames(seu_obj), value=T)
    neg_probes <- grep("Negative", rownames(seu_obj), value=T)
    seu_obj[["ControlProbe"]] <- CreateAssayObject(
      counts = seu_obj[["Nanostring"]]$counts[neg_probes,]
    )
    ## Make "Nanostring" assay
    tx_probes <- rownames(seu_obj)[!rownames(seu_obj) %in% c(sys_probes, neg_probes)]
    seu_obj[["RNA"]] <- CreateAssayObject(
      counts = seu_obj[["Nanostring"]]$counts[tx_probes,]
    )
    DefaultAssay(seu_obj) <- "RNA"
    seu_obj[["Nanostring"]] <- NULL
    
    ##### Cell Metadata
    print("Getting additional cell metadata")
    cell_meta <- data.table::fread(
      file.path(path, list.files(path, pattern="*metadata_file.csv.gz"))
    )
    #It's excessive, but we'll add all metadata into the object
    seu_obj@meta.data <- cbind(seu_obj@meta.data, cell_meta)
    seu_obj@meta.data$fov <- factor(paste0("FOV", seu_obj@meta.data$fov))
    seu_obj@meta.data$cell_area <- seu_obj$Area.um2
    seu_obj@meta.data$transcript_counts <- seu_obj$nCount_RNA
    seu_obj@meta.data$negprobe_counts <- seu_obj$nCount_ControlProbe
    
    ### Add tissue coordinates as an embedding for custom plotting
    print("Adding tissue coordinates as embedding")
    coords <- data.frame(
      Tissue_1 = cell_meta$CenterY_global_px,
      Tissue_2 = cell_meta$CenterX_global_px
    )
    coords <- as.matrix(coords[,1:2])
    colnames(coords) <- c("Tissue_1", "Tissue_2")
    rownames(coords) <- colnames(seu_obj)
    seu_obj[["tissue"]] <- CreateDimReducObject(coords, key="Tissue_", assay="RNA")

    
  } else if(paltform == "Merscope"){
    print("Working on support!")
    stop()
    
  } else{
    print("Not a supported platform")
    stop()
    
  }
  
  return(seu_obj)
}

#####
# readTxMeta() simply reads in the transcript localization/metadata table
# for each platform. This table will be used by subsequent functions
readTxMeta <- function(path, platform){
  if(platform == "Xenium"){
    df <- data.table::fread(file.path(path, "transcripts.csv.gz"))
  } else if(platform == "CosMx"){
    df <- data.table::fread(file.path(path, 
                                   list.files(path, pattern = "*tx_file.csv.gz")))
  } else if(platform == "Merscope"){
    print("Working on support!")
    stop()
  } else{
    print("Platform not supported")
    stop()
  }
}

#######
# QC
#######

### Transcripts per cell

getTxPerCell <- function(seu_obj, #features can be explicitly defined. Defaults to all targets
                         features=NULL){ 
  if(is.null(features)){
    features <- rownames(seu_obj)
  } else{
    features <- features
  }
  
  mean_tx <- mean(colSums(seu_obj[["RNA"]]$counts[features,]))
  
  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=mean_tx
  )
  return(res)
} 

### Transcripts per um2
getTxPerArea <- function(seu_obj, 
                         features=NULL){ 
  if(is.null(features)){
    features <- rownames(seu_obj)
  } else{
    features <- features
  }
  
  tx_count <- colSums(seu_obj[["RNA"]]$counts[features,])
  mean_tx_norm <- mean(tx_count / seu_obj$cell_area)
  
  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=mean_tx_norm
  )
  return(res)
}

### Transcripts per nucleus
getTxPerNuc <- function(seu_obj, 
                        features=NULL){
  if(is.null(features)){
    features <- rownames(seu_obj)
  } else{
    features <- features
  }
  
  path <- unique(seu_obj$path)
  platform <- unique(seu_obj$platform)
  
  # Read Tx localization data
  tx_df <- readTxMeta(path, platform)
  
  if(platform == "Xenium"){
    tx_df <- filter(tx_df, cell_id %in% colnames(seu_obj) & 
                      overlaps_nucleus == 1 &
                      feature_name %in% features) %>%
      group_by(cell_id) %>%
      summarize(nuc_counts = n())
    
    
  } else if(platform == "CosMx"){
    tx_df$cell_id <- paste(tx_df$cell_ID, tx_df$fov, sep="_")
    tx_df <- tx_df %>%
      filter(cell_id %in% colnames(seu_obj) & 
               CellComp == "Nuclear" &
               target %in% features) %>%
      group_by(cell_id) %>%
      summarize(nuc_counts = n())
    
  } else if(platform == "Merscope"){
    print("Working on support")
    
  } else{
    print("Platform not supported")
  }
  
  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=mean(tx_df$nuc_counts)
  )
  
  return(res)
}

### Transcripts per nucleus
getTxPerNuc <- function(seu_obj, 
                        features=NULL){
  if(is.null(features)){
    features <- rownames(seu_obj)
  } else{
    features <- features
  }
  
  path <- unique(seu_obj$path)
  platform <- unique(seu_obj$platform)
  
  # Read Tx localization data
  tx_df <- readTxMeta(path, platform)
  
  if(platform == "Xenium"){
    tx_df <- filter(tx_df, cell_id %in% colnames(seu_obj) & 
                      overlaps_nucleus == 1 &
                      feature_name %in% features) %>%
      group_by(cell_id) %>%
      summarize(nuc_counts = n())
    
    
  } else if(platform == "CosMx"){
    tx_df$cell_id <- paste(tx_df$cell_ID, tx_df$fov, sep="_")
    tx_df <- tx_df %>%
      filter(cell_id %in% colnames(seu_obj) & 
               CellComp == "Nuclear" &
               target %in% features) %>%
      group_by(cell_id) %>%
      summarize(nuc_counts = n())
    
  } else if(platform == "Merscope"){
    print("Working on support")
    
  } else{
    print("Platform not supported")
  }
  
  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=mean(tx_df$nuc_counts)
  )
  
  return(res)
}

### Per Probe Mean Expression
getMeanExpression <- function(seu_obj, 
                              features=NULL){
  if(is.null(features)){
    features <- rownames(seu_obj)
  } else{
    features <- features
  }
  
  target_df <- data.frame(
    target = features,
    value = rowMeans(seu_obj[["RNA"]]$counts[features,]),
    type = "Gene"
  )
  
  control_df <- data.frame(
    target = rownames(seu_obj[["ControlProbe"]]$counts),
    value = rowMeans(seu_obj[["ControlProbe"]]$counts),
    type = "Control"
  )
  
  res <- rbind(target_df, control_df)
  res$platform <- unique(seu_obj$platform)
  res$sample_id <- unique(seu_obj$sample_id)
  
  return(res)
} 

### log-ratio of mean gene counts to mean neg probe counts
getMeanSignalRatio <- function(seu_obj, 
                               features=NULL){
  if(is.null(features)){
    features <- rownames(seu_obj)
  } else{
    features <- features
  }
  
  tx_means <- rowMeans(seu_obj[["RNA"]]$counts[features,])
  neg_probe_means <- rowMeans(seu_obj[["ControlProbe"]]$counts)
  
  ratio <- log10(tx_means) - log10(mean(neg_probe_means))
  ratio <- mean(ratio)
  
  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=ratio
  )
  
  return(res)
}

### Fraction of transcripts in cells
getCellTxFraction <- function(seu_obj, 
                              features=NULL){
  if(is.null(features)){
    features <- rownames(seu_obj)
  } else{
    features <- features
  }
  
  path <- unique(seu_obj$path)
  platform <- unique(seu_obj$platform)
  
  tx_df <- readTxMeta(path, platform)
  
  if(platform == "Xenium"){
    tx_df <- filter(tx_df, feature_name %in% features)
    total_tx_count <- nrow(tx_df)
    unassigned_tx_count <- sum(tx_df$cell_id == "UNASSIGNED")
    
    cell_tx_fraction <- (total_tx_count - unassigned_tx_count) / total_tx_count
    
  } else if(platform == "CosMx"){
    tx_df <- filter(tx_df, target %in% features)
    total_tx_count <- nrow(tx_df)
    unassigned_tx_count <- sum(tx_df$CellComp == "None")
    
    cell_tx_fraction <- (total_tx_count - unassigned_tx_count) / total_tx_count
    
  } else if(platform == "Merscope"){
    print("Working on support")
    
  } else{
    print("Platform not supported")
  }
  
  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=cell_tx_fraction
  )
  
  return(res)
}

##### Dynamic Range 
# Log-ratio of highest mean exp vs. mean noise
getMaxRatio <- function(seu_obj, 
                        features=NULL){
  if(is.null(features)){
    features <- rownames(seu_obj)
  } else{
    features <- features
  }
  
  tx_means <- rowMeans(seu_obj[["RNA"]]$counts[features,])
  neg_probe_means <- rowMeans(seu_obj[["ControlProbe"]]$counts)
  
  ratio <- log10(max(tx_means)) - log10(mean(neg_probe_means))
  
  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=ratio
  )
  
  return(res)
}

# Distribution of maximal values
getMaxDetection <- function(seu_obj, 
                            features=NULL){
  if(is.null(features)){
    features <- rownames(seu_obj)
  } else{
    features <- features
  }
  
  max_vals <- matrixStats::rowMaxs(as.matrix(seu_obj[["RNA"]]$counts[features,]))
  
  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=max_vals,
    gene = features
  )
  
  return(res)
}

##### Mutually Exclusive Co-expression Rate (MECR) Implementation
getMECR <- function(seu_obj) {
  #This function comes from Hartman & Satija, bioRxiv, 2024
  #We are using a custom marker table. The original publication bases it on
  #scRNA-seq from matched tissue.
  marker_df <- data.frame(
    gene = c("EPCAM", "KRT19", "KRT8", 
             "CD3E", "CD3D", "CD8A", "NKG7",
             "MS4A1", "CD79A",
             "PECAM1", "CLDN5", "VWF",
             "C1QA", "C1QB", "CD14", "FCGR3A", "ITGAX", "ITGAM",
             "PDGFRA", "DPT", "COL1A1",
             "MYH11", "ACTG2"),
    cell_type = c("Epithelial", "Epithelial", "Epithelial",
                  "T", "T", "T", "T",
                  "B", "B",
                  "Endo", "Endo", "Endo",
                  "Macro", "Macro", "Macro", "Macro", "Macro", "Macro",
                  "Fibro", "Fibro", "Fibro",
                  "Muscle", "Muscle")
  )
  rownames(marker_df) <- marker_df$gene
  
  coexp.rates <- c()
  genes <- intersect(rownames(seu_obj), rownames(marker_df))
  print(paste0("Marker count: ", length(genes)))
  if (length(genes) > 25) { genes <- sample(genes, 25) }
  mtx <- as.matrix(seu_obj[['RNA']]$counts[genes, ])
  for (g1 in genes) {
    for (g2 in genes) {
      if ((g1 != g2) && (g1 > g2) && (marker_df[g1, "cell_type"] != marker_df[g2, "cell_type"])) {
        c1 <- mtx[g1, ]
        c2 <- mtx[g2, ]
        coexp.rates <- c(
          coexp.rates,
          sum(c1 > 0 & c2 > 0) / sum(c1 > 0 | c2 > 0)) # >0 too liberal of an expression threshold?
      }
    }
  }
  
  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=round(mean(coexp.rates), digits=3)
  )
  
  return(res)
}

##### Distribution spatial autocorrelation
getMorans <- function(seu_obj, 
                      features=NULL){
  #Requires SingleCellExperiment, SpatialFeatureExperiment, Voyager, scater
  
  if(is.null(features)){
    features <- rownames(seu_obj)
  } else{
    features <- features
  }
  
  #First run for gene-targeting probes
  print("Getting Moran's I for gene-targeting probes")
  sce <- SingleCellExperiment(list(counts=seu_obj[["RNA"]]$counts[features,]),
                              colData = seu_obj@meta.data)
  colData(sce) <- cbind(colData(sce), Embeddings(seu_obj, 'tissue'))
  spe <- toSpatialExperiment(sce, spatialCoordsNames = c("Tissue_1",
                                                         "Tissue_2"))
  sfe <- toSpatialFeatureExperiment(spe)
  sfe <- sfe[, colSums(counts(sfe)) > 0]
  rowData(sfe)$means <- rowMeans(counts(sfe))
  rowData(sfe)$vars <- rowVars(counts(sfe))
  sfe <- scater::logNormCounts(sfe)
  
  colGraph(sfe, "knn20") <- findSpatialNeighbors(sfe, method = "knearneigh", 
                                                 dist_type = "idw", k = 20, 
                                                 style = "W")
  
  sfe <- Voyager::runMoransI(sfe, colGraphName = "knn20", BPPARAM = MulticoreParam(8))
  
  spatial_cor <- as.data.frame(rowData(sfe))
  
  targeting <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=spatial_cor[,3], #morans I
    gene = rownames(spatial_cor),
    type = "Gene"
  )
  
  #Now run for control probes
  print("Getting Moran's I for non-targeting probes")
  sce <- SingleCellExperiment(list(counts=seu_obj[["ControlProbe"]]$counts),
                              colData = seu_obj@meta.data)
  colData(sce) <- cbind(colData(sce), Embeddings(seu_obj, 'tissue'))
  spe <- toSpatialExperiment(sce, spatialCoordsNames = c("Tissue_1",
                                                         "Tissue_2"))
  sfe <- toSpatialFeatureExperiment(spe)
  sfe <- sfe[, colSums(counts(sfe)) > 0]
  rowData(sfe)$means <- rowMeans(counts(sfe))
  rowData(sfe)$vars <- rowVars(counts(sfe))
  sfe <- scater::logNormCounts(sfe)
  
  #Nearest neighbor
  colGraph(sfe, "knn20") <- findSpatialNeighbors(sfe, method = "knearneigh", 
                                                 dist_type = "idw", k = 20, 
                                                 style = "W")
  #Moran's I
  sfe <- Voyager::runMoransI(sfe, colGraphName = "knn20", BPPARAM = MulticoreParam(8))
  
  spatial_cor <- as.data.frame(rowData(sfe))
  
  control <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=spatial_cor[,3], #morans I
    gene = rownames(spatial_cor),
    type = "Control"
  )
  
  res <- rbind(targeting, control)
  
  return(res)
}

##### Cluster evaluation: silhouette width
getSilhouetteWidth <- function(seu_obj){
  print("Clustering data")
  seu_obj <- seu_obj %>%
    NormalizeData() %>%
    ScaleData()
  VariableFeatures(seu_obj) <- rownames(seu_obj)
  seu_obj <- seu_obj %>%
    RunPCA(verbose=F) %>%
    FindNeighbors(dims=1:10) %>%
    FindClusters(resolution=0.2)
  
  #Downsample to 100 cells per cluster for silhouette calculation
  seu_obj <- subset(seu_obj, downsample = 10000)
  
  silhouette <- bluster::approxSilhouette(
    Embeddings(seu_obj, 'pca')[,1:10],
    clusters = seu_obj$seurat_clusters
  )
  
  silhouette <- as.data.frame(silhouette)
  
  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=round(mean(silhouette$width), digits=3)
  )
  
}



#######
# Plotting
#######
# All plots assume input is a tidy data frame with the following columns:
# 1) sample_id
# 2) platform
# 3) value (based on what is being plotted--from functions above)

plotSampleLabel <- function(sample_meta){
  df <- data.frame(
    samples = sample_meta$sample_id,
    platform = sample_meta$platform
  )
  
  p <- ggplot(df, aes(x="", y=samples)) +
    geom_text(aes(label=samples, color=platform), size=5, hjust=0.5) +
    scale_color_manual(values = c("#59C134", "#14B3E6")) + 
    theme_void() + theme(legend.position='none')
  return(p)
}

plotPanelSize <- function(df){
  p <- ggplot(df, aes(x="", y=sample_id)) +
    geom_point(shape=21, color='black', alpha=0.8, stroke=1,
               aes(size=value, fill=value)) +
    geom_shadowtext(color = "black", size = 4, #fontface = "bold", 
                    bg.colour = "white", bg.r = .2,
                    aes(label=scales::comma(value))) +
    scale_fill_gradientn(colours=viridis::mako(100)) +
    xlab("Panel size") + ylab("") +
    scale_size(range = c(6,12)) +
    scale_x_discrete(position='top',
                     labels = c("")) +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  return(p)
  
}

plotCellCount <- function(df){
  p <- ggplot(df, aes(x="", y=sample_id)) +
    geom_text(aes(label=scales::comma(value), color=platform), size=5, hjust=0.5) +
    scale_color_manual(values = c("#59C134", "#14B3E6")) + 
    scale_x_discrete(position='top') +
    xlab("Cell count") + ylab("") +
    theme_classic() + 
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=10, color="black"),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  return(p)
}

plotTxPerCell <- function(df){
  p <- ggplot(df, aes(x="", y=sample_id)) +
    geom_point(shape=21, color='black', alpha=0.8, stroke=1,
               aes(size=value, fill=value)) +
    geom_shadowtext(color = "black", size = 4, #fontface = "bold", 
                    bg.colour = "white", bg.r = .2,
                    aes(label=scales::comma(value))) +
    scale_fill_gradientn(colours=viridis::mako(100)) +
    xlab("Per cell") + ylab("") +
    scale_size(range = c(7,12)) +
    scale_x_discrete(position='top',
                     labels = c("")) +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_text(size=10),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  return(p)
}

plotTxPerArea <- function(df){
  p <- ggplot(df, aes(x="", y=sample_id)) +
    geom_point(shape=21, color='black', alpha=0.8, stroke=1,
               aes(size=value, fill=value)) +
    geom_shadowtext(color = "black", size = 4, #fontface = "bold", 
                    bg.colour = "white", bg.r = .2,
                    aes(label=scales::comma(value))) +
    xlab("Per cell\num^2") + ylab("") +
    scale_fill_gradientn(colours=viridis::mako(100)) +
    scale_size(range = c(7,12)) +
    scale_x_discrete(position='top') +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_text(size=10),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  return(p)
}

plotTxPerNuc <- function(df){
  df$column <- ""
  p <- ggplot(df, aes(x=column, y=sample_id)) +
    geom_point(shape=21, color='black', alpha=0.8, stroke=1,
               aes(size=value, fill=value)) +
    geom_shadowtext(color = "black", size = 4, #fontface = "bold", 
                    bg.colour = "white", bg.r = .2,
                    aes(label=scales::comma(value))) +
    xlab("Per\nnucleus") + ylab("") +
    scale_fill_gradientn(colours=viridis::mako(100)) +
    scale_size(range = c(7,12)) +
    scale_x_discrete(position='top') +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_text(size=10),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  return(p)
}

plotTxPerCellNorm <- function(df){
  p <- ggplot(df, aes(x="", y=sample_id)) +
    geom_point(shape=21, color='black', alpha=0.8, stroke=1,
               aes(size=value, fill=value)) +
    geom_shadowtext(color = "black", size = 4, #fontface = "bold", 
                    bg.colour = "white", bg.r = .2,
                    aes(label=scales::comma(value))) +
    xlab("Per cell\nper target") + ylab("") +
    scale_fill_gradientn(colours=viridis::mako(100)) +
    scale_size(range = c(7,12)) +
    scale_x_discrete(position='top') +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_text(size=10),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  return(p)
}

plotTxPerCellIntersect <- 

plotFractionTxInCell <- function(df){
  p <- ggplot(df, aes(x=value, y=sample_id)) +
    geom_col(color='black', fill='grey90', stroke=1) +
    geom_text(aes(label=scales::comma(value)),
              hjust = 1, nudge_x = -.05) +
    xlab("Fraction Tx in Cells") + ylab("") +
    scale_x_continuous(position='top',
                       expand = c(0,0),
                       breaks=c(0, 0.5, 1),
                       labels = c(0, 0.5, 1),
                       limits=c(0,1)) +
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=10, color="black"),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  return(p)
}

plotTxPerCell_Intersect <- function()

plotSignalRatio <- function(df){
  p <- ggplot(df, aes(x=value, y=sample_id)) +
    geom_col(color='black', fill='grey90') +
    geom_text(aes(label=scales::comma(value)),
              hjust = 1, nudge_x = -.05) +
    xlab("Mean log10-ratio\nexpression over noise") + ylab("") +
    scale_x_continuous(position='top',
                       expand = c(0,0)) +
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=10, color="black"),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  return(p)
}

plotMeanExpression <- function(df){
  p <- ggplot(df, aes(x=value, y=sample_id)) +
    geom_jitter(size=0.15, shape=16, aes(color=type),
                position = position_jitterdodge(jitter.width=0.25)) +
    geom_boxplot(color="black", 
                 alpha=0.5, outlier.size=0, outlier.colour = NA,
                 aes(fill=type)) +
    scale_fill_manual(values=c("lightgrey", "firebrick")) +
    scale_colour_manual(values=c("grey20", "firebrick")) +
    xlab("Mean gene\n detection per cell") + ylab("") +
    scale_x_log10(position='top', expand = c(0,0),
                  labels = label_log(digits = 2)) +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=10, color="black"),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  return(p)
}

plotMaxExpression <- function(df){
  p <- ggplot(df, aes(x=value, y=sample_id)) +
    geom_jitter(size=0.15, shape=16, aes(color=platform)) +
    geom_boxplot(color="black", fill="lightgrey",
                 alpha=0.5, outlier.size=0, outlier.colour = NA) +
    scale_colour_manual(values = c("#59C134", "#14B3E6")) +
    xlab("Maximal gene\ndetection per cell") + ylab("") +
    scale_x_continuous(position='top', expand = c(0,0),
                  limits = c(0,50), oob=squish) +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=10, color="black"),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
}

plotMECR <- function(df){
  p <- ggplot(df, aes(x="", y=sample_id)) +
    geom_point(shape=21, color='black', alpha=0.8, stroke=1,
               aes(size=value, fill=value)) +
    geom_shadowtext(color = "black", size = 4, #fontface = "bold", 
                    bg.colour = "white", bg.r = .2,
                    aes(label=scales::comma(value))) +
    scale_fill_gradientn(colours=RColorBrewer::brewer.pal(9, "YlOrRd"),
                         limits=c(0.01)) +
    xlab("MECR") + ylab("") +
    scale_size(range = c(7,12), limits = c(0, 0.1)) +
    scale_x_discrete(position='top',
                     labels = c("")) +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_text(size=10),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  return(p)
}

plotMorans <- function(df){
  p <- ggplot(df, aes(x=value, y=sample_id)) +
    geom_jitter(size=0.15, shape=16, aes(color=type),
                position = position_jitterdodge(jitter.width=0.25)) +
    geom_boxplot(color="black", 
                 alpha=0.5, outlier.size=0, outlier.colour = NA,
                 aes(fill=type)) +
    scale_fill_manual(values=c("lightgrey", "firebrick")) +
    scale_colour_manual(values=c("grey20", "firebrick")) +
    xlab("Mean gene\n spatial autocorrelation\n(Moran's I)") + ylab("") +
    scale_x_continuous(position='top', expand = c(0,0)) +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=10, color="black"),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  return(p)
}

plotSilhouette <- function(df){
  p <- ggplot(df, aes(x=value, y=sample_id)) +
    geom_col(color='black', fill='grey90') +
    geom_text(aes(label=scales::comma(value)),
              hjust = 1, nudge_x = -0.001) +
    xlab("Mean\nsilhouette width\n(Louvain res=0.2)") + ylab("") +
    scale_x_continuous(position='top',
                       expand = c(0,0)) +
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=10, color="black"),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  return(p)
}


#######
# Example
#######

# Set up metadata table for samples
# Two samples for testing
sample_meta <- data.frame(
  sample_id = c("Xen_BCa_4066", "CosMx_BCa_4066"),
  path = c(
    "/Volumes/Charmander/spatial_data/new_data/Xenium_Breast_4066_Run1/",
    "/Volumes/Charmander/spatial_data/new_data/CosMx_Breast_4066_Run1/"
  ),
  platform = c("Xenium", "CosMx")
)

# Full directory of samples
#path <- "/Volumes/Charmander/spatial_data/new_data/"
#sample_meta <- data.frame(
#  sample_id = list.files(path),
#  path = paste0(path, list.files(path)),
#  platform = stringr::word(list.files(path), 1, 1, sep="_")
#)

# Load objects from table
obj_list <- list()
for(i in 1:nrow(sample_meta)){
  obj_list[[i]] <- readSpatial(sample_id = sample_meta[i, "sample_id"],
                               path = sample_meta[i, "path"],
                               platform = sample_meta[i, "platform"])
}

# Probe set size
panel_size <- data.frame(
  sample_id = sample_meta$sample_id,
  platform = sample_meta$platform,
  value = unlist(lapply(obj_list, nrow))
)

# Total cells
cell_count <- data.frame(
  sample_id = sample_meta$sample_id,
  platform = sample_meta$platform,
  value = unlist(lapply(obj_list, ncol))
)

# Transcripts per cell
tx_per_cell <- do.call(rbind, lapply(obj_list, getTxPerCell))

# Transcripts per um2
tx_per_um2 <- do.call(rbind, lapply(obj_list, getTxPerArea))

# Transcripts per nucleus
tx_per_nuc <- do.call(rbind, lapply(obj_list, getTxPerNuc))

# Transcript per cell per 100 genes in panel
tx_per_cell_norm <- tx_per_cell
tx_per_cell_norm$value <- round(tx_per_cell_norm$value / panel_size$value,
                                digits = 2)

# Transcripts per cell with intesecting genes
common_genes <- Reduce(intersect, lapply(obj_list, rownames))
tx_per_cell_intersect <- do.call(rbind, lapply(obj_list, getTxPerCell, features=common_genes))

# Fraction transcripts in cells
tx_fraction_in_cell <- do.call(rbind, lapply(obj_list, getCellTxFraction))

# Mean signal-noise ratio
signal_ratio <- do.call(rbind, lapply(obj_list, getMeanSignalRatio))

# Expression
mean_expression <- do.call(rbind, lapply(obj_list, getMeanExpression))

# Max expression distribution
max_exp <- do.call(rbind, lapply(obj_list, getMaxDetection))

# MECR
mecr <- do.call(rbind, lapply(obj_list, getMECR))

# Moran's I
morans <- do.call(rbind, lapply(obj_list, getMorans))

# Silhouette Width
silhouette <- do.call(rbind, lapply(obj_list, getSilhouetteWidth))

# PLOT
p0 <- plotSampleLabel(sample_meta)
p1 <- plotPanelSize(panel_size)
p2 <- plotCellCount(cell_count)
p3 <- plotTxPerCell(tx_per_cell)
p4 <- plotTxPerArea(tx_per_um2) 
p5 <- plotTxPerNuc(tx_per_nuc)
p6 <- plotTxPerCellNorm(tx_per_cell_norm) 
p65 <- plotTxPerCell(tx_per_cell_intersect) + xlab("Per cell\ncommon\ngenes")
p7 <- plotFractionTxInCell(tx_fraction_in_cell)
p8 <- plotSignalRatio(signal_ratio) 
p9 <- plotMeanExpression(mean_expression)
p10 <- plotMECR(mecr)
p11 <- plotMorans(morans)
p12 <- plotSilhouette(silhouette)

p <- cowplot::plot_grid(p0, p1, p2, p3, p4, p5, p6, p65, p7, p8, p9, p10, p11, p12,
                        ncol=14, align='h',
                        rel_widths = c(0.75, 0.4, 0.5, 0.3, 0.3, 0.3, 0.3, 0.3, 0.75, 0.75, 0.8, 0.3, 0.8, 0.8))
cowplot::save_plot(p, filename="~/Downloads/benchmark_qc_small.pdf",
                   base_width=20, base_height=1.8)
