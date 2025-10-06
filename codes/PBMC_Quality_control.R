listOfPackages <- c("magrittr", "dplyr", "grid", "biomaRt", "ExperimentHub","ggrepel", "RColorBrewer",
                    "png", "plotly", "corpcor", "gplots", "VennDiagram", "SingleCellExperiment", 
                    # "scRNAseq",  "TabulaMurisData", "celldex", "singleCellNet",
                    "ggdendro" , "stringi", "sva", "scuttle", "EnsDb.Hsapiens.v86", "SeuratDisk",
                    "MASS", "data.table", "venn", "readxl", "tidyr","tibble", "writexl", "harmony", 
                    "edgeR", "ggplot2", "factoextra", "futile.logger", "psych", "Matrix", "scDblFinder",
                    "homologene", "tidyverse", "orca", "scales",  "ensembldb", "DESeq2", "SingleR",
                    "Seurat", "SeuratObject", "patchwork", "cowplot", "AnnotationHub")
new.packages <- listOfPackages[!(listOfPackages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,repos = "http://cran.us.r-project.org")
lapply(listOfPackages, require,character.only = TRUE)


source("Cryo_PBMC/Matrix.utils.R")

rm(list = ls(all = TRUE))
setwd("analysis/")

projectpath <- "data/"

sample <- c("PBMC1_fresh", "PBMC2_fresh", "PBMC3_fresh", "PBMC1_6M", "PBMC2_6M", "PBMC3_6M", "PBMC1_12M", "PBMC2_12M", "PBMC3_12M")

### initial QC for GEX
for (i in seq_along(sample)) {
  dataname <- sample[i]
  file <- paste0(projectpath, sample[i], "filtered_count/")
  assign(dataname, Read10X(data.dir = file))
  a <- CreateSeuratObject(counts = get(dataname), project = sample[i],  min.cells = 3, min.features = 200)
  
  # Mouse_mitochondrial_gene start with "^mt-", Human_mitochondrial_genes start with "^MT-"
  a[["percent.mt"]] <- PercentageFeatureSet(a, pattern = "^MT-")
  a$sample <- sample[i]
  a$group <- gsub("PBMC[1-3]_", "", sample[i])
  assign(sample[i], a)
}

# save(PBMC1_fresh, PBMC2_fresh, PBMC3_fresh, PBMC1_6M, PBMC2_6M, PBMC3_6M, "PBMC1_12M", "PBMC2_12M", "PBMC3_12M", file = "SampleObjects.RData")
load("/data/dcr_sp/DATA/Tom_PBMC/SampleObjects.RData")

merged_seurat <- merge(x = PBMC1_fresh, 
                       y = list(PBMC2_fresh, PBMC3_fresh, PBMC1_6M, PBMC2_6M, PBMC3_6M,  "PBMC1_12M", "PBMC2_12M", "PBMC3_12M"),
                       add.cell.id = sample)

# save(merged_seurat, file = "merged_seurat.RData")
load("analysis/merged_seurat.RData")
merged_seurat$Sample <- merged_seurat$orig.ident
head(merged_seurat@meta.data)
tail(merged_seurat@meta.data)

####QC
## Quality metrics
metadata <- merged_seurat@meta.data

# cell counts: the number of cell counts per sample
metadata %>% 
  ggplot(aes(x=Sample, fill=Sample)) + 
  ggplot(aes(x=Sample, fill=Sample)) + 
  #geom_bar(stat = "identity", position=position_dodge()) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  theme(legend.position = "none") +
  scale_x_discrete(limits = c("PBMC1_fresh", "PBMC1_6M", "PBMC1_12M", "PBMC2_fresh", "PBMC2_6M", "PBMC2_12M", "PBMC3_fresh", "PBMC3_6M", "PBMC3_12M")) +
  ggtitle("NCells")

metadata <- data.frame(sample = names(table(singlet_harmonized_PBMC$sample)))
metadata$filtered <- table(singlet_harmonized_PBMC$sample)[match(metadata$sample, names(table(singlet_harmonized_PBMC$sample)))]
metadata$original <- c( 8522, 11010, 12029, 9500, 11564, 11138, 6982, 8803, 13597)
metadata$Days <- str_replace(metadata$sample, "PBMC[1-3]_", "")
metadata$Days <- factor(metadata$Days, c("fresh", "6M", "12M"))
metadata$group <- gsub("_.*", "",  metadata$sample)
metadata$percentage <- round(metadata$filtered / metadata$original * 100, 2)
metadata <- metadata[match(c("PBMC1_fresh", "PBMC1_6M", "PBMC1_12M","PBMC2_fresh", "PBMC2_6M", "PBMC2_12M", "PBMC3_fresh", "PBMC3_6M", "PBMC3_12M"), metadata$sample),]
#metadata %>% dplyr::select(c(sample, group)) %>%
#  group_by(sample) %>%
#  mutate(Ncell =n(), 
#         Days = str_replace(sample, "PBMC[1-3]_", ""),
#         Days = factor(Days, c("fresh", "6M", "12M")),
#         group = gsub("_.*", "",  sample)) %>%
#  distinct() %>% 
metadata |>
  ggplot(aes(x=group, y = filtered, fill=Days)) + 
  geom_bar(stat = "identity", position=position_dodge()) +
  # geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
        axis.text.y = element_text(size = 10),
        axis.title  = element_text(size = 14),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  # theme(legend.position = "none") +
  # scale_x_discrete(limits = c("PBMC1_fresh", "PBMC1_6M", "PBMC1_12M","PBMC2_fresh", "PBMC2_6M", "PBMC2_12M", "PBMC3_fresh", "PBMC3_6M", "PBMC3_12M")) +
  ggtitle("Number of Cells")


## UMI counts (transcripts) per cell
# Visualize the number UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(color=Sample, x=nCount_RNA, fill= Sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  xlab("nUMI") +
  geom_vline(xintercept = 500)+
  geom_vline(xintercept = 200)

# Visualize the distribution of nUMI detected per cell via boxplot
metadata %>% 
  ggplot(aes(x=Sample, y=nCount_RNA, fill=Sample)) + 
  geom_violin() +
  geom_boxplot(width= 0.1, color="black", alpha=0.2, outlier.color = "red")+
  theme_classic() +
  scale_y_log10() + 
  ylab("nUMI") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  geom_hline(yintercept = 500)+
  geom_hline(yintercept = 200) +
  ggtitle("nCells vs nUMI")

## Genes detected per cell
#  Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color=Sample, x=nFeature_RNA, fill= Sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 200)+
  geom_vline(xintercept =  2500) +
  xlab("nGenes")

# Visualize the distribution of genes detected per cell via boxplot
metadata %>% 
  ggplot(aes(x=Sample, y=nFeature_RNA, fill=Sample)) + 
  geom_violin() +
  geom_boxplot(width= 0.1, color="black", alpha=0.2, outlier.color = "red")+
  theme_classic() +
  scale_y_log10() + 
  ylab("nGenes") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  geom_hline(yintercept = 300)+
  geom_hline(yintercept =  600) +
  geom_hline(yintercept =  1000) +
  ggtitle("nGenes per cell")

## UMIs vs. genes detected
# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  xlab("nUMI") +
  ylab("nGenes")+ 
  theme_classic() +
  geom_hline(yintercept = 200) +
  geom_hline(yintercept = 2500) +
  facet_wrap(~Sample)

## Mitochondrial counts ratio
# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=Sample, x=percent.mt, fill=Sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  xlab("mitoPercentile") +
  geom_vline(xintercept = 5) +
  geom_vline(xintercept = 20)

# Visualize the distribution of percentile of mitochondrial genes counts per cell via boxplot
metadata %>% 
  ggplot(aes(x=Sample, y=percent.mt, fill=Sample)) + 
  geom_violin() +
  geom_boxplot(width= 0.1, color="black", alpha=0.2, outlier.color = "red")+
  theme_classic() +
  scale_y_log10() + 
  ylab("mitoPercentile") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  geom_hline(yintercept = 10)+
#   geom_hline(yintercept =  20) +
  ggtitle("mitoPercentile per cell")


#### batch effect detection and correction
## cell cycle regression
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat as cc.genes.
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

convertHumanGeneList <- function(x){
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "www.ensembl.org")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "www.ensembl.org")
  genesV2 = getLDS(attributes = c("hgnc_symbol"),
                   filters = "hgnc_symbol",
                   values = x , 
                   mart = human, 
                   attributesL = c("mgi_symbol"), 
                   martL = mouse, 
                   uniqueRows=T)
  mousex <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(mousex))
  return(mousex)
}

# ms_genes <- homologene(s.genes, inTax = 9606, outTax = 10090)

mg2m_genes<- homologene(g2m.genes, inTax = 9606, outTax = 10090) 

merged_seurat  <- CellCycleScoring(merged_seurat, 
                                    g2m.features = g2m.genes,
                                    s.features = s.genes)

# save(merged_seurat, file = "merged_seurat.RData")
load("merged_seurat")

# normalize the counts
merged_seurat <- NormalizeData(merged_seurat, verbose =FALSE)

unique(merged_seurat$orig.ident)

# Identify the most variable genes
merged_seurat <- FindVariableFeatures(merged_seurat, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)
# Scale the counts
merged_seurat <- ScaleData(merged_seurat)

# Perform PCA
merged_seurat <- RunPCA(merged_seurat)
ElbowPlot(merged_seurat)

# Plot the PCA colored by cell cycle phase (batch effects)
DimPlot(merged_seurat,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")

# finding other batch effects
merged_seurat <- FindNeighbors(object = merged_seurat, dims= 1:20)
merged_seurat <- FindClusters(object = merged_seurat)
merged_seurat <- RunUMAP(object = merged_seurat, dims= 1:20)

#Plot
DimPlot(merged_seurat, reduction = "umap", group.by = "Sample")


#### prepare for integrate seurat object
# filter raw GEX data and merge VDJ data into GEX data with the determined QC metrics
for (i in seq_along(sample)) {
  dataname <- sample[i]
  file <- paste0(projectpath, sample[i], "/filtered_count/")
  assign(dataname, Read10X(data.dir = file))
  a <- CreateSeuratObject(counts = get(dataname), project = sample[i],  min.cells = 3, min.features = 200)
  
  # normalizing
  # Mouse_mitochondrial_gene start with "^mt-", Human_mitochondrial_genes start with "^MT-"
  a[["percent.mt"]] <- PercentageFeatureSet(a, pattern = "^MT-")
  a$sample <- sample[i]
  a$group <- gsub("PBMC[1-3]_", "", sample[i])
  
  ## filter cell on gene number per cell 
  a <- subset(a, subset = nFeature_RNA > 600  & percent.mt < 10)
  a <- NormalizeData(a, verbose = FALSE)
  a <- FindVariableFeatures(a, selection.method = "vst", nfeatures = 2000)
  
  # add cell_cycle score to Metadata
  a <- CellCycleScoring(a, 
                         g2m.features = g2m.genes, 
                         s.features = s.genes)
  assign(sample[i], a)
}

#### Integration of seurat object
## standard approach to integrate seurat object
## https://satijalab.org/seurat/articles/integration_introduction.html

#### CCA integration:
combined_data<- list(PBMC1_fresh = PBMC1_fresh,
                     PBMC2_fresh = PBMC2_fresh, 
                     PBMC3_fresh = PBMC3_fresh,
                     PBMC1_6M = PBMC1_6M,
                     PBMC2_6M = PBMC2_6M,
                     PBMC3_6M = PBMC3_6M,
                     PBMC1_12M = PBMC1_12M,
                     PBMC2_12M = PBMC2_12M,
                     PBMC3_12M = PBMC3_12M)

merged_seurat <- merge(x = combined_data[[1]], 
                       y = combined_data[-1],
                       add.cell.id = sample)

merged_seurat$Sample <- NULL
head(merged_seurat@meta.data)
tail(merged_seurat@meta.data)
save(merged_seurat, file = "merged_seurat.RData")

# select features that are repeatedly variable across datasets for integration
# https://satijalab.org/seurat/articles/integration_large_datasets.html
# preparing SCT normalization. 
load("merged_seurat.RData")
split_seurat <- vector(mode = 'list', length = 3)
options(future.globals.maxSize = 8000 * 1024^2)
for (i in 1:length(combined_data)) {
  split_seurat[[i]] <- SCTransform(combined_data[[i]], vars.to.regress = c("percent.mt"), vst.flavor = "v2")
}

features <- SelectIntegrationFeatures(object.list = split_seurat, nfeatures = 3000)
split_seurat <- PrepSCTIntegration(object.list = split_seurat, anchor.features = features)
  anchors30 <- FindIntegrationAnchors(object.list = split_seurat, 
                                              normalization.method = "SCT",
                                              anchor.features = features, 
                                              dim = 1:30) #save this object
#  save(anchors30, file = "anchors30.RData")}
combined30       <- IntegrateData(anchorset = anchors30, normalization.method = "SCT", dim = 1:30) #save this object
# save(combined30, file = "Combined_cca30.RData")
load("Combined_cca30.RData")

# specify that we will perform downstream analysis on the corrected data.
# note that the original unmodified data still resides in the 'RNA' assay
# Now we can run a single integrated analysis on all cells!
DefaultAssay(combined30) <- "integrated"

for (i in  1:1) {
  combined_scaled <- ScaleData(combined30, verbose = FALSE)
  combined_PCA    <- RunPCA(combined_scaled, npcs = 30, verbose = FALSE, features = VariableFeatures(object = combined30))
  DimHeatmap(combined_PCA, dims = 1, cells = 500, balanced = TRUE, fast=FALSE)
  
  #  combined_dim <- JackStraw(combined_PCA, num.replicate = 100)
  #  combined_dim <- ScoreJackStraw(combined_dim)  # save this object
  #  ElbowPlot(combined_dim, ndims=30)
  
  # Determine the K-nearest neighbor graph
  combined_cluster <- RunUMAP(combined_PCA, reduction = "pca", dims = 1:30)
  combined_cluster <- FindNeighbors(combined_cluster, reduction = "pca", dims = 1:30)
  
  # Determine the clusters for various resolutions                                
  combined_cluster <- FindClusters(object = combined_cluster,
                                   resolution = c(0.1, 0.2, 0.3, 0.4, 0.6, 0.8)) # save this object
  DefaultAssay(combined_cluster) <- "RNA"
  save(combined_cluster, file = "combined_cluster_cca30.RData")
}

# Visualization of clusters
DimPlot(combined_cluster, 
        reduction = "umap", 
        label =TRUE, 
        repel = TRUE,
        ncol = 7,
        split.by = "orig.ident")

#### modified approach to integrate large seurat object: rpca integration
# ref: https://satijalab.org/seurat/articles/integration_introduction.html
# https://satijalab.org/seurat/articles/integration_large_datasets.html
# https://satijalab.org/seurat/articles/integration_rpca.html


# preparing SCT normalization. 
#split_seurat <- vector(mode = 'list', length = 3)
options(future.globals.maxSize = 8000 * 1024^2)
#for (i in 1:length(combined_data)) {
#  split_seurat[[i]] <- SCTransform(combined_data[[i]], vars.to.regress = c("percent.mt"), vst.flavor = "v2")
#}
# or 
for (i in 1:1) {
split_seurat <- SplitObject(merged_seurat, split.by = "sample")
split_seurat <- lapply(X = split_seurat, FUN = SCTransform, vars.to.regress = c("percent.mt"), vst.flavor = "v2")

features <- SelectIntegrationFeatures(object.list = split_seurat, nfeatures = 3000)
split_seurat <- PrepSCTIntegration(object.list = split_seurat, anchor.features = features)
split_seurat <- lapply(X = split_seurat, FUN = RunPCA, features = features)
  
  rpca_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                         # reference = c(1), 
                                         normalization.method = "SCT",
                                         anchor.features = features,
                                         reduction = "rpca", 
                                         dims = 1:30,
                                         k.anchor = 20)
  # save(rpca_anchors, file = "rpca_anchors30.RData")
  combined_rpca <- IntegrateData(anchorset = rpca_anchors, dims = 1:30)
  DefaultAssay(combined_rpca) <- "integrated"
save(combined_rpca, file = "combined_rpca30.RData")
}


# specify that we will perform downstream analysis on the corrected data.
# note that the original unmodified data still resides in the 'RNA' assay
# Now we can run a single integrated analysis on all cells!


for (i in  1:1) {
  combined_scaled <- ScaleData(combined_rpca, verbose = FALSE)
  combined_PCA    <- RunPCA(combined_scaled, npcs = 30, verbose = FALSE, features = VariableFeatures(object = combined_rpca))
  DimHeatmap(combined_PCA, dims = 1, cells = 500, balanced = TRUE, fast=FALSE)
  
  #  combined_dim <- JackStraw(combined_PCA, num.replicate = 100)
  #  combined_dim <- ScoreJackStraw(combined_dim)  # save this object
  #  ElbowPlot(combined_dim, ndims=30)
  
  # Determine the K-nearest neighbor graph
  combined_cluster <- RunUMAP(combined_PCA, reduction = "pca", dims = 1:30)
  combined_cluster <- FindNeighbors(combined_cluster, reduction = "pca", dims = 1:30)
  
  # Determine the clusters for various resolutions                                
  combined_cluster <- FindClusters(object = combined_cluster,
                                   resolution = c(0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0)) # save this object
  DefaultAssay(combined_cluster) <- "RNA"
  save(combined_cluster, file = "combined_cluster_rpca30.RData")
}

# Visualization of clusters
DimPlot(combined_cluster, 
        reduction = "umap", 
        label =TRUE, 
        repel = TRUE,
        ncol = 3,
        split.by = "orig.ident")

#### Harmony integration
# ref: https://hbctraining.github.io/scRNA-seq_online/lessons/06a_integration_harmony.html
# preparing SCT normalization. 
split_seurat <- vector(mode = 'list', length = 3)
options(future.globals.maxSize = 8000 * 1024^2)
for (i in 1:length(combined_data)) {
  split_seurat[[i]] <- SCTransform(combined_data[[i]], vars.to.regress = c("percent.mt"), vst.flavor = "v2")
}
# or 
split_seurat <- lapply(X = combined_data, FUN = SCTransform, vars.to.regress = c("percent.mt"), vst.flavor = "v2")

# Find most variable features across samples to integrate
features <- SelectIntegrationFeatures(object.list = split_seurat, nfeatures = 3000)

# Merge normalized samples
merged_seurat <- merge(x = split_seurat[[1]],
                       y = split_seurat[2:length(split_seurat)],
                       merge.data = TRUE)

DefaultAssay(merged_seurat) <- "SCT"

# Manually set variable features of merged Seurat object
VariableFeatures(merged_seurat) <- features

# Calculate PCs using manually set variable features
merged_seurat <- RunPCA(merged_seurat, assay = "SCT", npcs = 50)

# Harmony integration
harmonized_seurat <- RunHarmony(merged_seurat, 
                                group.by.vars = "orig.ident", 
                                reduction = "pca",
                                assay.use = "SCT",
                                reduction.save = "harmony")

save(harmonized_seurat, file = "combined_harmony.RData")


# UMAP and clustering with top PCs
harmonized_seurat <- RunUMAP(harmonized_seurat, reduction = "harmony", assay = "SCT", dims = 1:30)
harmonized_seurat <- FindNeighbors(object = harmonized_seurat, reduction = "harmony")
harmonized_seurat <- FindClusters(harmonized_seurat, resolution = c(0.1, 0.2, 0.3, 0.4, 0.6, 0.8))

save(harmonized_seurat, file = "combined_cluster_harmony.RData")
load("combined_cluster_harmony.RData")

# Visualization of clusters
DimPlot(harmonized_seurat, 
        reduction = "umap", 
        label =TRUE, 
        repel = TRUE,
        ncol = 3,
        split.by = "orig.ident")



## predict celltype by Seurat mapping
# http://satijalab.org/seurat/articles/multimodal_reference_mapping.html
PBMC_seurat_ref <- LoadH5Seurat("pbmc_multimodal.h5seurat")
DimPlot(object = PBMC_seurat_ref, reduction = "wnn.umap",
        group.by = "celltype.l1",
        label = TRUE, label.size = 3, repel = TRUE) + NoLegend() |
DimPlot(object = PBMC_seurat_ref, reduction = "wnn.umap",
        group.by = "celltype.l2",
        label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
anchors <- FindTransferAnchors(reference = PBMC_seurat_ref,
                               query = combined_cluster,
                               normalization.method = "SCT",
                               reference.reduction = "spca",
                               dims = 1:50)

combined_cluster <- MapQuery(anchorset = anchors,
                              query = combined_cluster,
                              reference = PBMC_seurat_ref,
                              refdata = list(
                                PBMC_celltype.l1 = "celltype.l1",
                                PBMC_celltype.l2 = "celltype.l2"
                                # predicted_ADT = "ADT"
                              ),
                              reference.reduction = "spca",
                              reduction.model = "wnn.umap"
                              )

save(combined_cluster, file = "combined_cluster_rpca30.RData")

DimPlot(combined_cluster, reduction = "ref.umap", 
        group.by = "singler_PBMC_label", 
        label =TRUE, label.size = 3, repel =TRUE) + NoLegend() |
DimPlot(combined_cluster, reduction = "ref.umap", 
        group.by = "predicted.PBMC_celltype.l1", 
        label =TRUE, label.size = 3, repel =TRUE) + NoLegend() |
  DimPlot(combined_cluster, reduction = "ref.umap", 
          group.by = "predicted.PBMC_celltype.l2", 
          label =TRUE, label.size = 3, repel =TRUE) + NoLegend()


#### doublet/multiplet prediction and correction by scDblFinder 
## ref:https://github.com/plger/scDblFinder/blob/devel/vignettes/scDblFinder.Rmd
## ref:https://bioconductor.org/packages/devel/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html#:~:text=By%20default%2C%20scDblFinder%20will%20generate,lead%20to%20lower%20detection%20accuracy.
# load("combined_cluster_rpca30.RData")
sce <- as.SingleCellExperiment(combined_cluster)
table(sce$predicted.PBMC_celltype.l1, sce$sample)
sce <- scDblFinder(sce, clusters = "predicted.PBMC_celltype.l1", samples = "sample")
table(sce$scDblFinder.class)
dim(sce@colData); dim(combined_cluster@meta.data)

if (all(rownames(combined_cluster@meta.data) == rownames(sce@colData))) {
  combined_cluster$doublet_Class = sce$scDblFinder.class
  
  combined_cluster$condition <- gsub("PBMC[1-9]_", "", combined_cluster$sample)
  
  data.frame(a = combined_cluster$Sample, b = combined_cluster$condition, c = combined_cluster$Days) %>% distinct()
  
  DefaultAssay(combined_cluster) <- "RNA"
  save(combined_cluster, file = "combined_cluster_rpca30.RData")
}

DimPlot(combined_cluster,
        reduction = "umap",
        group.by = "predicted.PBMC_celltype.l1", 
        split.by = "doublet_Class",
        ncol = 2,
        label = TRUE,
        repel = TRUE,
        label.size = 4)

singlet_combined_cluster <- subset(combined_cluster, subset = doublet_Class == "singlet")
# save(singlet_combined_cluster, file = "batch123_singlet_cluster.RData")
load("singlet_harmonized_PBMC.RData")
DimPlot(singlet_harmonized_PBMC,
        reduction = "ref.umap",
        group.by = "predicted.PBMC_celltype.l2", 
        # split.by = "doublet_Class",
        # ncol = 2,
        label = TRUE,
        repel = TRUE,
        label.size = 4) & ggtitle("Cell Type")

#### IL-18 Expression level 
singlet_harmonized_PBMC$celltype <- singlet_harmonized_PBMC$predicted.PBMC_celltype.l1
singlet_harmonized_PBMC$celltype[singlet_harmonized_PBMC$predicted.PBMC_celltype.l2 == "gdT"] <- "gdT"
VlnPlot(singlet_harmonized_PBMC, c("IL18"), sort = TRUE, group.by = "celltype")
VlnPlot(singlet_harmonized_PBMC, c("CD28", "CCR5", "CXCR4"), sort = TRUE, 
        ncol = 1, pt.size = 0, group.by = "celltype")
singlet_harmonized_PBMC$celltype <- factor(singlet_harmonized_PBMC$celltype, 
                                           levels = c("Mono", "DC", "CD8 T", "NK", "CD4 T", "other", "other T", "gdT", "B"))
VlnPlot(singlet_harmonized_PBMC, c("IL10RA", "IL10RB", "IFNLR1"), 
        ncol = 1, pt.size = 0, group.by = "celltype")

#### Toll-Like Receptors (TLRs) Expression level 
VlnPlot(singlet_harmonized_PBMC, paste0("TLR", 1:10),
        sort = TRUE,
        ncol = 3,
        pt.size = 0,
        group.by = "celltype")



### immune cell-type markers
#https://www.jimmunol.org/content/208/2/396/tab-figures-data
#https://panglaodb.se/markers.html?cell_type=%27Stromal%20cells%27

# convert human genesymobl to mouse genesymbol:homologene package
homologene::taxData
homologene(c(""), 
inTax = 10090, outTax = 9606)

genelist <- list(CD4_T_cell = c("CD4", "IL7R", "CCR7", "CD3D"), 
                 CD8_T_cell = c("CD8A", "CD8B", "CD3D"),
                 B_cell = c("MS4A1", "CD79A", "CD79B", "CD19"),
                 NK_cell = c("NKG7", "GNLY", "KLRB1"),
                 Monocytes = c("CD14", "LYZ", "FCGR3A", "MS4A7"),
                 Dendritic_cells = c("FCER1A", "CST3")
)



## heatmap plot for top-upregulated/marker genes
genelist <-c(sig_res_12M$gene, sig_res_6M$gene) %>% unique() %>% sort()
Ave_expr <- Seurat::AverageExpression(object = singlet_harmonized_PBMC, layer = "counts", return.seurat = TRUE,
                                             group.by = "predicted.PBMC_celltype.l1", feature = unlist(genelist))

singlet_harmonized_PBMC$type = paste0(singlet_harmonized_PBMC$predicted.PBMC_celltype.l1, "_", singlet_harmonized_PBMC$group)
unique(singlet_harmonized_PBMC$type)

Ave_expr <- Seurat::AverageExpression(object = subset(singlet_harmonized_PBMC, subset = group == "fresh"),
                                      layer = "counts", return.seurat = TRUE,
                                      group.by = "predicted.PBMC_celltype.l1", feature = unlist(genelist))

Ave_expr <- Seurat::AverageExpression(object = singlet_harmonized_PBMC, layer = "counts", return.seurat = TRUE,
                                      group.by = "type", feature = unlist(genelist))

Ave_expr@meta.data$orig.ident

Ave_expr@active.ident <- factor(rownames(Ave_expr@meta.data),
                                # levels= c("CD4 T", "CD8 T", "B", "NK", "Mono", "DC", "other T", "other")
                                levels= c("CD4 T_fresh", "CD4 T_6M", "CD4 T_12M", "CD8 T_fresh", "CD8 T_6M", "CD8 T_12M",
                                          "B_fresh", "B_6M", "B_12M", "NK_fresh", "NK_6M", "NK_12M",
                                          "Mono_fresh", "Mono_6M", "Mono_12M", "DC_fresh", "DC_6M", "DC_12M",
                                          "other T_fresh", "other T_6M", "other T_12M", "other_fresh", "other_6M", "other_12M"))

Ave_expr@active.ident <- factor(Ave_expr$orig.ident, levels= c("CD4 T", "CD8 T", "B", "NK", "Mono", "DC", "other T", "other"))

dim(Ave_expr@assays$RNA@counts)
 
DoHeatmap(Ave_expr, draw.lines = TRUE, lines.width = 1, features = unlist(genelist), label=T, size = 3) + 
  guides(color="none") + 
  # ggtitle("fresh") +
  theme(plot.title = element_text(color="red", size=14, face="bold.italic", hjust = 0.5)) 

## ref: https://github.com/satijalab/seurat/issues/5629
cell.order <- c("CD4 T_fresh", "CD4 T_6M", "CD4 T_12M", "CD8 T_fresh", "CD8 T_6M", "CD8 T_12M",
                "B_fresh", "B_6M", "B_12M", "NK_fresh", "NK_6M", "NK_12M",
                "Mono_fresh", "Mono_6M", "Mono_12M", "DC_fresh", "DC_6M", "DC_12M",
                "other T_fresh", "other T_6M", "other T_12M", "other_fresh", "other_6M", "other_12M")

mat<- Ave_expr[["RNA"]]@data %>% as.matrix()
## scale the rows
mat<- t(scale(t(mat)))

cluster_anno<- Ave_expr@meta.data$orig.ident

## make the black color map to 0. the yellow map to highest and the purle map to the lowest
col_fun = circlize::colorRamp2(c(-1, 0, 3), c("#FF00FF", "black", "#FFFF00"))
sample_order<-  cell.order
# plot the heatmap
meheatmap1<-ComplexHeatmap::Heatmap(mat, name = "Expression",  
                    column_split = factor(cluster_anno, levels= c("CD4 T", "CD8 T", "B", "NK", "Mono", "DC", "other T", "other")),
                    cluster_columns = F,
                    show_column_dend = FALSE,
                    cluster_column_slices = F,
                    column_title_gp = gpar(fontsize = 10),
                    column_gap = unit(0.5, "mm"),
                    cluster_rows = F,
                    show_row_dend = FALSE,
                    col = col_fun,
                    row_names_gp = gpar(fontsize = 10),
                    column_title_rot = 90,
                    top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scales::hue_pal()(9)))),
                    show_column_names = T,
                    use_raster = TRUE,
                    raster_quality = 4,column_order = sample_order,row_names_side = "left")

draw(meheatmap1)

#### cell population change over time
a <- singlet_combined_cluster@meta.data

a <- a %>%
  group_by(predicted.PBMC_celltype.l1, orig.ident) %>%
  dplyr::summarise(cell.count = n()) %>%
  spread(., orig.ident, cell.count)
a[is.na(a)] <- 0

cell_counts <- a[, 1]
for (i in seq_len(ncol(a) - 1)) {
  cell_counts <- cbind(cell_counts, round(100 *a[, (i+1)]/sum(a[, (i+1)]), 2))
  colnames(cell_counts)[i + 1] <- colnames(a)[i+1]                    
}

cell_population <- a %>% gather(sample, value, -predicted.PBMC_celltype.l1)
# cell_population <- cell_counts %>% gather(sample, value, -predicted.PBMC_celltype.l1)
cell_population$sample <- factor(cell_population$sample, levels = c("PBMC1_fresh", "PBMC1_6M", "PBMC1_12M", "PBMC2_fresh", "PBMC2_6M", "PBMC2_12M","PBMC3_fresh", "PBMC3_6M", "PBMC3_12M"))
cell_population$group <- gsub("_.*", "", cell_population$sample)
cell_population$percentage <- cell_population$value/table(CD4T$sample)[match(cell_population$sample, names(table(CD4T$sample)))] * 100

cell_population %>%
  ggplot(aes(x = group, y = percentage)) +
  geom_col(aes(fill = sample), width = 0.7, position=position_dodge()) +
  xlab("CXCR4+/CCR5+ CD4T") +
  ylab("Percentage") +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+
 # facet_grid(~ group)


#### pseudo DE analysis for all cell with DESeq2
# ref: https://hbctraining.github.io/scRNA-seq_online/lessons/pseudobulk_DESeq2_scrnaseq.html
# Extract raw counts and metadata to create SingleCellExperiment object
table(singlet_harmonized_PBMC$group, singlet_harmonized_PBMC$sample)

counts <- singlet_harmonized_PBMC@assays$RNA@counts
# counts <- CD4_CCR5_Mono@assays$RNA@counts
head(counts)

metadata <- singlet_harmonized_PBMC@meta.data
# metadata <- CD4_CCR5_Mono@meta.data
dim(metadata)
all(colnames(counts) == rownames(metadata))

# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- singlet_harmonized_PBMC$sample
# metadata$cluster_id <- CD4_CCR5_Mono$sample
unique(metadata$cluster_id)

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), colData = metadata)

## Check the assays present
assays(sce)

## Check the counts matrix
dim(counts(sce))
counts(sce)[1:6, 1:6]

# Explore the cellular metadata for the dataset
dim(colData(sce))
head(colData(sce))

# Subset metadata to include only the variables you want to aggregate across (here, we want to aggregate by sample and by cluster)
groups <- colData(sce)[, c("sample", "predicted.PBMC_celltype.l1")]
dim(colData(sce))
unique(groups)

# Aggregate across cluster-sample groups
# transposing row/columns to have cell_ids as row names matching those of groups
dim(t(counts(sce)))

aggr_counts <- aggregate.Matrix(t(counts(sce)), 
                                groupings = groups,
                                fun = "sum") 

# Explore output matrix
class(aggr_counts)
dim(aggr_counts)
aggr_counts[1:6, 1:6]

# Splitting the counts matrix by condition
# Transpose aggregated matrix to have genes as rows and samples as columns
aggr_counts <- t(aggr_counts)
dim(aggr_counts)
aggr_counts[1:6, 1:6]

# batch <- c(1, 2, 1, 2, 1, 2, rep(1, 2), 2, rep(1, 3), rep(2, 3), rep(1, 3), rep(2, 2), rep(1, 3), rep(2 ,3), rep(1, 3), rep(2, 3))
# adjusted_aggr_counts <- ComBat_seq(as.matrix(aggr_counts), batch=batch, group=NULL)
# aggr_counts <- adjusted_aggr_counts
# aggr_counts[1:6, 1:6]

# As a reminder, we stored our cell types in a vector called cluster_names
cluster_names <- unique(colData(sce)[, "predicted.PBMC_celltype.l1"])

# Loop over all cell types to extract corresponding counts, and store information in a list
## Initiate empty list
counts_ls <- list()

for (i in 1:length(cluster_names)) {
  
  ## Extract indexes of columns in the global matrix that match a given cluster
  column_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[3]] == cluster_names[i])
  
  ## Store corresponding sub-matrix as one element of a list
  counts_ls[[i]] <- aggr_counts[, column_idx]
  names(counts_ls)[i] <- cluster_names[i]
  
}

# Explore the different components of the list
str(counts_ls)

## Initiate empty list
metadata_ls <- list()

# Number of cells per sample and cluster
t <- table(colData(sce)$sample,
           colData(sce)$predicted.PBMC_celltype.l1)
t[1:6, 1:6]

for (i in 1:length(counts_ls)) {
  
  ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
  df <- data.frame(cluster_sample_id = colnames(counts_ls[[i]]))
  
  ## Use tstrsplit() to separate cluster (cell type) and sample IDs
  df$cluster_id <- tstrsplit(df$cluster_sample_id, "_")[[3]]
  df$condition <- tstrsplit(df$cluster_sample_id, "_")[[2]]
  df$sample_id  <- gsub("_[^_]+$", "", df$cluster_sample_id)
  
  ## Retrieve cell count information for this cluster from global cell count table
  idx <- which(colnames(t) == unique(df$cluster_id))
  cell_counts <- t[, idx]
  
  ## Remove samples with zero cell contributing to the cluster
  cell_counts <- cell_counts[cell_counts > 0]
  
  ## Append cell_counts to data frame
  df$cell_count <- cell_counts[match(df$sample_id, names(cell_counts))]
  
  ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
  rownames(df) <- df$cluster_sample_id
  
  ## Store complete metadata for cluster i in list
  metadata_ls[[i]] <- df
  names(metadata_ls)[i] <- unique(df$cluster_id)
}

# Explore the different components of the list
str(metadata_ls)

# Create directories to save results if they don't already exist:
setwd("/hpcdata/dcr_sp/DATA/Tom_PBMC_081723")
if (!dir.exists("NK")) { dir.create("NK") }

# Function to run DESeq2 Wald Test and get results for any cluster:
## clustx is the name of the cluster (cell type) on which to run the function
## A is the sample group to compare (e.g. stimulated condition)
## B is the sample group to compare against (base/control level)
## padj_cutoff defines the ajusted p-value cutoff for significance (set to 0.05 by default)

## This function assumes the counts matrices and metadata for all clusters have been prepared
## and arranged in matching named lists (as illustrated in tutorial above)
## This function assumes the contrast (e.g. stim vs. control) is stored in a variable named "group_id"

get_dds_resultsAvsB <- function(clustx, A, B, padj_cutoff = 0.05) {
  for (i in c("6M", "12M")) {
  clustx = "NK"
  A = i
  B = "fresh"
  print(clustx) # useful for debugging
  
  # Extract counts matrix and metadata for cluster x
  idx <- which(names(counts_ls) == clustx)
  cluster_counts <- counts_ls[[idx]]
  cluster_counts <- counts_ls[[idx]][,grepl(paste0(A, "|", B), colnames(counts_ls[[idx]]))]
  cluster_metadata <- metadata_ls[[idx]]
  cluster_metadata <- cluster_metadata[grepl(paste0(A, "|", B), rownames(cluster_metadata)),]
  cluster_metadata$condition <- factor(cluster_metadata$condition, levels = c(B, A))
  
  # Print error message if sample names do not match
  if ( all(colnames(cluster_counts) != rownames(cluster_metadata)) ) {
    print("ERROR: sample names in counts matrix columns and metadata rows do not match!")
  }
  
  dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                colData = cluster_metadata, 
                                design = ~ condition)
  
  # Transform counts for data visualization
  rld <- rlog(dds, blind = TRUE)
  
  # Generate QC plots
  ## Plot and save PCA plot
  DESeq2::plotPCA(rld, intgroup = "condition")
  #if (!dir.exists("results")) { dir.create("results") }
  #ggsave(paste0("results/", clustx, "_specific_PCAplot.png"))
  
  ## Extract rlog matrix from the object and compute pairwise correlation values
  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)
  
  ## Plot and save heatmap
  #png(paste0("results/", clustx, "_specific_heatmap.png"),
  #    height = 6, width = 7.5, units = "in", res = 300)
  pheatmap(rld_cor, annotation = cluster_metadata[, c("condition"), drop = FALSE])
  #dev.off()
  
  # Run DESeq2 differential expression analysis
  dds <- DESeq(dds)
  
  ## Plot dispersion estimates
  #png(paste0("results/", clustx, "_dispersion_plot.png"),
  #    height = 5, width = 6, units = "in", res = 300)
  #plotDispEsts(dds)
  #dev.off()
  
  ## Output and shrink results of Wald test for contrast A vs B
  contrast <- paste(c("condition", A, "vs", B), collapse = "_")
  resultsNames(dds)
  
  res <- results(dds, name = contrast, alpha = 0.05)
  res <- lfcShrink(dds, coef = contrast, res = res)
  
  ## Turn the results object into a tibble for use with tidyverse functions
  res_tbl <- res %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble()
  if ( A == "6M") {
    res_tbl_6M = res_tbl
  } else {
    res_tbl_12M = res_tbl
  }
  # res_tbl$ENSEMBL <- mapIds(org.Mm.eg.db, keys = res_tbl$gene, keytype = "SYMBOL", column = "ENSEMBL")
  
  #write.csv(res_tbl,
  #          paste0( clustx, "_", contrast, "_all_genes.csv"),
  #          quote = FALSE, 
  #          row.names = FALSE)
  
  ## Subset the significant results
  padj_cutoff <- 0.05
  sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
    dplyr::arrange(padj)
  
  if ( A == "6M") {
  sig_res_6M = sig_res
  } else {
  sig_res_12M = sig_res
  }
  
  write.csv(sig_res,
            paste0(clustx, "_", A, "_signif_genes.csv"),
            quote = FALSE, 
            row.names = FALSE)
  
  # Generate results visualization plots
  ## Extract normalized counts from dds object
  normalized_counts <- counts(dds, normalized = TRUE)
  if ( A == "6M") {
    normalized_counts_6M = normalized_counts
  } else {
    normalized_counts_12M = normalized_counts
  }
  
  ## Extract top 20 DEG from resLFC (make sure to order by padj)
  top50_sig_genes <- sig_res %>%
    dplyr::arrange(padj) %>%
    dplyr::pull(gene) %>%
    head(n = dim(sig_res)[1])
  
  ## Extract matching normalized count values from matrix
  top50_sig_counts <- normalized_counts[rownames(normalized_counts) %in% top50_sig_genes, ]
  
  ## Convert wide matrix to long data frame for ggplot2
  top50_sig_df <- data.frame(top50_sig_counts)
  top50_sig_df$gene <- rownames(top50_sig_counts)
  
  top50_sig_df <- melt(setDT(top50_sig_df), 
                       id.vars = c("gene"),
                       variable.name = "Sample") %>% 
    data.frame()
  
  ## Replace "." by " " in cluster_sample_id variable (melt() introduced the ".")
  top50_sig_df$cluster_sample_id <- gsub("\\.", " ", top50_sig_df$Sample)
  
  ## Join counts data frame with metadata
  top50_sig_df <- plyr::join(top50_sig_df, as.data.frame(colData(dds)),
                             by = "cluster_sample_id")
  
  # Heatmap of top20-up-/top10-down-regulatated genes
  ## Extract top 20/1- DEG from resLFC (make sure to order by padj)

    top_sig_res <- sig_res[sig_res$gene %in% top50_sig_genes, ] %>% dplyr::arrange(desc(log2FoldChange))
      
    top_sig_genes <- top_sig_res %>% dplyr::pull(gene)
    
    ## Extract normalized counts for significant genes only
    sig_counts <- normalized_counts[rownames(normalized_counts) %in% top_sig_genes, ]
    sig_counts <- sig_counts[match(top_sig_genes, rownames(sig_counts)), c(2, 4, 6, 1, 3, 5)]
    
    ## Set a color-blind friendly palette
    heat_colors <- rev(brewer.pal(11, "PuOr"))
 
    ## Run pheatmap using the metadata data frame for the annotation
    pheatmap::pheatmap(sig_counts, 
             color = heat_colors, 
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             show_rownames = TRUE,
             # annotation = cluster_metadata[, c("condition", "sample")], 
             border_color = NA, 
             fontsize = 10, 
             scale = "row", 
             fontsize_row = 10, 
             height = 20) 
  }
    ####
  for (i in c("6M", "12M")) {
    print(i)
  if ( i == "12M") {
    sig_counts <- normalized_counts_12M[rownames(normalized_counts_12M) %in% sort(union(sig_res_12M$gene, sig_res_6M$gene)), ]
  } else {
    sig_counts <- normalized_counts_6M[rownames(normalized_counts_6M) %in% sort(union(sig_res_12M$gene, sig_res_6M$gene)), ]
  }
  
    sig_counts <- sig_counts[order(rownames(sig_counts)), c(2, 1, 4, 3, 6, 5)]
    sig_counts <- data.frame(sig_counts)
    
    if ( i == "12M") {
    sig_counts$padj = sig_res_12M$padj[match(rownames(sig_counts), sig_res_12M$gene)]
    } else {
    sig_counts$padj = sig_res_6M$padj[match(rownames(sig_counts), sig_res_6M$gene)]
    }
    
    sig_counts$colors=ifelse(!is.na(sig_counts$padj) , "red", "green")
    
    
    ## Run pheatmap using the metadata data frame for the annotation
    e <- 
    pheatmap::pheatmap(sig_counts[, 1:6], 
             color = heat_colors, 
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             show_rownames = TRUE,
             # annotation = cluster_metadata[, c("condition", "sample")], 
             border_color = NA, 
             fontsize = 10, 
             scale = "row", 
             fontsize_row = 10, 
             height = 20)
    cols = sig_counts[order(match(rownames(sig_counts), e$gtable$grobs[[3]]$label)), ]$colors
    e$gtable$grobs[[3]]$gp=gpar(col=cols)
    print(e)
  }
  
  for (i in 1:1){
  #### average expression for 
  ## ref: https://github.com/satijalab/seurat/issues/5629
  Ave_expr <- Seurat::AverageExpression(object = singlet_harmonized_PBMC, layer = "counts", return.seurat = TRUE,
                                        group.by = "type", feature = sort(union(sig_res_12M$gene, sig_res_6M$gene)))
  

  cell.order <- c("CD4 T_fresh", "CD4 T_6M", "CD4 T_12M", "CD8 T_fresh", "CD8 T_6M", "CD8 T_12M",
                  "B_fresh", "B_6M", "B_12M", "NK_fresh", "NK_6M", "NK_12M",
                  "Mono_fresh", "Mono_6M", "Mono_12M", "DC_fresh", "DC_6M", "DC_12M",
                  "other T_fresh", "other T_6M", "other T_12M", "other_fresh", "other_6M", "other_12M")
  
  mat<- Ave_expr[["RNA"]]@data %>% as.matrix()
  ## scale the rows
  mat<- t(scale(t(mat)))
  
  cluster_anno<- Ave_expr@meta.data$orig.ident
  
  ## make the black color map to 0. the yellow map to highest and the purle map to the lowest
  col_fun = circlize::colorRamp2(c(-1, 0, 3), c("#FF00FF", "black", "#FFFF00"))
  sample_order<-  cell.order
  # plot the heatmap
  meheatmap1<-ComplexHeatmap::Heatmap(mat, name = "Expression",  
                                      column_split = factor(cluster_anno, levels= c("CD4 T", "CD8 T", "B", "NK", "Mono", "DC", "other T", "other")),
                                      cluster_columns = F,
                                      show_column_dend = FALSE,
                                      cluster_column_slices = F,
                                      column_title_gp = gpar(fontsize = 10),
                                      column_gap = unit(0.5, "mm"),
                                      cluster_rows = F,
                                      show_row_dend = FALSE,
                                      col = col_fun,
                                      row_names_gp = gpar(fontsize = 10),
                                      column_title_rot = 90,
                                      top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scales::hue_pal()(9)))),
                                      show_column_names = T,
                                      use_raster = TRUE,
                                      raster_quality = 4,column_order = sample_order,row_names_side = "left")
  
  draw(meheatmap1)
  }
  
  #### DEG analysis on single-cell levels
  # pairwise comparison
  #dir.create("CD8 T/")
  for (i in 1:3) {
    type = "NK"
    A = paste0("PBMC", i, "_12M")
    B = paste0("PBMC", i, "_fresh")
  seu.object <- subset(singlet_harmonized_PBMC, subset = predicted.PBMC_celltype.l1 == type & sample %in% c(A, B))
  # seu.object <- subset(CXCR4_CCR5_CD4T, subset = predicted.PBMC_celltype.l1 == "CD4 T" & sample %in% c(A, B))
  DefaultAssay(seu.object) <- "RNA"
  seu.object$sample <- paste0(seu.object$sample, "_", type)
  Idents(seu.object) <- "sample"
  seu.object <- NormalizeData(seu.object)
  seu.object <- ScaleData(seu.object)
  de.markers <- FindMarkers(seu.object, assay = "RNA", slot = "data", ident.1 = paste0(A, "_", type), ident.2 = paste0(B, "_", type))
  de.markers <- de.markers[order(de.markers$avg_log2FC, decreasing = TRUE),]
  write.csv(de.markers,
            paste0(type, "/", type, "_", A, "_vs_", B,   "_sig_genes.csv"),
            quote = FALSE, 
            row.names = TRUE)
  
  # Single cell heatmap of feature expression
  
  pdf(file = paste0(type, "/", type, "_", A, "_vs_", B,   "_sig_genes.pdf"),
      width = 6,  
      height = 12) 
  
  print(DoHeatmap(seu.object, features = rownames(de.markers), size = 3) + NoLegend() + 
    theme(axis.text.y = element_text(size = 2)))
  
  dev.off()
  }
  
  #### multiple groups comparison
  for (i in 1:3) {
    type = "NK"
    A = paste0("PBMC", i, "_12M")
    B = paste0("PBMC", i, "_6M")
    C = paste0("PBMC", i, "_fresh")
    seu.object <- subset(singlet_harmonized_PBMC, subset = predicted.PBMC_celltype.l1 == type & sample %in% c(A, B, C))
    # seu.object <- subset(CXCR4_CCR5_CD4T, subset = predicted.PBMC_celltype.l1 == "CD4 T" & sample %in% c(A, B, C))
    DefaultAssay(seu.object) <- "RNA"
    seu.object$sample <- paste0(seu.object$sample, "_", type)
    Idents(seu.object) <- "sample"
    seu.object <- NormalizeData(seu.object)
    seu.object <- ScaleData(seu.object)
    de.markers <- FindAllMarkers(seu.object, assay = "RNA", slot = "data",
                                 only.pos = TRUE,
                                 logfc.threshold = 0.25)
    write.csv(de.markers,
              paste0(type, "/", type, "_", A, "_vs_", B, "_vs_", C, "_sig_genes.csv"),
              quote = FALSE, 
              row.names = TRUE)
    
    # Single cell heatmap of feature expression
    pdf(file = paste0(type, "/", type, "_", A, "_vs_", B, "_vs_", C, "_sig_genes.pdf"),
        width = 10,  
        height = 20)
    
    print(DoHeatmap(seu.object, features = rownames(de.markers), size = 3) + NoLegend() + 
      theme(axis.text.y = element_text(size = 5)))
    
    dev.off()
  }
  
  Mono_6M <- read.csv("/hpcdata/dcr_sp/DATA/Tom_PBMC_081723/CD4 T/CD4 T_PBMC1_6M_vs_PBMC1_fresh_sig_genes.csv", row.names = 1)
  Mono_12M <- read.csv("/hpcdata/dcr_sp/DATA/Tom_PBMC_081723/CD4 T/CD4 T_PBMC1_12M_vs_PBMC1_fresh_sig_genes.csv", row.names = 1)
  intersect(rownames(Mono_6M), rownames(Mono_12M))
  setdiff(rownames(Mono_6M), rownames(Mono_12M))
  setdiff(rownames(Mono_12M), rownames(Mono_6M))
  
  # Volcano plot
  # Set thresholds
  log2fc_cutoff <- 0.58
  
  res_table_thres <- res_tbl[!is.na(res_tbl$padj), ] %>% 
    mutate(threshold = case_when(
      padj < padj_cutoff & log2FoldChange >= log2fc_cutoff ~ "A",
      padj < padj_cutoff & log2FoldChange <= -log2fc_cutoff ~ "B",
      TRUE ~ "C"))
  
  ## Generate plot
  ggplot(res_table_thres) +
    geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
    ggtitle(paste0("Volcano plot of PBMC cells 6M to fresh: ",  clustx)) +
    xlab("log2 fold change") +
    # xlim(-4.5, 12) +
    ylab("-log10 adjusted p-value") +
    # scale_color_manual(values = c("red3", "grey70")) +
    scale_color_manual(values = c("red3", "blue1", "grey70")) +
    geom_text_repel(# data=res_table_thres[which(res_table_thres$gene %in% top50_sig_df$gene),],
      data=res_table_thres[res_table_thres$threshold %in% c("A", "B"),],
      aes(x = log2FoldChange, y = -log10(padj),label=gene)) +
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.3), hjust = 0.5),
          axis.title = element_text(size = rel(1.15))) +
    geom_hline(yintercept = -log10(0.05), linetype="dotted", alpha = 0.7) +
    geom_vline(xintercept = -0.58, linetype="dotted", alpha = 0.7) + 
    geom_vline(xintercept = 0.58, linetype="dotted", alpha = 0.7) +
    guides(colour = FALSE) +
    theme_minimal()
  
}


#### run DESeq2 on all cell types (clusters) and all levels of a condition - Likelihood Ratio Test
# ref: https://hbctraining.github.io/scRNA-seq_online/lessons/pseudobulk_DESeq2_scrnaseq.html
require(DEGreport)

# Function to run DESeq2 LRT and get results for any cluster:
## clustx is the name of the cluster (cell type) on which to run the function

## This function assumes the counts matrices and metadata for all clusters have been prepared
## and arranged in matching named lists (as illustrated in tutorial above)
## This function assumes the contrasted groups (e.g. stim A, stim B, control...) are stored in a variable named "group_id"


get_dds_LRTresults <- function(clustx){
  clustx = "Mono"
  print(clustx) # useful for debugging
  
  # Extract counts matrix and metadata for cluster x
  idx <- which(names(counts_ls) == clustx)
  cluster_counts <- counts_ls[[idx]]
  cluster_metadata <- metadata_ls[[idx]]
  
  # Print error message if sample names do not match
  if ( all(colnames(cluster_counts) != rownames(cluster_metadata)) ) {
    print("ERROR: sample names in counts matrix columns and metadata rows do not match!")
  }
  
  # Run DESeq2
  dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                colData = cluster_metadata, 
                                design = ~ condition)
  dds$condition <- factor(dds$condition, levels = c("fresh","6M", "12M"))
  dds$condition <- relevel(dds$condition, ref = "fresh")

  dds_lrt <- DESeq(dds, test = "LRT", reduced = ~1)
  
  # Extract results
  res_LRT <- results(dds_lrt)
  
  # Create a tibble for LRT results
  res_LRT_tb <- res_LRT %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>% 
    as_tibble()
  

  # Subset to return genes with padj < 0.05
  sigLRT_genes <- res_LRT_tb %>% 
    dplyr::filter(padj < 0.05)
  

  # Transform counts for data visualization
  rld <- rlog(dds_lrt, blind = TRUE)

  # Extract the rlog matrix from the object and compute pairwise correlation values
  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)
  
  # Obtain rlog values for those significant genes
  cluster_rlog <- rld_mat[sigLRT_genes$gene, ]
  cluster_meta_sig <- cluster_metadata[which(rownames(cluster_metadata) %in% colnames(cluster_rlog)), ]
  
  # Use the `degPatterns` function from DEGreport package to show gene clusters across sample groups
  cluster_groups <- degPatterns(cluster_rlog, metadata = cluster_meta_sig,
                                time = "condition", col = NULL)
  
  # Extract the Group 1 genes
  cluster <- cluster_groups$df

}

map(cluster_names, get_dds_LRTresults)


#### David_pathway analysis
require(stringr)
file <- list.files("DAVID_pathway/") |> str_subset("B.chartReport.txt")
data <- read.table(paste0("DAVID_pathway/", file), sep = "\t", header = T) |> filter(Benjamini < 0.05) |>
  # dplyr::filter(grepl("stress|MAPK|death|apoptosis|immun|inflama|NF-kappa", Term))
  dplyr::filter(Count >3)
data$Term <- gsub("[WRG].*~|[h|W|S|I].*:", "", data$Term)
data <- data |> dplyr::select(Category, Term, Count, Pvalue, Benjamini, FDR)
data$p.adjust <- -log(data$Benjamini, base=10)
data <- data[order(data$p.adjust, decreasing = T), ]
data <- data[!duplicated(data$Term), ]
p1 <- ggplot(data=data, aes(x= p.adjust, y=reorder(Term, p.adjust), fill = Category)) +
  geom_bar(stat="identity") +
  ylab("pathways") +
  xlab("-log10(p.adjust)") +
  ggtitle(gsub("\\()" , "", str_extract(f[1], "^.*\\)")))

data$Term |> str_subset("stress")
data$Term |> str_subset("stress")
