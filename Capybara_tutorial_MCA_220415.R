###############################################################
library( BiocParallel )
library( DESeq2 )
library( gplots )
library( ggplot2 )
library( AnnotationDbi )
library( org.Hs.eg.db )
library( tidyverse )
library( pheatmap )
library( dplyr )
library( ggrepel )

library( pathview )
library( gage )
library( gageData )
library( pathfindR )

library( RColorBrewer )
library( momr )

library("devtools")
library("Capybara")

library(Seurat)
library(tidyverse)
library(Matrix)
library(stringr)
library(RCurl)
library(cowplot)
library(edgeR)
library(scater)
library(scran)
library(monocle)
library(vcd)
library(cowplot)
library(SingleCellExperiment)
library(AnnotationHub)
library(ensembldb)

library(anndata)
library(rhdf5)
library('R.utils')
###############################################################
# File path
bulk.raw.path <- system.file("extdata", "Bulk Reference Raw.Rds", package = "Capybara")
bulk.rpkm.path <- system.file("extdata", "Bulk Reference RPKM.Rds", package = "Capybara")
# Read the matrices
bulk.raw <- readRDS(bulk.raw.path)
bulk.rpkm <- readRDS(bulk.rpkm.path)

# Read in the pancreatic data file that come with the package
fpath <- system.file("extdata", "baron_dataset.zip", package = "Capybara")
extract.dir <- "."
# Extract the dataset
unzip(fpath, overwrite = FALSE, exdir = ".")
# Identify the full path
full.fpath.meta <- paste0(extract.dir, "/", "baron_et_al_pancreatic_meta.csv")
full.fpath.raw <- paste0(extract.dir, "/", "baron_et_al_pancreatic.csv")
# Load the count matrix and the meta data
baron.expr <- read.csv(full.fpath.raw, header = T, row.names = 1, stringsAsFactors = F)
baron.meta <- read.csv(full.fpath.meta, header = T, row.names = 1, stringsAsFactors = F)

single.round.QP.analysis(bulk.raw, baron.expr, scale.bulk.sc = "scale", unix.par = TRUE, 
                         force.eq = 1, n.cores = 4, save.to.path = "./", 
                         save.to.filename = "baron_bulk_classification_qp")


# Read the meta data
mca.meta.fpath <- system.file("extdata", "MCA_CellAssignments.csv", package = "Capybara")
mca <- read.csv(mca.meta.fpath, row.names = 1, header = T, stringsAsFactors = F)
# Clean up the meta data
mca.meta <- data.frame(row.names = mca$Cell.name, 
                       tissue = mca$Tissue,
                       cell.type = mca$Annotation,
                       stringsAsFactors = F)

# List all possible files and tissues in the Mouse Cell Atlas
file.ls <- list.files("~/Stroke/Capybara_Cell_Stem_cell/MCA_data/rmbatch_dge/", full.names = T)
base.nms <- basename(file.ls)
# Identify the tissues
unq.tissue <- unique(base.nms)

# Set a path to save all QP files for all tissues
general.path.to.save <- "~/Stroke/Capybara_Cell_Stem_cell/MCA_All_Tissue_QP_rmbatch/"
for (k in 1:length(unq.tissue)) {
  curr.tissue <- unq.tissue[k]
  curr.tissue = unlist(strsplit(unq.tissue[k], split='_', fixed=TRUE))[1]
  curr.filename <- paste0("0", k, "_", curr.tissue, "_Bulk_ARCHS4")
  
  file.base.name <- base.nms[which(startsWith(base.nms, curr.tissue))][1]
  file.full <- file.ls[which(startsWith(base.nms, curr.tissue))][1]
  
  print(curr.tissue)
  
  sc.data <- read.table(gunzip(file.full,remove=FALSE))#, header = T, row.names = 1, stringsAsFactors = F))#, header = T, row.names = 0, stringsAsFactors = F)
  

 
if (all(is.na(sc.data))) {
    print("There is no data in this counting matrix!")
  } else {
    single.round.QP.analysis(bulk.raw, sc.data, scale.bulk.sc = "scale", unix.par = FALSE, 
                             force.eq = 1, n.cores = 1, save.to.path = general.path.to.save, 
                             save.to.filename = curr.filename)
  }
}

##select 90 cells
# Read the QP files from the directory
qp.files.to.read.clean <- list.files("~/Stroke/Capybara_Cell_Stem_cell/MCA_All_Tissue_QP_rmbatch", full.names = T)

full.qp.mtx.known.annotation <- data.frame()
full.qp.mtx.unknown.annotation <- data.frame()
for (i in 1:length(qp.files.to.read.clean)) {
  curr.file <- qp.files.to.read.clean[i]
  curr.qp.rslt <- read.csv(curr.file, header = T, row.names = 1, stringsAsFactors = F)
  
  cells.to.keep <- intersect(rownames(mca.meta), rownames(curr.qp.rslt))
  cells.unlabel <- setdiff(rownames(curr.qp.rslt), cells.to.keep)
  
  curr.sub.mtx.to.keep <- curr.qp.rslt[cells.to.keep, ]
  curr.sub.mtx.unlabel <- curr.qp.rslt[cells.unlabel, ]
  
  if (nrow(full.qp.mtx.known.annotation) <= 0) {
    full.qp.mtx.known.annotation <- curr.sub.mtx.to.keep
    full.qp.mtx.unknown.annotation <- curr.sub.mtx.unlabel
  } else {
    full.qp.mtx.known.annotation <- rbind(full.qp.mtx.known.annotation, curr.sub.mtx.to.keep)
    full.qp.mtx.unknown.annotation <- rbind(full.qp.mtx.unknown.annotation, curr.sub.mtx.unlabel)
  }
}

full.qp.mtx.known.annotation.qp.score.only <- full.qp.mtx.known.annotation[,c(1:(ncol(full.qp.mtx.known.annotation) - 2))]
# Create a map between MCA and ARCHS4
map.df <- data.frame(tm.tissue = c("Embryonic-Mesenchyme", "Embryonic-Stem-Cell", "Trophoblast-Stem-Cell", "Fetal_Brain",
                                    "Neonatal-Calvaria","Fetal_Intestine", "Fetal-Liver", "Fetal_Lung", "Fetal_Stomache",
                                    "Neonatal-Heart", "Neonatal-Muscle",
                                    "Neonatal-Rib", "Neonatal-Skin",  "NeonatalPancreas"),
                     corresponding = c("frxn_embryo", "frxn_embryo", "frxn_embryo", "frxn_brain","frxn_brain",
                                       "frxn_small.intestine", "frxn_liver", 
                                       "frxn_lung", "frxn_stomach",  "frxn_heart", "frxn_muscle", "frxn_muscle", 
                                       "frxn_skin", "frxn_pancreas"),
                     stringsAsFactors = F)

# Identify top 90 cells for each tissue
tm.tissue <- unique(map.df$tm.tissue)
cell.selector <- c()
n.sample <- 90
for (i in 1:length(tm.tissue)) {
  curr.tissue <- tm.tissue[i]
  cell.names <- rownames(mca.meta)[which(mca.meta$tissue == curr.tissue)]
  curr.qp.subset <- full.qp.mtx.known.annotation.qp.score.only[cell.names, ]
  curr.map <- map.df$corresponding[which(map.df$tm.tissue == curr.tissue)]
  if (length(curr.map) <= 1){
    curr.qp.subset.sub <- data.frame(score = curr.qp.subset[,curr.map], cell.name = cell.names, stringsAsFactors = F)
  } else {
    curr.qp.subset.sub <- data.frame(score = rowSums(curr.qp.subset[,curr.map]), cell.name = cell.names, stringsAsFactors = F)
  }
  curr.qp.subset.sub.sort <- curr.qp.subset.sub[order(-curr.qp.subset.sub$score), ]
  cells.to.incl <- curr.qp.subset.sub.sort$cell.name[1:n.sample]
  
  cell.selector <- c(cell.selector, cells.to.incl)
}
saveRDS(full.qp.mtx.known.annotation.qp.score.only[cell.selector, ], "./MCA_embryonic_background_rmbatch.RDS")
###################################################################
#load QP background matrix
background.qp.fpath <- system.file("extdata", "MCA Adult Background.Rds", package = "Capybara")
background.mtx <- readRDS(background.qp.fpath)

#Load the QP scores of the sample
## Load QP results
qp.rslt <- read.csv("~/Stroke/Capybara_Cell_Stem_cell/baron_bulk_classification_qp_scale.csv", row.names = 1, header = T, stringsAsFactors = F)

## Reshape the data
qp.rslt.sub <- qp.rslt[,c(1:(ncol(qp.rslt) - 2))]

#correlation calculation
mtx.test <- t(qp.rslt.sub[, colnames(background.mtx)])
ref.test <- t(background.mtx)

# Pearson's Correlation Calculation
corr.mtx <- WGCNA::cor(ref.test, mtx.test)

#Binarization based on correlation

# Setup a correlation cutoff to the 90th quantile of the correlation matrix
correlation.cutoff <- quantile(corr.mtx, 0.90)

# Binarization based on the correlation
new.corr.bin <- corr.mtx
new.corr.bin[which(new.corr.bin >= correlation.cutoff)] <- 1
new.corr.bin[which(new.corr.bin < correlation.cutoff)] <- 0
new.corr.bin <- as.data.frame(new.corr.bin)

# count tissue
# Count
count.in.cat <- c()
unique.cat <- unique(unlist(lapply(strsplit(rownames(new.corr.bin), "_"), function(x) x[1])))
for (uc in unique.cat) {
  curr.subset <- new.corr.bin[which(startsWith(rownames(new.corr.bin), uc)), c(1:1886)]
  count.in.cat[uc] <- sum(colSums(curr.subset) >= nrow(curr.subset) * 0.7)
}

count.in.cat <- as.data.frame(count.in.cat)
count.in.cat$perc <- round(count.in.cat$count.in.cat *100/sum(count.in.cat$count.in.cat), digits = 3)

# Check frequency
final.cell.types.adult <- rownames(count.in.cat)[which(count.in.cat$count.in.cat > 100)]
#####################################################
##Step 2: Generation of a High-Resolution Custom Reference, and Continuous Identity Measurement
#Get the counts of the cell types involved in the tissues selected

pancreatic.all.meta <- mca.meta[which(mca.meta$tissue %in% final.cell.types.adult), ]

stomach_ref = read.table("~/Stroke/Capybara_Cell_Stem_cell/MCA_data/rmbatch_dge/Stomach_rm.batch_dge.txt")
pancreas_ref = read.table("~/Stroke/Capybara_Cell_Stem_cell/MCA_data/rmbatch_dge/Pancreas_rm.batch_dge.txt")
stomach_ref$gene_id = rownames(stomach_ref)
pancreas_ref$gene_id = rownames(pancreas_ref)
stomach_ref = stomach_ref[,c(ncol(stomach_ref),1:(ncol(stomach_ref)-1))]
pancreas_ref = pancreas_ref[,c(ncol(pancreas_ref),1:(ncol(pancreas_ref)-1))]
tissues = merge(stomach_ref,pancreas_ref, by = "gene_id", all=T)
tissues[is.na(tissues)] <- 0
rownames(tissues) = tissues$gene_id
tissues = tissues[,-1]

mca.counts.all.involved = tissues

## Meta data filtering
pancreatic.all.meta$cell.type <- gsub("Ductal.cell", "ductal", pancreatic.all.meta$cell.type)
pancreatic.all.meta$cell.type <- gsub("β.cell", "beta", pancreatic.all.meta$cell.type)
pancreatic.all.meta$cell.type <- gsub("Endothelial.cell", "endothelial", pancreatic.all.meta$cell.type)
pancreatic.all.meta$cell.type.1 <- gsub("\\([^)]*\\)", "", pancreatic.all.meta$cell.type)
pancreatic.all.meta$cell.type.alone <- unlist(lapply(strsplit(pancreatic.all.meta$cell.type.1, "_"), function(x) x[1]))

## Filter out cell types with less than 30 cells
cell.type.alone.freq <- as.data.frame(table(pancreatic.all.meta$cell.type.alone))
cell.type.over.30 <- cell.type.alone.freq$Var1[which(cell.type.alone.freq$Freq >= 30)]
pancreatic.sub.meta <- pancreatic.all.meta[which(pancreatic.all.meta$cell.type.alone %in% as.character(cell.type.over.30)),]
coldata.df <- pancreatic.sub.meta

# Construction of a high-resolution reference
ref.list <- construct.high.res.reference(mca.counts.all.involved, coldata.df = coldata.df, criteria = "cell.type.alone")
# Get expression matrix and meta data of cells used to build the reference, as well as the constructed pseudo-bulk reference
ref.df <- ref.list[[3]]
ref.meta <- ref.list[[2]]
ref.sc <- ref.list[[1]]

# Measure cell identity in the reference dataset as a background 
single.round.QP.analysis(ref.df, ref.list[[1]], n.cores = 4, save.to.path = "./", save.to.filename = "01_MCA_Based_scClassifier_reference_mix90_normalize_select", unix.par = TRUE)

# Measure cell identity in the query dataset 
single.round.QP.analysis(ref.df, baron.expr, n.cores = 4, save.to.path = "./", save.to.filename = "02_MCA_Based_scClassifier_reference_mix90_test_normalize_select", unix.par = TRUE)

##############################################
##Discrete Cell Type Classification and Multiple Identity Scoring

#Empirical p-value calculation
# Read in background and testing identity scores
background.mtx <- read.csv("./01_MCA_Based_scClassifier_reference_mix90_normalize_select_scale.csv", header = T, row.names = 1, stringsAsFactors = F)
mtx.test <- read.csv("./02_MCA_Based_scClassifier_reference_mix90_test_normalize_select_scale.csv", header = T, row.names = 1, stringsAsFactors = F)

col.sub <- ncol(background.mtx) - 2

# Conduct reference randomization to get empirical p-value matrix
ref.perc.list <- percentage.calc(background.mtx[,c(1:col.sub)], background.mtx[,c(1:col.sub)])

# Conduct test randomization to get empirical p-value matrix
perc.list <- percentage.calc(as.matrix(mtx.test[,c(1:col.sub)]), as.matrix(background.mtx[,c(1:col.sub)]))

# Binarization of inference results
bin.count <- binarization.mann.whitney(mtx = mtx.test[,c(1:col.sub)], ref.perc.ls = ref.perc.list, ref.meta = ref.list[[2]], perc.ls = perc.list)

#Classification
classification <- binary.to.classification(bin.count[,c(1:col.sub)])
rownames(classification) <- classification$barcode

#Check with classification result
rownames(baron.meta) = gsub("-", ".", rownames(baron.meta))
classification$actual <- baron.meta[rownames(classification), "cell.type"]

table.freq <- table(classification$actual, classification$call)
#table.freq["Var1", 16] <- "beta"
table.freq.perc <- apply(table.freq, 1, function(x) round(x * 100/sum(x), digits = 3))



table.freq.sub <- as.data.frame(table.freq.perc[c("beta", "ductal", "endothelial",
                                                  "Macrophage", "T.cell", "Dendritic.cell", 
                                                  "Multi_ID", "Endocrine.cell"),
                                                c("B_cell", "beta", "ductal", "endothelial",
                                                  "macrophage", "T_cell", "immune_other", "activated_stellate", 
                                                  "alpha", "delta", "gamma")])
table.freq.sub$Capybara.Call <- rownames(table.freq.sub)
table.freq.melt <- melt(table.freq.sub)

table.freq.melt$Capybara.Call <- factor(table.freq.melt$Capybara.Call,
                                        levels = c("B_cell", "beta", "ductal", "endothelial",
                                                   "Macrophage", "T.cell", "Dendritic.cell", 
                                                   "Multi_ID", "Endocrine.cell"),
                                        ordered = T)
table.freq.melt$variable <- factor(table.freq.melt$variable,
                                   levels = c("B_cell", "beta", "ductal", "endothelial",
                                              "macrophage", "T_cell", "immune_other", "activated_stellate", 
                                              "alpha", "delta", "gamma"),
                                   ordered = T)

ggplot(table.freq.melt, aes(x = Capybara.Call, y = variable, size=ifelse(value==0, NA,  value))) +
  geom_point(aes(colour = variable)) +
  scale_size_area(name = "Percentage", max_size=12) +
  scale_color_viridis_d(option = "A", begin = 0.15, end = 0.85) +
  ggtitle("Mouse Pancreatic Dataset (Baron et al., 2016)") +
  guides(fill = FALSE, color = FALSE) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 12, face = "bold.italic", angle = 90),
        axis.text.y = element_text(size = 12, face = "bold.italic"),
        axis.ticks = element_blank(), 
        axis.title = element_blank(), 
        title = element_text(face = "bold.italic", size = 14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black", size = 1))


#################################################
###Analysis of Cells with Multiple Identities

#downloaddata
#download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE133452&format=file&file=GSE133452%5Fm1%5F1%5F2%5F3%5F7%5F14P%5Fpaper%2Ecsv%2Egz", "./cardiomyocyte_reprogramming_m1_14p.csv.gz")
#Go to terminal and unzip the file with command:
#gzip -d -k cardiomyocyte_reprogramming_m1_14p.csv.gz

##Preprocessing of the data with Seurat

# Read in the file path for all features and genes
feature.df <- read.csv("~/Stroke/Capybara_Cell_Stem_cell/features.csv", row.names = 1, header = F, stringsAsFactors = F)

# Load the data
stone.et.al <- read.csv("~/Stroke/Capybara_Cell_Stem_cell/cardiomyocyte_reprogramming_m1_14p.csv", row.names = 1, header = T)#, stringsAsFactors = F)

# Map the gene names fr
gene.name.subset <- feature.df[intersect(rownames(stone.et.al), rownames(feature.df)), ]
stone.et.al.subset <- stone.et.al[which(rownames(stone.et.al) %in% rownames(feature.df)), ]
stone.et.al.subset$gene.name <- gene.name.subset[rownames(stone.et.al.subset), "V2"]
stone.et.al.subset <- stone.et.al.subset[-which(duplicated(stone.et.al.subset$gene.name)), ]
rnm <- stone.et.al.subset$gene.name
stone.et.al.final <- stone.et.al.subset[, -c(1,ncol(stone.et.al.subset))]
rownames(stone.et.al.final) <- rnm

# Create Seurat object
sc.data.stone <- CreateSeuratObject(counts = stone.et.al.final, project = "cardiac.reprog", min.cells = 3, min.features = 200)

# Calculate mitochondria content
sc.data.stone[["percent.mt"]] <- PercentageFeatureSet(sc.data.stone, pattern = "mt-")

# Visualize QC metrics as a violin plot and scatter plot
VlnPlot(sc.data.stone, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(sc.data.stone, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sc.data.stone, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# Filter the dataset based on number of features
sc.data.stone <- subset(sc.data.stone, subset = nFeature_RNA > 200 & nFeature_RNA < 5500)

# Log normalize the data
sc.data.stone <- NormalizeData(sc.data.stone, normalization.method = "LogNormalize", scale.factor = 10000)

# Variable gene identification
sc.data.stone <- FindVariableFeatures(sc.data.stone, selection.method = "vst", nfeatures = 2000)

# Scale the data
all.genes <- rownames(sc.data.stone)
sc.data.stone <- ScaleData(sc.data.stone, features = all.genes)

# PCA
sc.data.stone <- RunPCA(sc.data.stone, features = VariableFeatures(object = sc.data.stone))

# JackStraw procedure and Elbow plot to select number of PCs
sc.data.stone <- JackStraw(sc.data.stone, num.replicate = 100)
sc.data.stone <- ScoreJackStraw(sc.data.stone, dims = 1:20)

JackStrawPlot(sc.data.stone, dims = 1:20)
ElbowPlot(sc.data.stone)

# Identify neighbors and clusters
sc.data.stone <- FindNeighbors(sc.data.stone, dims = 1:18)
sc.data.stone <- FindClusters(sc.data.stone, resolution = 0.8)

# UMAP embedding
sc.data.stone <- RunUMAP(sc.data.stone, dims = 1:18)
DimPlot(sc.data.stone)#, reduction = 'umap')

##Classification
#Tissue classification

# File path
bulk.raw.path <- system.file("extdata", "Bulk Reference Raw.Rds", package = "Capybara")
bulk.rpkm.path <- system.file("extdata", "Bulk Reference RPKM.Rds", package = "Capybara")
# Read the matrices
bulk.raw <- readRDS(bulk.raw.path)
bulk.rpkm <- readRDS(bulk.rpkm.path)
#QP
single.round.QP.analysis(bulk.raw, stone.et.al.final, scale.bulk.sc = "scale", unix.par = TRUE, 
                         force.eq = 1, n.cores = 4, save.to.path = "./", 
                         save.to.filename = "stone_bulk_classification_qp")
#Correction to analysist
## Load QP results
qp.rslt <- read.csv("./stone_bulk_classification_qp_scale.csv", row.names = 1, header = T, stringsAsFactors = F)

## Reshape the data
qp.rslt.sub <- qp.rslt[,c(1:(ncol(qp.rslt) - 2))]

## Background matrix
background.qp.fpath <- system.file("extdata", "MCA Embryonic Background.Rds", package = "Capybara")
background.mca <- readRDS(background.qp.fpath)
background.mtx <- background.mca[[2]]

## Correlation Analysis
mtx.test <- t(qp.rslt.sub[, colnames(background.mtx)])
ref.test <- t(background.mtx)

## Pearson's Correlation Calculation
corr.mtx <- WGCNA::cor(ref.test, mtx.test)

## Setup a correlation cutoff to the 90th quantile of the correlation matrix
correlation.cutoff <- quantile(corr.mtx, 0.90)

## Binarization based on the correlation
new.corr.bin <- corr.mtx
new.corr.bin[which(new.corr.bin >= correlation.cutoff)] <- 1
new.corr.bin[which(new.corr.bin < correlation.cutoff)] <- 0
new.corr.bin <- as.data.frame(new.corr.bin)

#Mapping to Tissues in MCA
# Count
count.in.cat <- c()
unique.cat <- unique(unlist(lapply(strsplit(rownames(new.corr.bin), "_"), function(x) x[1])))
for (uc in unique.cat) {
  curr.subset <- new.corr.bin[which(startsWith(rownames(new.corr.bin), uc)), c(1:30728)]
  count.in.cat[uc] <- sum(colSums(curr.subset) >= nrow(curr.subset) * 0.80)
}

count.in.cat <- as.data.frame(count.in.cat)
count.in.cat$perc <- round(count.in.cat$count.in.cat *100/sum(count.in.cat$count.in.cat), digits = 3)

final.cell.types.fetal <- rownames(count.in.cat)[which(count.in.cat$count.in.cat > 100)]


##High-resolution reference
#get background
# Background cells
mca <- read.csv("~/Stroke/Capybara_Cell_Stem_cell/MCA_data/MCA_CellAssignments.csv",
                row.names = 1, header = T, stringsAsFactors = F)
mca.meta <- data.frame(row.names = mca$Cell.name, 
                       tissue = mca$Tissue,
                       cell.bc.tissue = unlist(lapply(strsplit(mca$Cell.name, "_"), function(x) x[1])),
                       cell.type = mca$Annotation,
                       stringsAsFactors = F)

cardiac.rp.all.meta <- mca.meta[which(mca.meta$cell.bc.tissue %in% final.cell.types.fetal), ]

mca.counts.all.involved <- NULL

fetalstomach_ref = read.table("~/Stroke/Capybara_Cell_Stem_cell/MCA_data/rmbatch_dge/FetalStomach_rm.batch_dge.txt")
neonatalheart_ref = read.table("~/Stroke/Capybara_Cell_Stem_cell/MCA_data/rmbatch_dge/NeonatalHeart_rm.batch_dge.txt")
fetalstomach_ref$gene_id = rownames(fetalstomach_ref)
neonatalheart_ref$gene_id = rownames(neonatalheart_ref)
fetalstomach_ref = fetalstomach_ref[,c(ncol(fetalstomach_ref),1:(ncol(fetalstomach_ref)-1))]
neonatalheart_ref = neonatalheart_ref[,c(ncol(neonatalheart_ref),1:(ncol(neonatalheart_ref)-1))]
tissues = merge(fetalstomach_ref,neonatalheart_ref, by = "gene_id", all=T)
tissues[is.na(tissues)] <- 0
rownames(tissues) = tissues$gene_id
tissues = tissues[,-1]

mca.counts.all.involved = tissues

## meta data cleaning
cardiac.rp.all.meta$cell.type.1 <- gsub("\\([^)]*\\)", "", cardiac.rp.all.meta$cell.type)
cardiac.rp.all.meta$cell.type.alone <- unlist(lapply(strsplit(cardiac.rp.all.meta$cell.type.1, "_"), function(x) x[1]))

cardiac.rp.all.meta$cell.type.1 <- tolower(cardiac.rp.all.meta$cell.type.1)
coldata.df <- cardiac.rp.all.meta

# Construction of a high-resolution reference
ref.list <- construct.high.res.reference(mca.counts.all.involved, coldata.df = coldata.df, criteria = "cell.type.1")
# Get expression matrix and meta data of cells used to build the reference, as well as the constructed pseudo-bulk reference
ref.df = ref.list[[3]]
ref.meta = ref.list[[2]]
ref.sc = ref.list[[1]]

#QP
single.round.QP.analysis(ref.df, ref.list[[1]], n.cores = 4, save.to.path = "./", save.to.filename = "stone_et_al_reference_MCA", unix.par = TRUE)
single.round.QP.analysis(ref.df, stone.et.al.final, n.cores = 4, save.to.path = "./", save.to.filename = "stone_et_al_test_MCA", unix.par = TRUE)

##Discrete Cell type
# Read in background and testing identity scores
background.mtx <- read.csv("./stone_et_al_reference_MCA_scale.csv", header = T, row.names = 1, stringsAsFactors = F)
mtx.test <- read.csv("./stone_et_al_test_MCA_scale.csv", header = T, row.names = 1, stringsAsFactors = F)

col.sub <- ncol(background.mtx) - 2

# Conduct reference randomization to get empirical p-value matrix
ref.perc.list <- percentage.calc(background.mtx[,c(1:col.sub)], background.mtx[,c(1:col.sub)])

# Conduct test randomization to get empirical p-value matrix
perc.list <- percentage.calc(as.matrix(mtx.test[,c(1:col.sub)]), as.matrix(background.mtx[,c(1:col.sub)]))

# Binarization of inference results
bin.count <- binarization.mann.whitney(mtx = mtx.test[,c(1:col.sub)], ref.perc.ls = ref.perc.list, ref.meta = ref.list[[2]], perc.ls = perc.list)
# Classificationn
classification <- binary.to.classification(bin.count[,c(1:col.sub)])
rownames(classification) <- classification$barcode

#filter multiple ID cells
multi.classification.list <- multi.id.curate.qp(binary.counts = bin.count, classification = classification, qp.matrix = mtx.test)
# Reassign variables
actual.multi <- multi.classification.list[[1]]
new.classification <- multi.classification.list[[2]]

##calculate transition scores
score.df <- transition.score(actual.multi)

score_to_merge = score.df
score_to_merge$call = rownames(score_to_merge)
classification_with_score = merge(score_to_merge, new.classification, by = "call", all = F)
classification_with_score$tp = classification_with_score$barcode
classification_with_score$tp = substr(classification_with_score$tp,nchar(classification_with_score$tp),nchar(classification_with_score$tp))
classification_with_score_filtered = classification_with_score[which(classification_with_score$barcode %in% colnames(stone.et.al.final)), ]
classification_with_score_filtered$tp = as.numeric(classification_with_score_filtered$tp)

drawumap <- CreateSeuratObject(counts = stone.et.al.final, project = "cardiac.reprog", min.cells = 3, min.features = 200)
drawumap@meta.data = classification_with_score_filtered

# Filter the dataset based on number of features
drawumap <- subset(drawumap, subset = nFeature_RNA > 200 & nFeature_RNA < 5500)

# Log normalize the data
drawumap <- NormalizeData(drawumap, normalization.method = "LogNormalize", scale.factor = 10000)

# Variable gene identification
drawumap <- FindVariableFeatures(drawumap, selection.method = "vst", nfeatures = 2000)

# Scale the data
all.genes <- rownames(drawumap)
drawumap <- ScaleData(drawumap, features = all.genes)

# PCA
drawumap <- RunPCA(drawumap, features = VariableFeatures(object = drawumap))

# JackStraw procedure and Elbow plot to select number of PCs
drawumap <- JackStraw(drawumap, num.replicate = 100)
drawumap <- ScoreJackStraw(drawumap, dims = 1:20)

# Identify neighbors and clusters
drawumap <- FindNeighbors(drawumap, dims = 1:18)
drawumap <- FindClusters(drawumap, resolution = 0.8)
drawumap <- RunUMAP(drawumap, dims = 1:18)
DimPlot(drawumap, reduction = 'umap')
