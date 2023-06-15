# resource
# https://satijalab.org/seurat/articles/spatial_vignette.html

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.17")
BiocManager::install(c("Seurat", "ggplot2", "patchwork", "dplyr","devtools","svglite"))
devtools::install_github('satijalab/seurat-data') # not available with Bioconductor 3.17

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(svglite)

setwd("/home/jk/seurat/brainDemo")
# vignette data
InstallData("stxBrain")
brain <- LoadData("stxBrain", type = "anterior1")

# import .h5 output from spaceranger as input
data.dir <- "./data"
brain <- Load10X_Spatial(data.dir, filename=".h5",
                         assay="Spatial", 
                         slice="Visium_FFPE_Mouse_Brain_image.jpg"
  # This function expects two files in the following directory structure:
  # ./data/
  #   -.h5
  #   -spatial/
  #     -.jpg
)

# preprocessing: distribution of counts
plot1 <- VlnPlot(brain, features="nCount_Spatial", pt.size=0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features="nCount_Spatial")+ theme(legend.position="right")

pdf("./plot1.pdf", width=4, height=4)
print(plot1)
dev.off()

pdf("./plot2.pdf", width=4, height=4)
print(plot2)
dev.off()


# "sctransform" normalization: apparently better than log transform 
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)


# visualization 1: some gene of interest
pdf("plot3.pdf", width=4, height=4)
SpatialFeaturePlot(brain, features=c("Hpca", "Ttr"), alpha=c(1, 1)) 
dev.off()

theme(legend.text = element_text(size = 0),
        legend.title = element_text(size = 20), legend.key.size = unit(1, "cm"))


# dimension reduction, clustering and more viz
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p1 + p2

# separate by identify 
SpatialDimPlot(brain, cells.highlight = 
               CellsByIdentities(object = brain, idents = c(1:14)),
               facet.highlight = TRUE, ncol = 5)

# interactive plots
SpatialDimPlot(brain, interactive = TRUE)
SpatialFeaturePlot(brain, features = c("Ttr","Hpca"), interactive = TRUE)

# spatially variable features: with prior clustering knowledge
de_markers <- FindMarkers(brain, ident.1 = 5, 
                          ident.2 = 6) 
# ident.? refer to clusters identified above with FindCLusters, RunUMAP and DimPlot/SpatialDimPlot
SpatialFeaturePlot(object = brain, 
                   features = rownames(de_markers)[1:3], 
                   alpha = c(0.1, 1), ncol = 3)

# spatially variable features: without prior knowledge
brain <- FindSpatiallyVariableFeatures(
  brain, assay = "SCT", features = VariableFeatures(brain)[1:1000],
  selection.method = "moransi")

top.features <- head(SpatiallyVariableFeatures(
  brain, selection.method = "moransi"), 6)
SpatialFeaturePlot(brain, features = top.features, 
                   ncol = 3, alpha = c(0.1, 1))

