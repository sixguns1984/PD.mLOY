
# scRNA-seq analysis

library(Seurat)
library(cowplot)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(dplyr)
library(parallel)

dim.usage <-  30
pc.usage <- 50
seed.usage <-  123
k.usage <-  20
maxdim.usage <-  "2L"

### input RDS files

files <- list.dirs( recursive = F)

ReadData <- function(files,i){
  seurat_data <- Read10X(data.dir = files[i])
  sample <- gsub("./","",files[i])
  colnames(seurat_data) <- paste0(sample,"_",colnames(seurat_data))
  seurat_data <- CreateSeuratObject(counts = seurat_data,
                                    min.features = 100,
                                    project = files[i])
  MTgene <- length(grep ("^MT-", rownames(seurat_data[["RNA"]]),value = T))
  seurat_data[["percent.mt"]] <- PercentageFeatureSet(seurat_data, pattern = "^MT-")
  seurat_data
}

objectlist <- mclapply(1:length(files), function(i) try(ReadData(files = files,i = i), TRUE),
                       mc.cores = 8)
names(objectlist) <- gsub("./","",files)

### Seurat Integrate
objectlist <- lapply(X = objectlist, FUN = function(x) {
  x <- subset(x, subset = nCount_RNA >= 800 & nFeature_RNA >= 350 & percent.mt < 10)
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

object.anchors <- FindIntegrationAnchors(object.list = objectlist, dims = 1:dim.usage)
object.combined <- IntegrateData(anchorset = object.anchors, dims = 1:dim.usage)
DefaultAssay(object.combined) <- "integrated"

object.combined <- CellCycleScoring(
  object.combined,
  g2m.features = cc.genes$g2m.genes,
  s.features = cc.genes$s.genes
)


all.genes <- rownames(object.combined)

object.combined <- ScaleData(object.combined, features = all.genes,
                             vars.to.regress = c("nFeature_RNA", "percent.mt","S.Score", "G2M.Score"))

object.combined <- RunPCA(object.combined, features = VariableFeatures(object = object.combined),
                          npcs = pc.usage, seed.use = seed.usage,verbose = FALSE)


## 
ElbowPlot(object.combined,ndims = 40, reduction = "pca") 
# pc = 22

object.combined <- RunUMAP(object.combined, reduction = "pca", dims = 1:22, seed.use = seed.usage, max.dim = maxdim.usage)
object.combined <- RunTSNE(object.combined, reduction = "pca", dims = 1:22, seed.use = seed.usage)

DimPlot(object.combined, reduction = "pca")

object.combined <- FindNeighbors(object.combined, reduction = "pca", dims = 1:22, k.param = k.usage, annoy.metric = "euclidean")

obj <- FindClusters(object.combined, resolution = seq(0.4,1.2,by=0.2))

library("clustree")
clustree(obj)
##  resolution=0.7

object.combined <- FindClusters(object.combined, resolution = 0.7, algorithm = 1, random.seed = seed.usage)
DefaultAssay(object.combined) <- "RNA"

seurat_markers <- FindAllMarkers(object.combined,logfc.threshold = 0.25,min.pct = 0.5,only.pos = T)

object.combined_for_annotation <- GetAssayData(object.combined, slot="data")

library("scHCL")
object.combined_cell_type <- scHCL(scdata = object.combined_for_annotation, numbers_plot = 3)

out=as.data.frame(unlist(object.combined_cell_type$scHCL))

object.combined@meta.data$scHCL = out[match(rownames(object.combined@meta.data),rownames(out)),1]
saveRDS(object.combined,file = "object.combined.RDS")

object.combined <- readRDS("object.combined.RDS")

DefaultAssay(object.combined) <- "RNA"


Idents(object.combined) <- object.combined$split
levels(object.combined)
new.ids <- c('mLOY','nonmLOY','mLOY','mLOY',
             'nonmLOY','mLOY','nonmLOY','nonmLOY')
names(new.ids) <- levels(object.combined)
object.combined <- RenameIdents(object.combined, new.ids)
object.combined$condition <- Idents(object.combined)

Idents(object.combined) <- object.combined$split
levels(object.combined)

new.ids <- c(74,52,59,56,
             53,74,70,76)
names(new.ids) <- levels(object.combined)
object.combined <- RenameIdents(object.combined, new.ids)
object.combined$Age <- Idents(object.combined)
# 1

object.combined <- ScaleData(object.combined, features = VariableFeatures(object.combined),
                                  vars.to.regress = c("nFeature_RNA", "percent.mt","S.Score", "G2M.Score"))
saveRDS(object.combined,"object.combined_V2.RDS")


## 计算LOY的另一种方法
object.combined <- readRDS("object.combined_V2.RDS")

Y_scLOY_genes <- read.csv("chrY_gene2.csv")

Y_scLOY_genes$gene_name -> Y_genes
Y_genes[Y_genes %in% rownames(object.combined)] -> Y_genes

object.combined@meta.data -> df
object.combined@assays$RNA@counts -> counts
counts[Y_genes,,drop=F] -> counts_chrY
object.combined@assays$RNA@data -> norm
norm[Y_genes,,drop=F] -> norm_chrY

as.data.frame(colMeans(as.matrix(norm_chrY))) -> chrY_norm; names(chrY_norm) <- "chrY_norm"
as.data.frame(colSums(as.matrix(counts_chrY))) -> chrY_UMI; names(chrY_UMI) <- "chrY_UMI"
as.data.frame(colSums(as.matrix(counts_chrY > 0))) ->  chrY_genes; names(chrY_genes) <- "chrY_genes"
as.data.frame(colSums(as.matrix(counts_chrY > 0))) ->  PAR_genes; names(PAR_genes) <- "PAR_genes"

df$chrY_genes <- chrY_genes$chrY_genes
df$chrY_UMI <- chrY_UMI$chrY_UMI
df$chrY_norm <- chrY_norm$chrY_norm

# tibble::rownames_to_column(df, "CB") -> df
dplyr::mutate(df, LOY = ifelse(test = chrY_UMI == 0 , yes = "LOY", no = "NORMAL"))  -> df
dplyr::mutate(df, LOY_levels = ifelse(test = chrY_UMI > 2 & chrY_genes > 2 , yes = 4,
                                      no = ifelse(test = chrY_UMI > 1 & chrY_genes > 1, yes = 3,
                                                  no = ifelse(test = chrY_UMI > 1 | chrY_genes > 1, yes = 2,
                                                              no = ifelse(test = chrY_UMI >0, yes = 1, no = 0)))))  -> df
ifelse(df$LOY_levels == 0, yes = 0, no = ifelse(df$LOY_levels %in% c(1,2), yes = 1, no = 2) ) -> df$LOY_levels_condensed

dplyr::mutate(df, normY = log1p((chrY_UMI / nCount_RNA) * 10000)) -> df
object.combined@meta.data <- df

##  gene marker注释cell
library(Seurat)
library(data.table)
library(RColorBrewer)
library(dplyr)

object.combined_left <- readRDS("object.combined_V2.RDS")

features = c("MRC1","ITGAX","FCGR1A", "CD14",
             "CD163")

VlnPlot(object.combined_left, features = features,pt.size = 0)

features = c("LYZ","C1QB","CD3D", "CD3G",
             "MS4A1","CD79A", "KLRF1","KLRD1")


VlnPlot(object.combined_left, features = features,pt.size = 0)


object.combined_left$anno <- as.integer(as.character(object.combined_left$seurat_clusters))

object.combined_left$anno[object.combined_left$anno == 0] = "NK_cell"
object.combined_left$anno[object.combined_left$anno == 1] = "T_cell"
object.combined_left$anno[object.combined_left$anno == 2] = "Myeloid_cell"
object.combined_left$anno[object.combined_left$anno == 3] = "T_cell"
object.combined_left$anno[object.combined_left$anno == 4] = "NK_cell"
object.combined_left$anno[object.combined_left$anno == 5] = "T_cell"
object.combined_left$anno[object.combined_left$anno == 6] = "T_cell"
object.combined_left$anno[object.combined_left$anno == 7] = "Myeloid_cell"
object.combined_left$anno[object.combined_left$anno == 8] = "B_cell"
object.combined_left$anno[object.combined_left$anno == 9] = "T_cell"
object.combined_left$anno[object.combined_left$anno == 10] = "T_cell"
object.combined_left$anno[object.combined_left$anno == 11] = "B_cell"
object.combined_left$anno[object.combined_left$anno == 12] = "T_cell"
object.combined_left$anno[object.combined_left$anno == 13] = "T_cell"
object.combined_left$anno[object.combined_left$anno == 14] = "T_cell"
object.combined_left$anno[object.combined_left$anno == 16] = "Myeloid_cell"

object.combined_left$anno[object.combined_left$anno == 15] = "Myeloid_cell"
object.combined_left$anno[object.combined_left$anno == 17] = "T_cell"
object.combined_left$anno[object.combined_left$anno == 18] = "Red_blood_cell"
object.combined_left$anno[object.combined_left$anno == 19] = "Myeloid_cell"
object.combined_left$anno[object.combined_left$anno == 20] = "Megakaryocyte"

saveRDS(object.combined_left,"object.combined_V2.RDS")



#####   大类分群   #####

## T cell 



############ T cell  ########

library(Seurat)
library(data.table)
library(RColorBrewer)
library(dplyr)
library(cowplot)
library(tidyr)
library(ggplot2)
# library(scMCA)
library(patchwork)
library(parallel)
library(harmony)

object.t_cell <- readRDS("object.combined_V2.RDS")
object.t_cell <- subset(object.t_cell, anno == "T_cell" )

count <- GetAssayData(object.t_cell, assay = "RNA", slot = "count") # normalized data matrix

meta = object.t_cell@meta.data[,c("nFeature_RNA", "percent.mt","S.Score",
                                  "G2M.Score","split","condition","Age","anno")] # a dataframe with rownames containing cell mata data


object.t_cell <- CreateSeuratObject(
  counts = count)
object.t_cell@meta.data <- meta

object.t_cell@meta.data$plate <- object.t_cell@meta.data$split

object.t_cell@meta.data$plate[object.t_cell@meta.data$plate == "mLOY_sample1"] <- 2
object.t_cell@meta.data$plate[object.t_cell@meta.data$plate == "mLOY_sample2"] <- 1
object.t_cell@meta.data$plate[object.t_cell@meta.data$plate == "mLOY_sample3"] <- 2
object.t_cell@meta.data$plate[object.t_cell@meta.data$plate == "mLOY_sample4"] <- 1
object.t_cell@meta.data$plate[object.t_cell@meta.data$plate == "nonLOY_sample1"] <- 2
object.t_cell@meta.data$plate[object.t_cell@meta.data$plate == "nonLOY_sample2"] <- 1
object.t_cell@meta.data$plate[object.t_cell@meta.data$plate == "nonLOY_sample3"] <- 1
object.t_cell@meta.data$plate[object.t_cell@meta.data$plate == "nonLOY_sample4"] <- 2

object.t_cell <- NormalizeData(object.t_cell)
object.t_cell <- FindVariableFeatures(object.t_cell, selection.method = "vst", nfeatures = 2000)

all.gene <- rownames(object.t_cell)
object.t_cell <- ScaleData(object.t_cell)

object.t_cell <- RunPCA(object.t_cell, features = VariableFeatures(object.t_cell),
                        npcs = 50, seed.use = 123,verbose = FALSE)

# ElbowPlot(object.t_cell,ndims = 20, reduction = "pca") 
# pc = 9

object.t_cell <- object.t_cell %>%  RunHarmony(group.by.vars="plate", plot_convergence = F,project.dim = F)


## umap and tsne
object.t_cell <- RunUMAP(object.t_cell, reduction = "harmony", dims = 1:9, seed.use = 123)
object.t_cell <- RunTSNE(object.t_cell, reduction = "harmony", dims = 1:9, seed.use = 123)

DimPlot(object.t_cell, reduction = "umap",group.by = "plate")

## t cell annotation

# DimPlot(object.t_cell, reduction = "pca")

object.t_cell <- FindNeighbors(object.t_cell, reduction = "harmony", dims = 1:9, k.param = 20, annoy.metric = "euclidean")

obj <- FindClusters(object.t_cell, resolution = seq(0.4,1.2,by=0.1))

# library("clustree")
# clustree(obj)

## resolution =0.7

object.t_cell <- FindClusters(object.t_cell, resolution = 0.7, algorithm = 1, random.seed = 123)

saveRDS(object.t_cell,"object.t_cell.rds")

seurat_markers <- FindAllMarkers(object.t_cell,logfc.threshold = 0.25,min.pct = 0.5,only.pos = T)

seurat_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) -> top5

DoHeatmap(object.t_cell, features = top5$gene) + NoLegend()


##############   CellTypist 


library(Seurat)
library(data.table)
library(RColorBrewer)
library(dplyr)
library(cowplot)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(parallel)

object.t_cell <- readRDS("object.t_cell.rds")

library(reticulate)
use_python("/home/bin/python")

scanpy = import("scanpy")
celltypist = import("celltypist")
pandas <- import("pandas")
numpy = import("numpy")

adata = scanpy$AnnData(X = numpy$array(t(as.matrix(object.t_cell[['RNA']]@counts))),
                       obs = pandas$DataFrame(object.t_cell@meta.data),
                       var = pandas$DataFrame(data.frame(gene = rownames(object.t_cell[['RNA']]@counts),
                                                         row.names = rownames(object.t_cell[['RNA']]@counts))))
model = celltypist$models$Model$load(model = 'Immune_All_Low.pkl')

scanpy$pp$normalize_total(adata, target_sum=1e4)
scanpy$pp$log1p(adata)

predictions = celltypist$annotate(adata, model = 'Immune_All_Low.pkl', majority_voting = T)

object.t_cell = AddMetaData(object.t_cell, predictions$predicted_labels)

fwrite(object.t_cell@meta.data,"t.txt",row.names=T,sep = "\t",quote = F)

cell <- fread("t.txt",data.table = F)

object.t_cell@meta.data$celltype <- as.character(object.t_cell@meta.data$seurat_clusters)

object.t_cell@meta.data$celltype[object.t_cell@meta.data$celltype == "0"] = "TcmNaive_helper_T_cells"
object.t_cell@meta.data$celltype[object.t_cell@meta.data$celltype == "1"] = "TcmNaive_helper_T_cells"
object.t_cell@meta.data$celltype[object.t_cell@meta.data$celltype == "2" & 
                                   object.t_cell@meta.data$majority_voting == "MAIT cells"] = "MAIT_cells"
object.t_cell@meta.data$celltype[object.t_cell@meta.data$celltype == "2"] = "TemTemra_cytotoxic_T_cells"

object.t_cell@meta.data$celltype[object.t_cell@meta.data$celltype == "3"] = "TemTemra_cytotoxic_T_cells"
object.t_cell@meta.data$celltype[object.t_cell@meta.data$celltype == "4"] = "TemTemra_cytotoxic_T_cells"
object.t_cell@meta.data$celltype[object.t_cell@meta.data$celltype == "5"] = "TemTemra_cytotoxic_T_cells"
object.t_cell@meta.data$celltype[object.t_cell@meta.data$celltype == "6" &
                                   object.t_cell@meta.data$majority_voting == "Regulatory T cells"] = "Regulatory_T_cells"
object.t_cell@meta.data$celltype[object.t_cell@meta.data$celltype == "6"] = "TcmNaive_helper_T_cells"

object.t_cell@meta.data$celltype[object.t_cell@meta.data$celltype == "7"] = "TcmNaive_helper_T_cells"

object.t_cell@meta.data$celltype[object.t_cell@meta.data$celltype == "8"] = "TcmNaive_helper_T_cells"
object.t_cell@meta.data$celltype[object.t_cell@meta.data$celltype == "9"] = "TemTemra_cytotoxic_T_cells"

object.t_cell <- subset(object.t_cell, majority_voting != "CD16+ NK cells")

saveRDS(object.t_cell,"object.t_cell_final.rds")

library(Seurat)
library(data.table)
library(RColorBrewer)
library(dplyr)
library(cowplot)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(parallel)

object.t_cell <- readRDS("object.t_cell_final.rds")

object.t_cell$celltype <- factor(object.t_cell$celltype,levels = c("TcmNaive_helper_T_cells","TemTemra_cytotoxic_T_cells",
                                                                   "MAIT_cells","Regulatory_T_cells"))

DimPlot(object.t_cell,reduction="tsne", label = F,pt.size=0.4,group.by = "celltype",
        cols = c("#4A8260","#65A6AB","#49277A","#804825"))+
  theme(legend.position="none")

DimPlot(object.t_cell,reduction="tsne",cols = "Set1",pt.size=0.1,group.by = "condition")+
  theme(legend.position="none")
DimPlot(object.t_cell,reduction="tsne",cols = "Set1",pt.size=0.1,group.by = "split")+
  theme(legend.position="none")

Idents(object.t_cell) <- object.t_cell$celltype
seurat_markers <- FindAllMarkers(object.t_cell,logfc.threshold = 0.25,min.pct = 0.5,only.pos = T)

seurat_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) -> top5

DefaultAssay(object.t_cell)<-"RNA"
DotPlot(object.t_cell,features=unique(top5$gene),cols=c("white","red"))+
  RotatedAxis()+
  theme_test()+
  theme(axis.text=element_text(size=6,face="bold"),
        axis.title=element_blank(),legend.position="bottom",
        legend.text=element_text(size=5.5),
        legend.title=element_text(size=6,face="bold"),
        axis.text.x=element_text(angle = 90, hjust = 1, vjust = .5))


#########   GSEA  
library(Seurat)
library(UCell)
library(irGSEA)
library(doMC)

object.t_cell <- readRDS("object.t_cell_final.rds")

Idents(object.t_cell) <- object.t_cell$celltype

pbmc3k.final <- irGSEA.score(object = object.t_cell, assay = "RNA", 
                             slot = "data", seeds = 123, ncores = 10,
                             min.cells = 3, min.feature = 0,
                             custom = F, geneset = NULL, msigdb = T, 
                             species = "Homo sapiens", category = "H",  
                             subcategory = NULL, geneid = "symbol",
                             method = c("AUCell", "UCell", "singscore", 
                                        "ssgsea"),
                             aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                             kcdf = 'Gaussian')

result.dge <- irGSEA.integrate(object = pbmc3k.final, 
                               group.by = "cell.subtypes",
                               metadata = NULL, col.name = NULL,
                               method = c("AUCell","UCell","singscore",
                                          "ssgsea"))

saveRDS(result.dge,"irGSEA.rds")

# 1. Global show
irGSEA.heatmap.plot <- irGSEA.heatmap(object = result.dge, 
                                      method = "RRA",
                                      top = 50, 
                                      show.geneset = NULL)
irGSEA.heatmap.plot2 <- irGSEA.heatmap(object = result.dge, 
                                       method = "ssgsea",
                                       top = 50, 
                                       show.geneset = NULL)

irGSEA.bubble.plot <- irGSEA.bubble(object = result.dge, 
                                    method = "RRA", 
                                    top = 50)

irGSEA.upset.plot <- irGSEA.upset(object = result.dge,
                                  method = "RRA")

irGSEA.barplot.plot <- irGSEA.barplot(object = result.dge,
                                      method = c("AUCell", "UCell", "singscore",
                                                 "ssgsea"))
print(irGSEA.heatmap.plot)
print(irGSEA.heatmap.plot2)
print(irGSEA.bubble.plot)
print(irGSEA.upset.plot)
print(irGSEA.barplot.plot)

# 2. local show
# geneset <- rownames(pbmc3k.final[["UCell"]])


scatterplot <- irGSEA.density.scatterplot(object = pbmc3k.final,
                                          method = "AUCell",
                                          show.geneset = "HALLMARK-TGF-BETA-SIGNALING",
                                          reduction = "umap")

halfvlnplot <- irGSEA.halfvlnplot(object = pbmc3k.final,
                                  method = "UCell",
                                  show.geneset = "HALLMARK-TGF-BETA-SIGNALING")

ridgeplot <- irGSEA.ridgeplot(object = pbmc3k.final,
                              method = "UCell",
                              show.geneset = "HALLMARK-TGF-BETA-SIGNALING")

densityheatmap <- irGSEA.densityheatmap(object = pbmc3k.final,
                                        method = "UCell",
                                        show.geneset = "HALLMARK-TGF-BETA-SIGNALING")

print(scatterplot)
print(halfvlnplot)
print(ridgeplot)
print(densityheatmap)


############ NK cell  ########

library(Seurat)
library(data.table)
library(RColorBrewer)
library(dplyr)
library(cowplot)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(parallel)
library(harmony)

object.nk_cell <- readRDS("object.combined_V2.RDS")

object.nk_cell <- subset(object.nk_cell, anno == "NK_cell" )

count <- GetAssayData(object.nk_cell, assay = "RNA", slot = "count") # normalized data matrix

meta = object.nk_cell@meta.data[,c("nFeature_RNA", "percent.mt","S.Score","seurat_clusters",
                                   "G2M.Score","split","condition","Age","anno")] # a dataframe with rownames containing cell mata data


object.nk_cell <- CreateSeuratObject(counts = count) 
object.nk_cell@meta.data <- meta

object.nk_cell@meta.data$plate <- object.nk_cell@meta.data$split

object.nk_cell@meta.data$plate[object.nk_cell@meta.data$plate == "mLOY_sample1"] <- 2
object.nk_cell@meta.data$plate[object.nk_cell@meta.data$plate == "mLOY_sample2"] <- 1
object.nk_cell@meta.data$plate[object.nk_cell@meta.data$plate == "mLOY_sample3"] <- 2
object.nk_cell@meta.data$plate[object.nk_cell@meta.data$plate == "mLOY_sample4"] <- 1
object.nk_cell@meta.data$plate[object.nk_cell@meta.data$plate == "nonLOY_sample1"] <- 2
object.nk_cell@meta.data$plate[object.nk_cell@meta.data$plate == "nonLOY_sample2"] <- 1
object.nk_cell@meta.data$plate[object.nk_cell@meta.data$plate == "nonLOY_sample3"] <- 1
object.nk_cell@meta.data$plate[object.nk_cell@meta.data$plate == "nonLOY_sample4"] <- 2

object.nk_cell <- NormalizeData(object.nk_cell)
object.nk_cell <- FindVariableFeatures(object.nk_cell, selection.method = "vst", nfeatures = 2000)

all.gene <- rownames(object.nk_cell)

object.nk_cell <- ScaleData(object.nk_cell)
object.nk_cell <- ScaleData(object.nk_cell, features = all.gene)

object.nk_cell <- RunPCA(object.nk_cell, features = VariableFeatures(object = object.nk_cell),
                         npcs = 50, seed.use = 123,verbose = FALSE)

## 
# ElbowPlot(object.nk_cell,ndims = 20, reduction = "pca") #取碎石图拐角处的PC个数
# pc = 9

object.nk_cell <- object.nk_cell %>%  RunHarmony(group.by.vars="plate", plot_convergence = F,project.dim = F)

## umap and tsne 
object.nk_cell <- RunUMAP(object.nk_cell, reduction = "harmony", dims = 1:9, seed.use = 123)
object.nk_cell <- RunTSNE(object.nk_cell, reduction = "harmony", dims = 1:9, seed.use = 123)


# 
# DimPlot(object.nk_cell, reduction = "pca")

# object.nk_cell <- FindNeighbors(object.nk_cell, reduction = "pca", dims = 1:9, k.param = 20, annoy.metric = "euclidean")
object.nk_cell <- FindNeighbors(object.nk_cell, reduction = "harmony", dims = 1:9, k.param = 20, annoy.metric = "euclidean")

obj <- FindClusters(object.nk_cell, resolution = seq(0.4,1.2,by=0.1),random.seed = 123)

library("clustree")
clustree(obj)

##  resolution =0.8

# object.nk_cell <- FindClusters(object.nk_cell, resolution = 0.6, algorithm = 1, random.seed = 123)
object.nk_cell <- FindClusters(object.nk_cell, resolution = 0.8, algorithm = 1, random.seed = 123)
DimPlot(object.nk_cell, reduction = "umap",group.by = "plate")

saveRDS(object.nk_cell,"object.nk_cell.rds")

###  cell types annotation

library(Seurat)
library(data.table)
library(RColorBrewer)
library(dplyr)
library(cowplot)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(parallel)

object.nk_cell <- readRDS("object.nk_cell.rds")

library(reshape2)
library(ggplot2)
library(ggpubr) #用于统计分析添加统计指标
library(gridExtra) #用于组合图

# housekeeping
features = c("CD7", "NCAM1", "KLRD1", "KLRF1", "NKG7", "GNLY")
features = unique(c("NKG7",	"GZMB",	"FGFBP2",	"PRF1",	"CST7",	"GNLY",	"KLRF1",	"KLRB1",
                    "CD3D",	"CXCR4",	"CD3E",	"CD52",	"IL32",	"GNLY",	"KLRF1",	"PRF1",	
                    "CD3D",	"VIM",	"CD52",	"IL7R",	"GZMK",	"LTB",	"JUNB",	"DUSP1",
                    "JUN",	"SMC4",	"FCGR3A",	"IL2RB",	"FGFBP2",	"FCGR3A",	"IL2RB",
                    "FGFBP2",	"GNLY",	"GZMB",	"KLRF1",	"CST7",	"IGHA1",	"IGKC",	"VIM"))

## precursor NK、immature NK、(mature) NK cells

# precursor NK CD56bright CD127+
#   immature  CD56bright CD127-
#   (mature) NK cells  CD56dim CD16+ 
# CD16A (FCGR3A)
# CD56 (NCAM1)
# CD127(IL7R)
# CD94 (KLRD1)

Idents(object.nk_cell) <- object.nk_cell$seurat_clusters

features = c("IL7R","NCAM1","FCGR3A",
             "KLRD1","KLRA1","KLRC2")

cols = c("#4A8260","#A04F1E","#5D5594","#A3256D",
         "#668E2E","#BC9726","#7F6324","#535358",
         "#A0B8D0","#3D5B92","#A8C978","#4F8630",
         "#CC8382","#A02120","#D8AD64","#C36D26",
         "#B09CC3","#49277A","#804825","#6376B7",
         "#5A086B")

DimPlot(object.nk_cell,reduction="tsne", label = F,pt.size=0.4,group.by = "seurat_clusters",
        label.color="#E64B35FF",label.size=5,
        cols = c("#4A8260","#A04F1E","#5D5594","#A3256D",
                 "#668E2E","#BC9726","#7F6324","#535358"))
VlnPlot(object.nk_cell, features = features,
        cols = c("#4A8260","#A04F1E","#5D5594","#A3256D",
                                  "#668E2E","#BC9726","#7F6324","#5D5594"),pt.size = 0)


object.nk_cell@meta.data$celltype <- as.character(object.nk_cell@meta.data$seurat_clusters)
object.nk_cell@meta.data$celltype[object.nk_cell@meta.data$celltype == "0"] = "Mature_NK2"
object.nk_cell@meta.data$celltype[object.nk_cell@meta.data$celltype == "1"] = "Mature_NK2"
object.nk_cell@meta.data$celltype[object.nk_cell@meta.data$celltype == "2"] = "Mature_NK1"
object.nk_cell@meta.data$celltype[object.nk_cell@meta.data$celltype == "3"] = "Mature_NK1"
object.nk_cell@meta.data$celltype[object.nk_cell@meta.data$celltype == "4"] = "Mature_NK1"
object.nk_cell@meta.data$celltype[object.nk_cell@meta.data$celltype == "5"] = "Mature_NK1"
object.nk_cell@meta.data$celltype[object.nk_cell@meta.data$celltype == "6"] = "Mature_NK2"
object.nk_cell@meta.data$celltype[object.nk_cell@meta.data$celltype == "7"] = "Mature_NK1"

saveRDS(object.nk_cell,"object.nk_cell.rds")

#########   GSEA  
library(Seurat)
library(UCell)
library(irGSEA)
library(doMC)
library(data.table)

BiocManager::install("SeuratData")

object.nk_cell <- readRDS("object.nk_cell.rds")

Idents(object.nk_cell) <- object.nk_cell$celltype

pbmc3k.final <- irGSEA.score(object = object.nk_cell, assay = "RNA", 
                             slot = "data", seeds = 123, ncores = 10,
                             min.cells = 3, min.feature = 0,
                             custom = F, geneset = NULL, msigdb = T, 
                             species = "Homo sapiens", category = "H",  
                             subcategory = NULL, geneid = "symbol",
                             method = c("AUCell", "UCell", "singscore", 
                                        "ssgsea"),
                             aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                             kcdf = 'Gaussian')

result.dge <- irGSEA.integrate(object = pbmc3k.final, 
                               group.by = "cell.subtypes",
                               metadata = NULL, col.name = NULL,
                               method = c("AUCell","UCell","singscore",
                                          "ssgsea"))

saveRDS(result.dge,"irGSEA.rds")


######  myeloid 

library(Seurat)
library(data.table)
library(RColorBrewer)
library(dplyr)
library(cowplot)
library(tidyr)
library(ggplot2)
library(patchwork)
library(parallel)
library(harmony)

 #########   myeloid  cell clusters   ############

object.mye_cell <- readRDS("object.combined_V2.RDS")

object.mye_cell <- subset(object.mye_cell, anno == "Myeloid_cell")

count <- GetAssayData(object.mye_cell, assay = "RNA", slot = "count") # normalized data matrix

meta = object.mye_cell@meta.data[,c("nFeature_RNA", "percent.mt","S.Score", 
                                    "G2M.Score","split","condition","Age","anno")] # a dataframe with rownames containing cell mata data

object.mye_cell <- CreateSeuratObject(
  counts = count) 
object.mye_cell@meta.data <- meta

object.mye_cell@meta.data$plate <- object.mye_cell@meta.data$split

object.mye_cell@meta.data$plate[object.mye_cell@meta.data$plate == "mLOY_sample1"] <- 2
object.mye_cell@meta.data$plate[object.mye_cell@meta.data$plate == "mLOY_sample2"] <- 1
object.mye_cell@meta.data$plate[object.mye_cell@meta.data$plate == "mLOY_sample3"] <- 2
object.mye_cell@meta.data$plate[object.mye_cell@meta.data$plate == "mLOY_sample4"] <- 1
object.mye_cell@meta.data$plate[object.mye_cell@meta.data$plate == "nonLOY_sample1"] <- 2
object.mye_cell@meta.data$plate[object.mye_cell@meta.data$plate == "nonLOY_sample2"] <- 1
object.mye_cell@meta.data$plate[object.mye_cell@meta.data$plate == "nonLOY_sample3"] <- 1
object.mye_cell@meta.data$plate[object.mye_cell@meta.data$plate == "nonLOY_sample4"] <- 2

object.mye_cell <- NormalizeData(object.mye_cell)
object.mye_cell <- FindVariableFeatures(object.mye_cell, selection.method = "vst", nfeatures = 2000)

# object.mye_cell <- ScaleData(object.mye_cell)
all.gene <- rownames(object.mye_cell)
object.mye_cell <- ScaleData(object.mye_cell)
object.mye_cell <- RunPCA(object.mye_cell, features = VariableFeatures(object = object.mye_cell),
                          npcs = 50, seed.use = 123,verbose = FALSE)

# ElbowPlot(object.mye_cell,ndims = 20, reduction = "pca") #取碎石图拐角处的PC个数
# pc = 9

object.mye_cell <- object.mye_cell %>%  RunHarmony(group.by.vars="plate", plot_convergence = F,project.dim = F)

## umap and tsne 
object.mye_cell <- RunUMAP(object.mye_cell, reduction = "harmony", dims = 1:9, seed.use = 123)
object.mye_cell <- RunTSNE(object.mye_cell, reduction = "harmony", dims = 1:9)

# DimPlot(object.mye_cell, reduction = "pca")

object.mye_cell <- FindNeighbors(object.mye_cell, reduction = "harmony", dims = 1:9, k.param = 20, annoy.metric = "euclidean")

obj <- FindClusters(object.mye_cell, resolution = seq(0.4,1.2,by=0.1))

# library("clustree")
# clustree(obj)

##  resolution =0.7

object.mye_cell <- FindClusters(object.mye_cell, resolution = 0.7, algorithm = 1, random.seed = 123)

DimPlot(object.mye_cell, reduction = "umap",group.by = "plate")

saveRDS(object.mye_cell,"object.mye_cell.rds")


##############    CellTypist 


library(Seurat)
library(data.table)
library(RColorBrewer)
library(dplyr)
library(cowplot)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(parallel)

object.mye_cell <- readRDS("object.mye_cell.rds")

library(reticulate)
use_python("home/bin/python")

scanpy = import("scanpy")
celltypist = import("celltypist")
pandas <- import("pandas")
numpy = import("numpy")

adata = scanpy$AnnData(X = numpy$array(t(as.matrix(object.mye_cell[['RNA']]@counts))),
                       obs = pandas$DataFrame(object.mye_cell@meta.data),
                       var = pandas$DataFrame(data.frame(gene = rownames(object.mye_cell[['RNA']]@counts),
                                                         row.names = rownames(object.mye_cell[['RNA']]@counts))))
model = celltypist$models$Model$load(model = 'Immune_All_Low.pkl')
scanpy$pp$normalize_total(adata, target_sum=1e4)
scanpy$pp$log1p(adata)

predictions = celltypist$annotate(adata, model = 'Immune_All_Low.pkl', majority_voting = T)

object.mye_cell = AddMetaData(object.mye_cell, predictions$predicted_labels)


fwrite(object.mye_cell@meta.data,"mye.txt",row.names=T,sep = "\t",quote = F)

cell <- fread("mye.txt",data.table = F)


object.mye_cell@meta.data$celltype <- as.character(object.mye_cell@meta.data$seurat_clusters)

object.mye_cell@meta.data$celltype[object.mye_cell@meta.data$celltype == "0"] = "Classical_monocytes"
object.mye_cell@meta.data$celltype[object.mye_cell@meta.data$celltype == "1"] = "Classical_monocytes"
object.mye_cell@meta.data$celltype[object.mye_cell@meta.data$celltype == "2"] = "NonClassical_monocytes"
object.mye_cell@meta.data$celltype[object.mye_cell@meta.data$celltype == "3"] = "Classical_monocytes"
object.mye_cell@meta.data$celltype[object.mye_cell@meta.data$celltype == "4"] = "NonClassical_monocytes"
object.mye_cell@meta.data$celltype[object.mye_cell@meta.data$celltype == "5"] = "Classical_monocytes"
object.mye_cell@meta.data$celltype[object.mye_cell@meta.data$celltype == "6"] = "DC2"
object.mye_cell@meta.data$celltype[object.mye_cell@meta.data$celltype == "7"] = "TemTemra_cytotoxic_T_cells"
object.mye_cell@meta.data$celltype[object.mye_cell@meta.data$celltype == "8"] = "pDC cells"
object.mye_cell@meta.data$celltype[object.mye_cell@meta.data$celltype == "9"] = "Classical_monocytes"

object.mye_cell <- subset(object.mye_cell, seurat_clusters != "7")

saveRDS(object.mye_cell,"object.mye_cell_final.rds")


library(Seurat)
library(data.table)
library(RColorBrewer)
library(dplyr)
library(cowplot)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(parallel)

object.mye_cell <- readRDS("object.mye_cell_final.rds")

DimPlot(object.mye_cell,reduction="tsne", label = F,pt.size=0.4,group.by = "celltype",
        label.color="#E64B35FF",label.size=5,
        cols = c("#A0B8D0","#D8AD64","#BB533F","#4F8630"))+
  theme(legend.position="none")

DimPlot(object.mye_cell,reduction="tsne",cols = "Set1",pt.size=0.1,group.by = "condition")+
  theme(legend.position="none")
DimPlot(object.mye_cell,reduction="tsne",cols = "Set1",pt.size=0.1,group.by = "split")+
  theme(legend.position="none")




Idents(object.mye_cell) <- object.mye_cell$celltype
seurat_markers <- FindAllMarkers(object.mye_cell,logfc.threshold = 0.25,min.pct = 0.5,only.pos = T)

seurat_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) -> top5

DefaultAssay(object.mye_cell)<-"RNA"
DotPlot(object.mye_cell,features=unique(top5$gene),cols=c("white","red"))+
  RotatedAxis()+
  theme_test()+
  theme(axis.text=element_text(size=6,face="bold"),
        axis.title=element_blank(),legend.position="bottom",
        legend.text=element_text(size=5.5),
        legend.title=element_text(size=6,face="bold"),
        axis.text.x=element_text(angle = 90, hjust = 1, vjust = .5))



#######   scDC  ##########
library("Seurat")
library(scDC)
library(ggplot2)
library(dplyr)
library(data.table)

object.combined_left <- readRDS("object.combined_V2.RDS")

object.combined_left  <- subset(object.combined_left, anno2 != "Myeloid_cell"  & anno2 != "T_cell" &
                                 anno2 != "Red_blood_cell" & anno2 != "Megakaryocyte"  )
object.combined_left$anno2 <- factor(object.combined_left$anno2,
                                     levels = c("TcmNaive_helper_T_cells","TemTemra_cytotoxic_T_cells",
                                                "MAIT_cells","Regulatory_T_cells",
                                                "Mature_NK1","Mature_NK2",
                                                "Classical_monocytes","NonClassical_monocytes",
                                                "DC2","pDC cells","B_cell"))



subject <- object.combined_left@meta.data$split
cellTypes <- object.combined_left@meta.data$anno2

sample <- data.frame(unique(object.combined_left@meta.data[,c("split","condition")]))
rownames(sample) <- sample$split
colnames(sample)[1] <- "subject"
# 
res_scDC_noClust <- scDC_noClustering(cellTypes, subject, calCI = TRUE,
                                      calCI_method = "BCa",ncores = 8,verbose = F,
                                      nboot = 10000)
saveRDS(res_scDC_noClust,"scDC_subject2.rds")

res_scDC_noClust <- readRDS("scDC_subject2.rds")

df_toPlot <- res_scDC_noClust$results
df_toPlot$median <- apply(res_scDC_noClust$thetastar, 1, median)

df_toPlot <- left_join(df_toPlot,sample)
# df_toPlot$method <- factor(df_toPlot$method, levels = c("BCa", "percentile", "multinom"))
n_method <- length(unique(df_toPlot$method))

n_celltype = length(unique(df_toPlot$cellTypes))
              
g_bar <- ggplot2::ggplot(df_toPlot, aes(x = subject, y = median, fill = condition)) +
  ggplot2::geom_bar(stat="identity", position = "dodge", alpha = 0.8) +
  ggplot2::theme_bw() +
  ggplot2::scale_fill_brewer(palette = "Set2") +
  ggplot2::ylab("Proportion") +
  ggplot2::geom_errorbar(aes(ymin=conf_low, ymax=conf_high, color = method), width=.3,lwd = 1,
                         position=position_dodge(width = 0.5)) +
  ggplot2::theme(axis.text.x = element_text(angle = 90), text = element_text(size = 12)) +
  ggplot2::scale_color_manual(values = "black") +
  ggplot2::facet_wrap(~cellTypes, ncol = n_celltype,
                      labeller = labeller(cellTypes  = label_wrap_gen(width = 10,  multi_line = TRUE))) +
  ggplot2::coord_flip()+
  ggplot2::ylim(c(0,1))
  
  
print(g_bar)


df_toPlot2 <- df_toPlot[,c("subject","cellTypes","median")]

p <- ggplot(df_toPlot2, aes(x=subject, y=median,fill=cellTypes)) + 
  geom_bar(stat="identity",position = "stack")+ 
  scale_fill_brewer(palette = "Set2")+      
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        strip.background = element_blank(),panel.border = element_blank(),
        axis.line=element_line())+
  # labs(title="housekeeping gene") +
  ylab("Ratio of cell types")
print(p)


cellTypes <- object.combined_left@meta.data$anno2

cond <- object.combined_left@meta.data$condition

# 
res_scDC_noClust_cond <- scDC_noClustering(cellTypes, cond, calCI = TRUE,
                                      calCI_method = "BCa",ncores = 8,verbose = F,
                                      nboot = 10000)
# saveRDS(res_scDC_noClust_cond,"scDC_cond_main.rds")
saveRDS(res_scDC_noClust_cond,"scDC_cond2.rds")

res_scDC_noClust_cond <- readRDS("scDC_cond2.rds")
# res_scDC_noClust_cond <- readRDS("scDC_cond_main.rds")

df_toPlot_cond <- res_scDC_noClust_cond$results
df_toPlot_cond$median <- apply(res_scDC_noClust_cond$thetastar, 1, median)
df_toPlot_cond$conf_low <- res_scDC_noClust_cond$results$conf_low[,1]
df_toPlot_cond$conf_high <- res_scDC_noClust_cond$results$conf_high[,1]

df_toPlot_cond$subject <- factor(df_toPlot_cond$subject,levels = c("nonmLOY","mLOY"))

df_toPlot_cond$conf_low <- df_toPlot_cond$conf_low * 100
df_toPlot_cond$conf_high <- df_toPlot_cond$conf_high * 100
df_toPlot_cond$median <- df_toPlot_cond$median * 100

df_toPlot_cond$cellTypes <- factor(df_toPlot_cond$cellTypes,levels =c("TcmNaive_helper_T_cells","TemTemra_cytotoxic_T_cells",
                                                            "MAIT_cells","Regulatory_T_cells",
                                                            "Mature_NK1","Mature_NK2",
                                                            "Classical_monocytes","NonClassical_monocytes",
                                                            "DC2","pDC cells","B_cell"))
library(ggplot2)
library(parallel)
# mean
pp1 <- ggplot(df_toPlot_cond,aes(x=cellTypes,y=median))+
  geom_bar(mapping=aes(fill = subject),position=position_dodge(0.6),width = 0.5,stat = "identity")+
  geom_errorbar(aes(ymin=conf_low,ymax=conf_high,fill = subject),position=position_dodge(0.6),width=0.2)+
  # geom_text(aes(x=factor(cellTypes),y=conf_high+0.2),
  #           size=3,position= position_dodge(0.6))+
  labs(x = "",y = "Nutrient (mg/L)")+
  scale_y_continuous(breaks = seq(0,100,10))+
  theme_classic()+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25), 
        axis.line=element_line(colour="black",size=0.25), 
        axis.title=element_text(size=13,color="black"),
        axis.text = element_text(size=12,color="black"), 
        legend.position="none"
        )
print(pp1)

pp <- function(i,df_toPlot_cond) {
  
  cellTypes = as.character(unique(df_toPlot_cond$cellTypes))[i]
  d1 <- df_toPlot_cond[df_toPlot_cond$cellTypes == cellTypes,]
  d1 <- d1[order(d1$subject,decreasing=F),]
  
  ymax <- max(ceiling(d1$conf_high))+1
  ymin <- min(floor(d1$conf_low))-1
  
  plot<-ggplot()+
          geom_ribbon(data=d1, aes(x=subject, y=median, group = 1,
                                   ymin=conf_low,ymax=conf_high),fill='#3E86B5', alpha=0.3)+
          geom_line(data=d1, aes(x=subject, y=median, group = 1),color='#3E86B5')+theme_bw()+
          xlab("Conditions")+ylab("Percentage(Booststrap)(%)")+
          scale_y_continuous(breaks = seq(0,100,1),limits = c(ymin,ymax))+
          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                strip.background = element_blank(),panel.border = element_blank(),
                axis.line=element_line())+
          labs(title=cellTypes)
  
  plot
}

P_plot <- mclapply(1:length(as.character(unique(df_toPlot_cond$cellTypes))), function(i) try(pp(i,df_toPlot_cond), TRUE))

plot1 <- P_plot[[1]]
plot1


#### TF DoRothEA   ####

library(Seurat)
library(dorothea)
library(tidyverse)
library(data.table)
object.combined_left <- readRDS("object.combined_V2.RDS")

object.combined_left  <- subset(object.combined_left, anno2 != "Myeloid_cell"  & anno2 != "T_cell" &
                                  anno2 != "Red_blood_cell" & anno2 != "Megakaryocyte"  )

object.combined_left$anno2 <- factor(object.combined_left$anno2,
                                     levels = c("TcmNaive_helper_T_cells","TemTemra_cytotoxic_T_cells",
                                                "MAIT_cells","Regulatory_T_cells",
                                                "Mature_NK1","Mature_NK2",
                                                "Classical_monocytes","NonClassical_monocytes",
                                                "DC2","pDC cells","B_cell"))

dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))
regulon <- dorothea_regulon_human %>%
  dplyr::filter(confidence %in% c("A","B","C"))

sce <- run_viper(object.combined_left, regulon,
                 options = list(method = "scale", minsize = 4,
                                eset.filter = FALSE, cores = 8,
                                verbose = FALSE))

DefaultAssay(object = sce) <- "dorothea"
Idents(sce) <- sce$anno2

library(future)
# check the current active plan
plan()
plan("multiprocess", workers = 4)
plan()


sce.markers <- FindAllMarkers(object = sce, 
                              only.pos = TRUE, 
                              min.pct = 0.25, 
                              thresh.use = 0.25)
DT::datatable(sce.markers)
pro='dorothea-markers-for-pbmc3k'
write.csv(sce.markers,file=paste0(pro,'_sce.markers.csv'))
save(sce.markers,file = paste0(pro, '_sce.markers.Rdata'))

## We transform Viper scores, scaled by seurat, into a data frame to better 
## handling the results
viper_scores_df <- GetAssayData(sce, slot = "data", 
                                assay = "dorothea") %>%
  data.frame(check.names = F) %>%
  t()
# viper_scores_df[1:4,1:4]


## We create a data frame containing the cells and their clusters
CellsClusters <- data.frame(cell = rownames(object.combined_left@meta.data), 
                            cell_type = as.character(object.combined_left@meta.data$anno2),
                            condition = as.character(object.combined_left@meta.data$condition),
                            check.names = F)
head(CellsClusters)

## We create a data frame with the Viper score per cell and its clusters
viper_scores_clusters <- viper_scores_df  %>%
  data.frame() %>% 
  rownames_to_column("cell") %>%
  gather(tf, activity, -cell) %>%
  inner_join(CellsClusters)

saveRDS(viper_scores_clusters,"combind_viper_scores_df.rds")

# viper_scores_clusters <- left_join(viper_scores_clusters,CellsClusters)
viper_scores_clusters <- readRDS("combind_viper_scores_df.rds")

viper_scores_clusters <- viper_scores_clusters[,c("cell","tf","activity","cell_type")]

## We summarize the Viper scores by cellpopulation
summarized_viper_scores <- viper_scores_clusters %>% 
  group_by(tf, cell_type) %>%
  summarise(avg = mean(activity),
            std = sd(activity))
# For visualization purposes, we select the 20 most variable TFs across clusters according to our scores.
head(summarized_viper_scores)

## We select the 20 most variable TFs. (20*11 populations = 220)
highly_variable_tfs <- summarized_viper_scores %>%
  group_by(tf) %>%
  mutate(var = var(avg))  %>%
  ungroup() %>%
  top_n(220, var) %>%
  distinct(tf)
highly_variable_tfs

## We prepare the data for the plot
summarized_viper_scores_df <- summarized_viper_scores %>%
  semi_join(highly_variable_tfs, by = "tf") %>%
  dplyr::select(-std) %>%   
  spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE) 
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)
colnames(summarized_viper_scores_df)


my_breaks <- c(seq(min(summarized_viper_scores_df), 0, 
                   length.out=ceiling(palette_length/2) + 1),
               seq(max(summarized_viper_scores_df)/palette_length, 
                   max(summarized_viper_scores_df), 
                   length.out=floor(palette_length/2)))
library(pheatmap)
viper_hmap <- pheatmap(t(summarized_viper_scores_df),fontsize=14, 
                       fontsize_row = 10, 
                       color=my_color, breaks = my_breaks, 
                       main = "DoRothEA (ABC)", angle_col = 45,
                       treeheight_col = 0,  border_color = NA) 

print(viper_hmap)


library(parallel)
library(data.table)
library(dplyr)
viper_scores_clusters <- readRDS("combind_viper_scores_df.rds")

viper_scores_clusters$type <- paste0(viper_scores_clusters$tf,"_",viper_scores_clusters$cell_type)
# viper_scores_clusters$acti_zscore <- scale(viper_scores_clusters$activity)



#########  irGSEA hallmark GSEA   ############
library(Seurat)
# library(SeuratData)
library(UCell)
library(irGSEA)
library(doMC)

object.combined_left <- readRDS("object.combined_V2.RDS")

object.combined_left  <- subset(object.combined_left, anno2 != "Myeloid_cell"  & anno2 != "T_cell" &
                                  anno2 != "Red_blood_cell" & anno2 != "Megakaryocyte"  )

object.combined_left$anno2 <- factor(object.combined_left$anno2,
                                     levels = c("TcmNaive_helper_T_cells","TemTemra_cytotoxic_T_cells",
                                                "MAIT_cells","Regulatory_T_cells",
                                                "Mature_NK1","Mature_NK2",
                                                "Classical_monocytes","NonClassical_monocytes",
                                                "DC2","pDC cells","B_cell"))
Idents(object.combined_left) <- object.combined_left$anno2

pbmc3k.final <- irGSEA.score(object = object.combined_left, assay = "RNA", 
                             slot = "data", seeds = 123, ncores = 10,
                             min.cells = 3, min.feature = 0,
                             custom = F, geneset = NULL, msigdb = T, 
                             species = "Homo sapiens", category = "H",  
                             subcategory = NULL, geneid = "symbol",
                             method = c("AUCell", "UCell", "singscore", 
                                        "ssgsea"),
                             aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                             kcdf = 'Gaussian')

object <- pbmc3k.final

object <- SeuratObject::SetIdent(object, value = "anno2")
anno.ident <- SeuratObject::Idents(object)
anno.ident <- as.factor(as.character(anno.ident))

method = c("AUCell","UCell","singscore","ssgsea")
# group1 <- levels(Idents(object))
deg.geneset <- list()
avg_diff <- NULL
for (i in seq_along(method)) {
  message(paste0("Calculate differential gene set", " : ", method[i]))
  
  marker.geneset <- lapply(levels(anno.ident), function(x){
    a <- Seurat::FindMarkers(object = object,
                             assay = method[i],
                             slot = "scale.data",
                             subset.ident = x,
                             group.by	= "condition",
                             
                             ident.1 = "mLOY",
                             ident.2 = "nonmLOY",
                             test.use = "wilcox",
                             min.pct = -Inf,
                             logfc.threshold = 0,
                             min.cells.group = 0,
                             min.diff.pct = -Inf,
                             verbose = F,
                             min.cells.feature = 0)
    a <- a %>% tibble::rownames_to_column(var = "gene") %>%
      dplyr::mutate(cluster = x, direction = dplyr::if_else(avg_diff >0, "up", "down")) %>%
      dplyr::select(-c("pct.1", "pct.2"))
  })
  marker.geneset <- do.call(rbind, marker.geneset)
  deg.geneset[[i]] <- cbind(marker.geneset, methods = method[i])
  
}
names(deg.geneset) <- method
deg.geneset.list <- deg.geneset %>% purrr::map( ~.x %>% dplyr::rename(Name = gene))
deg.geneset <- do.call(rbind, deg.geneset)
p_val_adj <- NULL
cluster <- NULL
methods <- NULL
deg.cluster <- deg.geneset %>%
  dplyr::filter(p_val_adj <= 0.05) %>%
  dplyr::select(c("avg_diff", "cluster", "gene","methods"))

deg.cluster$cluster <- as.factor(as.character(deg.cluster$cluster))
deg.cluster <- deg.cluster %>%
  dplyr::group_split(cluster) %>%
  purrr::set_names(levels(deg.cluster$cluster))

deg.cluster.postive <- lapply(deg.cluster, function(x){
  a <- x %>%
    dplyr::filter(avg_diff > 0) %>%
    dplyr::arrange(methods, dplyr::desc(avg_diff))
  a$methods <- as.factor(as.character(a$methods))
  b <- a %>%
    dplyr::group_split(methods) %>%
    purrr::set_names(levels(a$methods))
})

deg.cluster.negative <- lapply(deg.cluster, function(x){
  a <- x %>%
    dplyr::filter(avg_diff < 0) %>%
    dplyr::arrange(methods, avg_diff)
  a$methods <- as.factor(as.character(a$methods))
  b <- a %>%
    dplyr::group_split(methods) %>%
    purrr::set_names(levels(a$methods))
  
})

if (! identical(names(deg.cluster.postive), names(deg.cluster.postive %>% purrr::compact()))) {
  a <- setdiff(names(deg.cluster.postive), names(deg.cluster.postive %>% purrr::compact()))
  a <- stringr::str_c(a, collapse = ", ")
  message(paste0("No sigficant genesets in cluster : ", a, " after wilicox test"))
  deg.cluster.postive <- deg.cluster.postive %>% purrr::compact()
}
if (! identical(names(deg.cluster.negative), names(deg.cluster.negative %>% purrr::compact()))) {
  a <- setdiff(names(deg.cluster.negative), names(deg.cluster.negative %>% purrr::compact()))
  a <- stringr::str_c(a, collapse = ", ")
  message(paste0("No sigficant genesets in cluster : ", a, " after wilicox test"))
  deg.cluster.negative <- deg.cluster.negative %>% purrr::compact()
}

# perform RRA analysis
sig.genesets.postive <- lapply(deg.cluster.postive, function(x){
  x <- lapply(x,function(y){y$gene})
  x <- RobustRankAggreg::aggregateRanks(glist = x, N = length(unique(unlist(x))))
})
sig.genesets.negative <- lapply(deg.cluster.negative,function(x){
  x <- lapply(x,function(y){y$gene})
  x <- RobustRankAggreg::aggregateRanks(glist = x, N = length(unique(unlist(x))))
})
Name <- NULL
sig.genesets.postive <- sig.genesets.postive %>%
  reshape2::melt(id.vars = "Name") %>%
  dplyr::mutate(direction = "up")
sig.genesets.negative <- sig.genesets.negative %>%
  reshape2::melt(id.vars = "Name") %>%
  dplyr::mutate(direction = "down")
sig.genesets <- rbind(sig.genesets.postive, sig.genesets.negative)
colnames(sig.genesets)[4] <- "cluster"
value <- NULL
sig.genesets <- sig.genesets %>%
  dplyr::select(c("Name", "value","cluster","direction")) %>%
  dplyr::rename(pvalue = value) %>%
  dplyr::mutate(method = "RRA")
deg.geneset.list[["RRA"]] <- sig.genesets


result.dge <- deg.geneset.list

saveRDS(pbmc3k.final,"irGSEA.rds")


  
  
#######  cellchat  cell-cell communication  #######

### cellchat soft 
library(Seurat)
library(CellChat)
library(patchwork)
library(data.table)

object.combined_left <- readRDS("object.combined_V2.RDS")

object.combined_left  <- subset(object.combined_left, anno2 != "Myeloid_cell"  & anno2 != "T_cell" &
                                  anno2 != "Red_blood_cell" & anno2 != "Megakaryocyte"  )

## 3 v 3
object.combined_left  <- subset(object.combined_left, split != "mLOY_sample4"  & split != "nonLOY_sample1")


object.combined_left$anno2 <- factor(object.combined_left$anno2,
                                     levels = c("TcmNaive_helper_T_cells","TemTemra_cytotoxic_T_cells",
                                                "MAIT_cells","Regulatory_T_cells",
                                                "Mature_NK1","Mature_NK2",
                                                "Classical_monocytes","NonClassical_monocytes",
                                                "DC2","pDC cells","B_cell"))

data.input = GetAssayData(object.combined_left, assay = "RNA", slot = "data") # normalized data matrix
meta = object.combined_left@meta.data # a dataframe with rownames containing cell mata data

cell.LOY = rownames(meta)[meta$condition == "mLOY"] # extract the cell names from disease data
cell.nonLOY = rownames(meta)[meta$condition == "nonmLOY"] # extract the cell names from disease data

data.LOY = data.input[, cell.LOY]
meta.LOY = meta[cell.LOY, ]

data.nonLOY = data.input[, cell.nonLOY]
meta.nonLOY = meta[cell.nonLOY, ]

## LOY samples cell-cell communication
cellchat <- createCellChat(object = data.LOY, meta = meta.LOY, group.by = "anno2")

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
# subset the expression data of signaling genes for saving computation cost
cellchat@DB <- CellChatDB

cellchat <- subsetData(cellchat)


# future::plan("multiprocess", workers = 4) # do parallel
#> Warning: Strategy 'multiprocess' is deprecated in future (>= 1.20.0). Instead,
#> explicitly specify either 'multisession' or 'multicore'. In the current R
#> session, 'multiprocess' equals 'multisession'.
#> Warning in supportsMulticoreAndRStudio(...): [ONE-TIME WARNING] Forked
#> processing ('multicore') is not supported when running R from RStudio
#> because it is considered unstable. For more details, how to control forked
#> processing or not, and how to silence this warning in future R sessions, see ?
#> parallelly::supportsMulticore
future::plan("multiprocess", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)


cellchat <- computeCommunProb(cellchat)
#> triMean is used for calculating the average gene expression per cell group. 
#> [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2022-11-26 08:34:38]"
#> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2022-11-26 08:35:36]"
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

# groupSize <- as.numeric(table(cellchat@idents))
# par(mfrow = c(1,2), xpd=TRUE)
# netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
# netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
# netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

cellchat.LOY <- cellchat






### nonLOY samples cell-cell communication
cellchat <- createCellChat(object = data.nonLOY, meta = meta.nonLOY, group.by = "anno2")

cellchat <- addMeta(cellchat, meta = meta.nonLOY)
cellchat <- setIdent(cellchat, ident.use = "anno2") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
# subset the expression data of signaling genes for saving computation cost
cellchat@DB <- CellChatDB

cellchat <- subsetData(cellchat)


future::plan("multiprocess", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)


cellchat <- computeCommunProb(cellchat)

cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
# netAnalysis_signalingRole_network(cellchat.LOY, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)


cellchat.nonLOY <- cellchat

saveRDS(cellchat.nonLOY,"cellchat.nonLOY.rds")
saveRDS(cellchat.LOY,"cellchat.LOY.rds")

saveRDS(cellchat.nonLOY,"cellchat.nonLOY_3v3.rds")
saveRDS(cellchat.LOY,"cellchat.LOY_3v3.rds")

### LOY vs nonLOY diff cell-cell communication
library(Seurat)
library(CellChat)
library(patchwork)

# object.combined_left <- readRDS("object.combined_V2.RDS")

cellchat.nonLOY <- readRDS("cellchat.nonLOY.rds")
cellchat.LOY <- readRDS("cellchat.LOY.rds")

cellchat.nonLOY <- readRDS("cellchat.nonLOY_3v3.rds")
cellchat.LOY <- readRDS("cellchat.LOY_3v3.rds")


object.list <- list( nonLOY = cellchat.nonLOY, LOY = cellchat.LOY)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

setwd("/home/wangjunhao/test/")
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

# levels(object.list$nonLOY@idents)

cellchat@netP$nonLOY$pathways
cellchat@netP$LOY$pathways


num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)

library(reticulate)
use_python('home/bin/python3',required = T)

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat,type = "functional")

# cellchat <- netEmbedding(cellchat, umap.method = 'uwot',type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2

# netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)


rankSimilarity(cellchat, type = "functional")
#> Compute the distance of signaling networks between datasets 1 2

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2


library(ComplexHeatmap)

i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))



target <-  c("TcmNaive_helper_T_cells","TemTemra_cytotoxic_T_cells",
             "MAIT_cells","Regulatory_T_cells")

library(RColorBrewer)
Incre_sig <- c()
Decre_sig <- c()
Incre_path <- c()
Decre_path <- c()

for (targets.use in target) {
  # targets.use <- target[1]

  p <- netVisual_bubble(cellchat, sources.use = c(1:11), targets.use = targets.use,  
                        color.heatmap="Spectral",n.colors = 3,comparison = c(1, 2), angle.x = 45)
  
  print(p)

  gg1 <- netVisual_bubble(cellchat, sources.use = c(1:11), targets.use = targets.use,  
                          color.heatmap="Spectral",n.colors = 3,comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in LOY", angle.x = 45, remove.isolate = T)
  #> Comparing communications on a merged object
  gg2 <- netVisual_bubble(cellchat, sources.use = c(1:11), targets.use = targets.use,  
                          color.heatmap="Spectral",n.colors = 3,comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in LOY", angle.x = 45, remove.isolate = T)
  #> Comparing communications on a merged object
  print(gg1 + gg2)
  Incre_sig <- c(Incre_sig, unique(as.character(gg1$data$interaction_name_2)))
  Decre_sig <- c(Decre_sig, unique(as.character(gg2$data$interaction_name_2)))
  Incre_path <- c(Incre_path, unique(as.character(gg1$data$pathway_name)))
  Decre_path <- c(Decre_path, unique(as.character(gg2$data$pathway_name)))
}




# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "LOY"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LOY
net.up <- subsetCommunication(cellchat, net = net, datasets = "LOY",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in nonLOY, i.e.,downregulated in LOY
net.down <- subsetCommunication(cellchat, net = net, datasets = "nonLOY",ligand.logFC = -0.1, receptor.logFC = -0.1)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use =c(1:9), targets.use = 4, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use =c(1:9), targets.use = 4, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2


# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = 4, targets.use = c(1:9), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Note: The first link end is drawn out of sector 'MIF'.
netVisual_chord_gene(object.list[[1]], sources.use = 4, targets.use = c(1:9), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]),legend.pos.y = 30)
# FCER2A_ITGAM_ITGB2
# FCER2A_ITGAX_ITGB2

netVisual_chord_gene(object.list[[1]], sources.use = 4, targets.use = c(1:9), net = net.up,lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(object.list[[1]], sources.use = 4, targets.use = c(1:9),slot.name = 'net', lab.cex = 0.5,small.gap = 3.5,legend.pos.y = 30)



pathways.show <- Incre_path[1]
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.



cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("nonLOY", "LOY")) # set factor level
plotGeneExpression(cellchat, signaling = Incre_path[1], split.by = "datasets", colors.ggplot = T)

saveRDS(cellchat, file = "cellchat_LOY_vs_nonLOY.rds")


#### useful result
# 1 Pathway distance
t1 <- rankSimilarity(cellchat, type = "functional")
pathway1 <- as.character(t1$data$name)

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
pathway2 <- unique(as.character(gg1$data$name))

cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("nonLOY", "LOY")) # set factor level

for (num in 1:length(pathway1)) {
  # num =1
  weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathway1[num]) # control the edge weights across different datasets
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(object.list)) {
    netVisual_aggregate(object.list[[i]], signaling = pathway1[num], layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathway1[num], names(object.list)[i]))
  }
  par(mfrow = c(1,2), xpd=TRUE)
  ht <- list()
  for (i in 1:length(object.list)) {
    ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling =  pathway1[num], color.heatmap = "Reds",title.name = paste( pathway1[num], "signaling ",names(object.list)[i]))
  }
  #> Do heatmap based on a single object 
  #> 
  #> Do heatmap based on a single object
  ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
  
  print(plotGeneExpression(cellchat, signaling = pathway1[num], split.by = "datasets", colors.ggplot = T))
  
}



pathway3 <- c("BTLA", "CD23","TNF")

for (num in 1:length(pathway3)) {
  netVisual_aggregate(cellchat.LOY, signaling = pathway3[num], layout = "circle")
}


