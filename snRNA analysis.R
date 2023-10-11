

##  public snRNA-seq PD data

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
require(ensembldb)

setwd("/GSE184950/")

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
  x <- subset(x,  nFeature_RNA > 200 & nFeature_RNA < 2500)
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

object.combined@meta.data -> df
df$cell <- rownames(df)
sample <- fread("/GSE184950_sample.csv",data.table = F)
colnames(sample)[1] <- "orig.ident"

df <- left_join(df,sample,by="orig.ident")

object.combined@meta.data <- df

all.genes <- rownames(object.combined@assays$RNA@counts)
object.combined <- ScaleData(object.combined, features = all.genes,
                             vars.to.regress = c("nFeature_RNA", "Age","Gender", "postmortem_interval_hours", "orig.ident"))

object.combined2 <- NormalizeData(object.combined)
object.combined2 <- FindVariableFeatures(object.combined2, selection.method = "vst", nfeatures = 2000)

object.combined2 <- RunPCA(object.combined2,
                          npcs = 50, seed.use = 123,verbose = FALSE)

object.combined2 <- object.combined2 %>%  RunHarmony(group.by.vars="orig.ident", plot_convergence = F,project.dim = F)
object.combined2 <- RunUMAP(object.combined2, reduction = "harmony", dims = 1:20, seed.use = 123)
object.combined2 <- RunTSNE(object.combined2, reduction = "harmony", dims = 1:20, seed.use = 123)

object.combined2@assays$RNA$scale.data <- matrix(c(1:4),ncol=2)

object.combined2 <- FindNeighbors(object.combined2, reduction = "harmony", dims = 1:20, k.param = 20, annoy.metric = "euclidean")

object.combined2 <- FindClusters(object.combined2, resolution = 0.2, algorithm = 1, random.seed = 123)

DimPlot(object.combined2, reduction = "umap",group.by = "orig.ident")
DimPlot(object.combined2, reduction = "umap",group.by = "seurat_clusters")
dev.off()

Idents(object.combined2) <- object.combined2$seurat_clusters

features <-c("RBFOX3", "GAD1", "NRGN","RIT2", "GALNTL6", "SLC17A6", "GAD2", "SLC6A3", "TH", "SLC18A2",
             "AQP4", "GFAP",
             "MOG",
             "C3", "CSF1R", "CD74", "TYROBP",
             "VCAN","DCC",
             "FLT1",
             "PDGFRB",
             "SKAP1",
             "COL1A2","ACTA2")
  
pdf("f1.pdf",width = 12,height = 15)

VlnPlot(object.combined2, features = features,pt.size = 0)
dev.off()

Idents(object.combined2) <- object.combined2$seurat_clusters
seurat_markers <- FindAllMarkers(object.combined2,logfc.threshold = 0.25,min.pct = 0.5,only.pos = T)

seurat_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) -> top5

# neurons (RBFOX3, GAD1, NRGN RIT2 GALNTL6 SLC17A6 GAD2 SLC6A3 TH SLC18A2),
# astrocytes (AQP4, GFAP),
# oligodendrocytes (MOG), 
# microglia (C3, CSF1R, CD74, TYROBP), 
# oligodendrocyte progenitor cells (VCAN DCC), 
# endothelial cells (FLT1),
# and pericytes (PDGFRB) 
#  T (SKAP1) 
# FIb (COL1A2)   ACTA2

object.combined2@meta.data$celltype <- as.character(object.combined2@meta.data$seurat_clusters)
object.combined2@meta.data$celltype[object.combined2@meta.data$celltype == "0"] = "oligodendrocytes"
object.combined2@meta.data$celltype[object.combined2@meta.data$celltype == "1"] = "oligodendrocytes"
object.combined2@meta.data$celltype[object.combined2@meta.data$celltype == "2"] = "microglia"
object.combined2@meta.data$celltype[object.combined2@meta.data$celltype == "3"] = "astrocytes"
object.combined2@meta.data$celltype[object.combined2@meta.data$celltype == "4"] = "endothelial"
object.combined2@meta.data$celltype[object.combined2@meta.data$celltype == "5"] = "oligodendrocyte_progenitor"
object.combined2@meta.data$celltype[object.combined2@meta.data$celltype == "6"] = "neurons"
object.combined2@meta.data$celltype[object.combined2@meta.data$celltype == "7"] = "neurons"
object.combined2@meta.data$celltype[object.combined2@meta.data$celltype == "8"] = "pericytes"
object.combined2@meta.data$celltype[object.combined2@meta.data$celltype == "9"] = "neurons"
object.combined2@meta.data$celltype[object.combined2@meta.data$celltype == "10"] = "oligodendrocytes"
object.combined2@meta.data$celltype[object.combined2@meta.data$celltype == "11"] = "FIb"
object.combined2@meta.data$celltype[object.combined2@meta.data$celltype == "12"] = "T"
object.combined2@meta.data$celltype[object.combined2@meta.data$celltype == "13"] = "neurons"
object.combined2@meta.data$celltype[object.combined2@meta.data$celltype == "14"] = "FIb"
object.combined2@meta.data$celltype[object.combined2@meta.data$celltype == "15"] = "oligodendrocytes"


DimPlot(object.combined2, reduction = "umap",label = T,group.by = "celltype")
DimPlot(object.combined2, reduction = "tsne",label = T,group.by = "celltype")

object.combined2_male <- subset(object.combined2, Gender == "Male")


Y_gene <- read.csv("chrY_gene2.csv")

Y_gene$gene_name -> Y_genes
Y_genes[Y_genes %in% rownames(object.combined2_male@assays$RNA@counts)] -> Y_genes

object.combined2_male@meta.data -> df
object.combined2_male@assays$RNA@counts -> counts
counts[Y_genes,,drop=F] -> counts_chrY
object.combined2_male@assays$RNA@data -> norm
norm[Y_genes,,drop=F] -> norm_chrY

as.data.frame(colMeans(as.matrix(norm_chrY))) -> chrY_norm; names(chrY_norm) <- "chrY_norm"
as.data.frame(colSums(as.matrix(counts_chrY))) -> chrY_UMI; names(chrY_UMI) <- "chrY_UMI"
as.data.frame(colSums(as.matrix(counts_chrY > 0))) ->  chrY_genes; names(chrY_genes) <- "chrY_genes"
as.data.frame(colSums(as.matrix(counts_chrY > 0))) ->  PAR_genes; names(PAR_genes) <- "PAR_genes"

df$chrY_genes <- chrY_genes$chrY_genes
df$chrY_UMI <- chrY_UMI$chrY_UMI
df$chrY_norm <- chrY_norm$chrY_norm

dplyr::mutate(df, LOY = ifelse(test = chrY_UMI == 0 , yes = "LOY", no = "NORMAL"))  -> df
dplyr::mutate(df, LOY_levels = ifelse(test = chrY_UMI > 2 & chrY_genes > 2 , yes = 4,
                                      no = ifelse(test = chrY_UMI > 1 & chrY_genes > 1, yes = 3,
                                                  no = ifelse(test = chrY_UMI > 1 | chrY_genes > 1, yes = 2,
                                                              no = ifelse(test = chrY_UMI >0, yes = 1, no = 0)))))  -> df
ifelse(df$LOY_levels == 0, yes = 0, no = ifelse(df$LOY_levels %in% c(1,2), yes = 1, no = 2) ) -> df$LOY_levels_condensed

dplyr::mutate(df, normY = log1p((chrY_UMI / nCount_RNA) * 10000)) -> df

object.combined2_male@meta.data <- df
object.combined2_male@assays$RNA@scale.data <- matrix()

df <- fread("/public2/data/wangjunhao/data/MGI_SingCell/public_singlecell/PD_brain/GSE184950_male.txt",data.table = F)

table(df$LOY)

sample <- unique(df[,c("orig.ident","Sample_title","Class","Age","Gender",
                       "Ethnicity","postmortem_interval_hours", "braak_stage"   )])
table(sample$Class)

## Extract top 5% of LOY cells and top 5% of non LOY cells
quantile(object.combined2_male$nFeature_RNA[object.combined2_male$LOY == "LOY"],0.95)
quantile(object.combined2_male$nFeature_RNA[object.combined2_male$LOY != "LOY"],0.95)

object.combined2_male2 <- subset(object.combined2_male, (LOY == "LOY"&nFeature_RNA>1262)|(LOY != "LOY"&nFeature_RNA>2230))
sample <- unique(object.combined2_male2@meta.data[,c("orig.ident","Sample_title","Class","Age","Gender",
                       "Ethnicity","postmortem_interval_hours", "braak_stage"   )])

### Perform transcription factor and pathway analysis on PD, PDD, and HC respectively

all.genes <- rownames(object.combined2_male2)
object.combined2_male2 <- ScaleData(object.combined2_male2, features = all.genes,
                             vars.to.regress = c("nFeature_RNA", "Age", "postmortem_interval_hours", "orig.ident"))


# base plot
library(Seurat)
library(tidyverse)
library(data.table)

Idents(object.combined2_male2) <- object.combined2_male2$celltype

DimPlot(object.combined2_male2,reduction="umap", label = F,pt.size=0.4,group.by = "celltype",
        label.color="#E64B35FF",label.size=5,
        cols = c("#CC8382","#A02120","#65A6AB","#C36D26",
                 "#B09CC3","#49277A","#4A8260","#6376B7",
                 "#804825")) 
DimPlot(object.combined2_male2,reduction="umap", label = F,pt.size=0.4,group.by = "LOY",
        label.color="#E64B35FF",label.size=5,
        cols = "Set1") 
DimPlot(object.combined2_male2,reduction="umap", label = F,pt.size=0.4,group.by = "Class",
        label.color="#E64B35FF",label.size=5,
        cols = c(
                 "#ee695b","#4da6a0","#f6b853"
                 )) 
DimPlot(object.combined2_male2,reduction="umap", label = F,pt.size=0.4,group.by = "braak_stage",
        label.color="#E64B35FF",label.size=5,
        cols = c("#A3256D","#5D5594", "#ee695b","#4da6a0","#f6b853")) 


fwrite(object.combined2_male2@meta.data,"/public2/data/wangjunhao/data/MGI_SingCell/public_singlecell/PD_brain/male_top5_metadata.txt",
       row.names = F,quote = F,sep = "\t")

meta <- object.combined2_male2@meta.data

cellLOY <- data.frame(table(meta$celltype,meta$LOY))
colnames(cellLOY) <- c("celltype","LOY","Freq")

num <- data.frame(table(meta$celltype))
colnames(num) <- c("celltype","Sum")
num <- num[order(num$Sum,decreasing = T),]
num$celltype <- factor(num$celltype,levels = c("oligodendrocytes","neurons","oligodendrocyte_progenitor", "microglia",
                                               "endothelial","astrocytes","pericytes","T",
                                               "FIb" ))

cellLOY <- left_join(cellLOY,num)
cellLOY$percent <- cellLOY$Freq/cellLOY$Sum

p <- ggplot(cellLOY, aes(x = celltype, y = percent, fill = LOY)) +
  ggplot2::geom_bar(stat="identity", position = "stack", alpha = 0.8) + 
  scale_fill_brewer(palette = "Set2")+      
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        strip.background = element_blank(),panel.border = element_blank(),
        axis.line=element_line())+
  # labs(title="housekeeping gene") +
  ylab("Ratio of LOY cells")
p2 <- ggplot(num, aes(x = "", y = Sum, fill = celltype)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Pie Chart", fill = "Category") +
  theme(plot.title = element_text(hjust = 0.5))

print(p2)
print(p)


## TF 

library(Seurat)
library(dorothea)
library(tidyverse)
library(data.table)

dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))
regulon <- dorothea_regulon_human %>%
  dplyr::filter(confidence %in% c("A","B","C"))

sce <- run_viper(object.combined2_male2, regulon,
                 options = list(method = "scale", minsize = 4,
                                eset.filter = FALSE, cores = 8,
                                verbose = FALSE))

DefaultAssay(object = sce) <- "dorothea"
Idents(sce) <- sce$celltype

library(future)

## We transform Viper scores, scaled by seurat, into a data frame to better 
## handling the results
viper_scores_df <- GetAssayData(sce, slot = "data", 
                                assay = "dorothea") %>%
  data.frame(check.names = F) %>%
  t()
# viper_scores_df[1:4,1:4]

## We create a data frame containing the cells and their clusters
CellsClusters <- data.frame(cell = rownames(object.combined2_male2@meta.data), 
                            cell_type = as.character(object.combined2_male2@meta.data$celltype),
                            LOY = as.character(object.combined2_male2@meta.data$LOY),
                            class = as.character(object.combined2_male2@meta.data$Class),
                            braak_stage = as.character(object.combined2_male2@meta.data$braak_stage),
                            check.names = F)

## We create a data frame with the Viper score per cell and its clusters
viper_scores_clusters <- viper_scores_df  %>%
  data.frame() %>% 
  rownames_to_column("cell") %>%
  gather(tf, activity, -cell) %>%
  inner_join(CellsClusters)

saveRDS(viper_scores_clusters,"maletop5_viper_scores.rds")





##  scheirerRayHare test   HC vs PD vs PDD  and stages

# if(!require(rcompanion)){install.packages("rcompanion")}
# if(!require(FSA)){install.packages("FSA")}
library(rcompanion)
library(FSA)

viper_scores_clusters <- readRDS("maletop5_viper_scores.rds")
viper_scores_clusters$type <- paste0(viper_scores_clusters$tf,"_",viper_scores_clusters$cell_type)
viper_scores_clusters$stage <- viper_scores_clusters$braak_stage
viper_scores_clusters$stage[viper_scores_clusters$stage == "N/A"] <- "Stage0"
viper_scores_clusters$stage[viper_scores_clusters$stage == "Stage II" | viper_scores_clusters$stage == "Stage III" ] <- "Stage1+2"
viper_scores_clusters$stage[viper_scores_clusters$stage == "Stage IV" | viper_scores_clusters$stage == "Stage VI"] <- "Stage3+4"
viper_scores_clusters$stage <- factor(viper_scores_clusters$stage,levels = c("Stage0","Stage1+2","Stage3+4"))

viper_scores_clusters$class[viper_scores_clusters$class=="Unaffected Control"] = "HC"
viper_scores_clusters$class[viper_scores_clusters$class=="Parkinson's Disease"] = "PD"
viper_scores_clusters$class[viper_scores_clusters$class=="Parkinson's Disease Dementia"] = "PDD"
viper_scores_clusters$class <- factor(viper_scores_clusters$class,levels = c("HC","PD","PDD"))

scheirer <- function(i,dat,condition){
  type <- unique(dat$type)
  data85 <- dat[dat$type == type[i],]
  data85$condition = data85[,which(colnames(data85) == condition)]
  p = scheirerRayHare(activity ~ LOY*condition,data = data85,verbose=F)
  p_aov = aov(activity ~ LOY * condition, data = data85)
  # summary(p_aov)[[1]]["LOY:condition","Pr(>F)"]
  
  value <- data.frame(p)
  value$type <- type[i]
  value$tf <- gsub("_.*","",type[i])
  value$condition <- condition
  value$effect <- rownames(value)
  value$aov <- summary(p_aov)[[1]][3,"Pr(>F)"]
  value
}

library(parallel)
a1 <- Sys.time()
Pvalue_class <- mclapply(1:length(unique(viper_scores_clusters$type)), 
                       function(i) try(scheirer(i,dat = viper_scores_clusters,condition="class"), TRUE),
                       mc.cores = 16)
a2 <- Sys.time() -a1
Pvalue_class <- data.frame(do.call(rbind,Pvalue_class),stringsAsFactors = F)
rownames(Pvalue_class) <- NULL
# sig_type <- Pvalue_LOY$type[Pvalue_LOY$p.value < 0.05 & Pvalue_LOY$effect=="LOY:condition"]
sig_type <- Pvalue_class[Pvalue_class$type %in% Pvalue_class$type[Pvalue_class$p.value < 0.05 & Pvalue_class$effect=="LOY:condition"] , ]

dim(Pvalue_class[Pvalue_class$p.value < 0.05 & Pvalue_class$effect=="LOY:condition",])




## GSEA
library(Seurat)
# library(SeuratData)
library(UCell)
library(irGSEA)
library(doMC)


object.combined2_male2 <- readRDS("/public2/data/wangjunhao/data/MGI_SingCell/public_singlecell/PD_brain/GSE184950_male_top5.rds")
Idents(object.combined2_male2) <- object.combined2_male2$celltype


object <- irGSEA.score(object = object.combined2_male2, assay = "RNA", 
                             slot = "data", seeds = 123, ncores = 10,
                             min.cells = 3, min.feature = 0,
                             custom = F, geneset = NULL, msigdb = T, 
                             species = "Homo sapiens", category = "H",  
                             subcategory = NULL, geneid = "symbol",
                             method = c("AUCell", "UCell", "singscore", 
                                        "ssgsea"),
                             aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                             kcdf = 'Gaussian')



object <- SeuratObject::SetIdent(object, value = "celltype")
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
                             group.by	= "LOY",
                             
                             ident.1 = "LOY",
                             ident.2 = "NORMAL",
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


# 1. Global show
irGSEA.heatmap.plot <- irGSEA.heatmap(object = result.dge, 
                                      method = "RRA",
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

scatterplot <- irGSEA.density.scatterplot(object = object,
                                          method = "UCell",
                                          show.geneset = "HALLMARK-INTERFERON-ALPHA-RESPONSE",
                                          
                                          reduction = "tsne")


pdf("brain5.pdf",width=12)
print(irGSEA.heatmap.plot)
print(irGSEA.bubble.plot)
print(irGSEA.upset.plot)
# print(upset.plot)
print(irGSEA.barplot.plot)
# print(scatterplot)
dev.off()

pdf("ttt5.pdf",width=12)
ComplexHeatmap::draw(ht)
dev.off()





result.dge_up <- result.dge$RRA[result.dge$RRA$direction == "up" & result.dge$RRA$pvalue < 0.05,]
result.dge_up$value <- -log10(result.dge_up$pvalue)
summary(result.dge_up$value)
result.dge_up <- result.dge_up[,c("Name","cluster","value")]
result.dge_up <- tidyr::spread(
  data=result.dge_up,
  key=cluster,
  value=value )
rownames(result.dge_up) <- result.dge_up$Name
result.dge_up <- result.dge_up[,-1]
result.dge_up[is.na(result.dge_up)] <- 0

result.dge_down <- result.dge$RRA[result.dge$RRA$direction == "down" & result.dge$RRA$pvalue < 0.05,]
result.dge_down$value <- -log10(result.dge_down$pvalue)
summary(result.dge_down$value)
result.dge_down <- result.dge_down[,c("Name","cluster","value")]
result.dge_down <- tidyr::spread(
  data=result.dge_down,
  key=cluster,
  value=value )
rownames(result.dge_down) <- result.dge_down$Name
result.dge_down <- result.dge_down[,-1]
result.dge_down[is.na(result.dge_down)] <- 0

library(ComplexHeatmap)
library(circlize)

col_fun <- colorRamp2(c(1, 3),   c("white", "#e64b35ff"))
up <- Heatmap(result.dge_up, name = "expression", col = col_fun)

col_fun <- colorRamp2(c(1, 6),   c("white", "#00468bff"))
down <- Heatmap(result.dge_down, name = "expression", col = col_fun)


pdf("brain_pathway.pdf",width = 12)
print(up)
print(down)
dev.off()
