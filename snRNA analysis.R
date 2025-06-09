
#############  AMP-PD TDEseq ##########

##   LOY analyse  

library(Matrix)
library(Rcpp)
library(reticulate)
use_python('/public/labdata/soft/Miniconda/envs/python3',required = T)

library(dplyr)
library(patchwork)
library(limma)
library(purrr)
library(data.table)
library(RColorBrewer)
library(tidyr)
library(ggplot2)
library(parallel)
library(Seurat)
library(openxlsx)
library(harmony)
library(fastmap)
library(tibble)
library(TDEseq)
library(ggrepel)
library(ggsci)
library(ggpubr)
library(magick)


setwd("/public2/data/AMP_PD/")

file <- list.files(pattern = "rds")

ProcessData <- function(i,file){
  brainData <- readRDS(file[i])
  
  ## for male
  brainData <- subset(brainData, SEX == "Male")
  
  brainData$age_at_baseline <- as.numeric(brainData$age_at_baseline)
  
  brainData <- Seurat::NormalizeData(brainData)

  Y_scLOY_genes <- readRDS("/public2/data/Y_scLOY_genes.RDS")

  PAR_scLOY_genes <- readRDS("/public2/data//PAR_scLOY_genes.RDS")
  
  PAR_scLOY_genes$gene_name -> PAR_genes
  unique(PAR_genes[PAR_genes %in% rownames(brainData)]) -> PAR_genes
  
  Y_scLOY_genes$gene_name -> Y_genes
  Y_genes[Y_genes %in% rownames(brainData)] -> Y_genes
  
  brainData@meta.data -> df
  brainData@assays$RNA@counts -> counts
  counts[Y_genes,,drop=F] -> counts_chrY
  counts[PAR_genes,,drop=F] -> counts_PAR
  brainData@assays$RNA@data -> norm
  norm[Y_genes,,drop=F] -> norm_chrY
  
  
  as.data.frame(colMeans(as.matrix(norm_chrY))) -> chrY_norm; names(chrY_norm) <- "chrY_norm"
  as.data.frame(colSums(as.matrix(counts_chrY))) -> chrY_UMI; names(chrY_UMI) <- "chrY_UMI"
  as.data.frame(colSums(as.matrix(counts_chrY > 0))) ->  chrY_genes; names(chrY_genes) <- "chrY_genes"
  as.data.frame(colSums(as.matrix(counts_PAR))) -> PAR_UMI;  names(PAR_UMI) <- "PAR_UMI"
  as.data.frame(colSums(as.matrix(counts_chrY > 0))) ->  PAR_genes; names(PAR_genes) <- "PAR_genes"
  
  df$chrY_genes <- chrY_genes$chrY_genes
  df$chrY_UMI <- chrY_UMI$chrY_UMI
  df$chrY_norm <- chrY_norm$chrY_norm
  df$PAR_UMI <- PAR_UMI$PAR_UMI
  df$PAR_genes <- PAR_genes$PAR_genes
  
  
  dplyr::mutate(df, LOY = ifelse(test = chrY_UMI == 0 , yes = "LOY", no = "NORMAL"))  -> df
  dplyr::mutate(df, LOY_levels = ifelse(test = chrY_UMI > 2 & chrY_genes > 2 , yes = 4,
                                        no = ifelse(test = chrY_UMI > 1 & chrY_genes > 1, yes = 3,
                                                    no = ifelse(test = chrY_UMI > 1 | chrY_genes > 1, yes = 2,
                                                                no = ifelse(test = chrY_UMI >0, yes = 1, no = 0)))))  -> df
  ifelse(df$LOY_levels == 0, yes = 0, no = ifelse(df$LOY_levels %in% c(1,2), yes = 1, no = 2) ) -> df$LOY_levels_condensed
  
  dplyr::mutate(df, normY = log1p((chrY_UMI / nCount_RNA) * 10000)) -> df
  dplyr::mutate(df, normPAR = log1p((PAR_UMI / nCount_RNA) * 10000)) -> df
   
  brainData@meta.data <- df

  name <- gsub(".rds","",file[i])
  pos <- paste0("/public2/data/AMP_PD/",name,"_male.rds")
  
  metadata1 <- brainData@meta.data
  pos1 <- paste0("/public2/data/AMP_PD/",name,"_metadata.txt")
  
  fwrite(metadata1,pos1,sep = "\t",quote = F,row.names = T)
  
  saveRDS(brainData,pos)
  pos
  
}
LOY_Process <- mclapply(1:length(file),  function(i) ProcessData(i=i,file=file), mc.cores = 3)
LOY_Process <- mclapply(12,  function(i) ProcessData(i=i,file=file), mc.cores = 3)



setwd("/public2/data/AMP_PD/")

files <- c("Immune_combinedData_male.rds")

d2 <- readRDS(files[1])
d2@meta.data <- d2@meta.data[,1:76]
d2_LOY <- subset(d2, LOY == "LOY")

LOY_data <- d2_LOY

table(LOY_data$subtype,LOY_data$region)


# filter Mic P2RY12
LOY_data_mir <- subset(LOY_data, subtype == "Mic P2RY12")

LOY_data_mir$path_braak_lb[is.na(LOY_data_mir$path_braak_lb)] <- 0


LOY_data_mir <- NormalizeData(LOY_data_mir)
LOY_data_mir$stage <- as.numeric(as.character(LOY_data_mir$path_braak_lb))
LOY_data_mir$stage[LOY_data_mir$stage == 0] <- "CTRL"
LOY_data_mir$stage[LOY_data_mir$stage %in% c(1,2)] <- "early"
LOY_data_mir$stage[LOY_data_mir$stage %in% c(3,4)] <- "medium"
LOY_data_mir$stage[LOY_data_mir$stage %in% c(5,6)] <- "last"

LOY_data_mir$region <- factor(LOY_data_mir$region,levels = c( "DMNX","GPI","PMC","PFC","PVC"))

LOY_data_mir$path_braak_lb <- factor(LOY_data_mir$path_braak_lb,levels = c(0,1,2,3,4,5,6))
LOY_data_mir$stage <- factor(LOY_data_mir$stage,levels = c("CTRL","early","medium","last"))


counts<-Seurat::GetAssayData(LOY_data_mir,slot ='counts')  ##raw counts data
data.norm<-Seurat::GetAssayData(LOY_data_mir,slot ='data') ##log normalized data
meta.data<-LOY_data_mir@meta.data                    ##metadata


tde <- CreateTDEseqObject(counts = counts,data=data.norm,meta.data=meta.data)

## tde for regions
tde_param5 <- list(sample.var = "batch",
                   stage.var = "region",
                   fit.model = "lmm",
                   tde.thr = 0.05,
                   pct = 0.1,
                   lfc = 0.1,
                   max.gcells = Inf,
                   min.tcells = 3)
tde5 <- tdeseq(object = tde, tde.method = "cell", tde.param = tde_param5, num.core = 3)

## tde for braak stages
tde_param6 <- list(sample.var = "batch",
                   stage.var = "stage",
                   fit.model = "lmm",
                   tde.thr = 0.05,
                   pct = 0.1,
                   lfc = 0.1,
                   max.gcells = Inf,
                   min.tcells = 3)
tde6 <- tdeseq(object = tde, tde.method = "cell", tde.param = tde_param6, num.core = 10)

objects_to_save <- c(  "tde5","tde6")

save(list = objects_to_save, file = paste0("/public2/data/AMP_PD/result/tde2.RData"))


####  markGenes plot
PatternHeatmap1 <- function (obj, stage.id, features = NULL, features.show = NULL, 
                             features.num = 50, cols = c("navy", "white", "firebrick3")){
  suppressPackageStartupMessages(library("Seurat"))
  suppressPackageStartupMessages(library("ggplot2"))
  suppressPackageStartupMessages(library("ComplexHeatmap"))
  suppressPackageStartupMessages(library("circlize"))
  data.norm <- GetTDEseqAssayData(obj, "data")
  mat <- Seurat::ScaleData(data.norm)
  metadata <- GetTDEseqAssayData(obj, "meta.data")
  
  res_dat = GetTDEseqAssayData(obj, "tde")
  res_dat = res_dat[order(res_dat$padj), ]
  feature_annot = c()
  if (is.null(features)) {
    features_plot = c()
    for (pattern in c("Growth", "Recession", "Peak", "Trough")) {
      idx = which(res_dat$pattern == pattern)
      if (length(idx) >= features.num) {
        features_length <- features.show[features.show %in% res_dat$gene[idx]]
        
        fe <- unique(c(res_dat$gene[idx[1:features.num]],features_length))
        
        features_plot = c(features_plot, fe)
        feature_annot = c(feature_annot, rep(pattern, 
                                             length(fe)))
        
      }      else {
        features_plot = c(features_plot, res_dat$gene[idx])
        feature_annot = c(feature_annot, rep(pattern, 
                                             length(idx)))
      }
    }
  }  else {
    features_plot = c()
    idx = match(features, res_dat$gene)
    if (any(is.na(idx))) {
      idx = idx[-which(is.na(idx))]
    }
    features = features[idx]
    subres_dat = res_dat[idx, ]
    for (pattern in c("Growth", "Recession", "Peak", "Trough")) {
      idx = which(subres_dat$pattern == pattern)
      features_plot = c(features_plot, res_dat$gene[idx])
      feature_annot = c(feature_annot, rep(pattern, length(idx)))
    }
  }
  group_info <- metadata$stage
  col_fun = colorRamp2(c(-2, 0, 2), cols)
  mat = mat[features_plot, ]
  if (!is.null(features.show)) {
    gene_pos <- match(features.show, rownames(mat))
    row_anno <- rowAnnotation(gene = anno_mark(at = gene_pos, 
                                               labels = features.show, labels_gp = gpar(fontface = 3)))
    f1 = Heatmap(mat, name = "Expression", col = col_fun, 
                 cluster_rows = FALSE, cluster_columns = F, show_column_dend = F, 
                 show_row_dend = F, show_row_names = FALSE, show_column_names = F, 
                 column_split = group_info, row_split = feature_annot, 
                 right_annotation = row_anno, border_gp = gpar(col = "black", 
                                                               lty = 2))
  }
  else {
    f1 = Heatmap(mat, name = "Expression", col = col_fun, 
                 cluster_rows = FALSE, cluster_columns = F, show_column_dend = F, 
                 show_row_dend = F, show_row_names = FALSE, show_column_names = F, 
                 column_split = group_info, row_split = feature_annot, 
                 border_gp = gpar(col = "black", lty = 2))
  }
  return(f1)
}


load("/public2/data/AMP_PD/result/tde2.RData")

t5 <- tde5@assays$RNA@tde
t6 <- tde6@assays$RNA@tde

Growth_genes <- sort(intersect(t5$gene[t5$pattern == "Growth"],t6$gene[t6$pattern == "Growth"]))
Recession_genes <- sort(intersect(t5$gene[t5$pattern == "Recession"],t6$gene[t6$pattern == "Recession"]))

G_genes <- c("ANP32A","ARHGEF7",
             "ATP13A3","BCL6",
             "FLT1","HIF1A",
             "PFKFB3","PTPN1")

R_genes <- c("ACTG1","C1QA",
             "EEF1A1","GPX1",
             "IL6ST","NPC2",
             "S100B","TREM2")


p31<-PatternHeatmap1(obj=tde6,stage.id='stage',features.show=c(G_genes,R_genes),
                   cols = c("#5Dbfe9", "white", "#c72228"))
print(p31)

tde5@assays[["RNA"]]@meta.data$stage = tde5@assays[["RNA"]]@meta.data$region
p41<-PatternHeatmap1(obj=tde5,stage.id='stage',features.show=c(G_genes,R_genes),
                   cols = c("#5Dbfe9", "white", "#c72228"))
print(p41)

pdf("/public2/data/AMP_PD/result/tde.pdf",height = 12,width = 12)
print(p31)
print(p41)
dev.off()


