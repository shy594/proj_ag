#write by Sun Haiayng at 2024.12.18
#01、搭环境
{
  rm(list=ls())
  for (i in c("ReactomePA","ggplot2","Seurat","SingleR","celldex","SingleCellExperiment",
              "SeuratObject","scRNAtoolVis","scDblFinder","presto","SeuratWrappers","ChIPseeker","stringr","DESeq2","clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db","TxDb.Mmusculus.UCSC.mm10.knownGene")) {library(i, character.only = T)};rm(i)
}
#02、给索引
{
  path = "/disk5/aaSHY/scAging/"
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/fa/GRCm38.dna.primary_assembly.fa"
  inx2 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/crREf_Run/chrNo"#inx4 = "/Reference/aaSHY/AppTools/crREf_Iso/chrNo"
  inx5 = "/Reference/aaSHY/AppTools/scTE/Data/shy.exclusive.idx"
  inxA = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MERVL-int::ERVL::LTR.bed"
  inxB = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MT2_Mm::ERVL::LTR.bed"
  inxC = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/IAPEY_LTR::ERVK::LTR.bed"
  inxD = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/RLTR1B-int::ERV1::LTR.bed"
  samp = "Ovary"
}
#03、写函数
{
  volcan <- function(f_Inputs, f_Output = "1111"){
    ress <- f_Inputs
    ress$cuts <- "theOthers"
    ress$cuts[ress$avg_log2FC >  log2(2.0) & ress$p_val_adj <0.05] = "Upup__Features"
    ress$cuts[ress$avg_log2FC < -log2(2.0) & ress$p_val_adj <0.05] = "Down__Features"
    lab_ <- ress[which(ress$cuts != "theOthers"),]
    lab_$rank <- abs(lab_$avg_log2FC)
    lab_ <- lab_[order(lab_$rank, decreasing = T),]
    lab_ <- lab_[1:20,]
    ress$cuts <- factor(ress$cuts, 
                        levels = c("Upup__Features","Down__Features","theOthers"))
    vaue <- ceiling(max(abs(ress$avg_log2FC)))
    vauf <- ceiling(max(-log10(ress$p_val_adj[which(ress$p_val_adj!=0)])))
    ress$p_val_adj <- ifelse(ress$p_val_adj < 10^(-300),10^(-300),ress$p_val_adj)
    p___01 <<- ggplot(data = ress) +
      geom_point(aes(x= avg_log2FC, y = -log10(p_val_adj), color = cuts), size=0.8) +
      labs(x = "Log2(FoldChange)", 
           y = "-Log10(adjusted p-value)", color = "") + 
      scale_x_continuous(limits = c(-vaue, vaue), breaks = seq(-vaue, vaue, 1)) +
      #scale_y_continuous(limits = c(0,vauf*1.2)) +
      geom_hline(aes(yintercept=-log10(0.05)), linetype="dashed", colour="royalblue") +
      geom_vline(aes(xintercept=-log2(2.0)), linetype="dashed", colour="royalblue") +
      geom_vline(aes(xintercept= log2(2.0)), linetype="dashed", colour="royalblue") +
      theme_bw() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      scale_color_manual(values = c("red","blue","black"),
                         labels = c(paste0(names(table(ress$cuts)[1])," (",unname(table(ress$cuts)[1]),")"),
                                    paste0(names(table(ress$cuts)[2])," (",unname(table(ress$cuts)[2]),")"),
                                    paste0(names(table(ress$cuts)[3])," (",unname(table(ress$cuts)[3]),")"))) +
      theme(plot.margin = margin(unit = "cm", c(1,1,1,1)),
            legend.position = c(0.2,0.9),
            legend.background = element_rect(fill = NA),
            legend.spacing.x = unit(0,"mm"),
            legend.spacing.y = unit(0,"mm"),
            legend.key = element_blank(),
            legend.text = element_text(size = 13),
            axis.title = element_text(size = 13)) +
      geom_text_repel(data = lab_,
                      aes(x= avg_log2FC, y = -log10(p_val_adj), label = rownames(lab_))) +
      guides(color = guide_legend(override.aes = list(size = 3)))
    #ggsave(plot = p___01, f_Output, units = "cm", width = 18, height = 16)
  }
}
################################################################################
#
#
#
#
#----------------------------------质控和预处理---------------------------------
#
#
#
#
################################################################################
#01、[merge----整合样本] [By scTE的结果]
{
  #options(Seurat.object.assay.version = "v5")
  prex <- c("Ocon","Ohucs","Ohucp")
  for (op in seq_along(prex)) {
    print(op)
    dirr <- paste0(path,"scTE_Iso/",prex[op],"/",prex[op],".csv")
    assign(paste0("tmp",op),as.data.frame(fread(dirr, sep=",", header = T)))
  }
  tmps <- list(tmp1,tmp2,tmp3)
  tmps <- lapply(tmps, function(x) {
    rownames(x) <- x[,1]
    x <- x[,-1]
    x <- t(x)
    x <- as(as.matrix(x), "dgCMatrix")
    x <- Seurat::CreateSeuratObject(counts=x, min.cells=3, min.features=200)
  })
  tmps[[1]]@meta.data$group <- prex[1]
  tmps[[2]]@meta.data$group <- prex[2]
  tmps[[3]]@meta.data$group <- prex[3]
  objs <- merge(tmps[[1]], y = c(tmps[[2]],tmps[[3]]),add.cell.ids = prex)
  #
  objs_1 <- objs
  rmGG <- c("Malat1", "Gm42418", "Gm26917", "AY036118")
  objs_1 = subset(objs_1, features = setdiff(rownames(objs_1),rmGG))
  rmMI <- rownames(objs_1)[grep("^mt",rownames(objs_1))]
  objs_1[["percMito"]] = PercentageFeatureSet(objs_1, features = rmMI)
  objs_1 = subset(objs_1, subset = percMito < 25)
  #
  options(future.globals.maxSize = 32 * 1e9)
  objs_2 <- objs_1
  objs_2 <- Seurat::NormalizeData(objs_2)
  objs_2 <- Seurat::FindVariableFeatures(objs_2, selection.method="vst", nfeatures=3000)
  objs_2 <- Seurat::ScaleData(objs_2)
  objs_2 <- JoinLayers(objs_2)
  #objs_2 <- Seurat::SCTransform(objs_2)
  #
  objs_3 <- objs_2
  objs_3 <- Seurat::RunPCA(objs_3)
  objs_3 <- Seurat::FindNeighbors(objs_3)
  objs_3 <- Seurat::FindClusters(objs_3, resolution = 0.5)
  data.used <- Seurat::Stdev(objs_3,reduction = "pca")
  vaue <- 0
  for (i in seq_along(data.used)) {
    vaue <- vaue + data.used[i]
    if(vaue > sum(data.used)*0.85) {
      print(i)
      vaue <- i
      break
    }
  }
  #
  objs_4 <- objs_3
  objs_4 <- Seurat::RunUMAP(objs_4, dims = 1:vaue)
  objs_4 <- Seurat::RunTSNE(objs_4, dims = 1:vaue)
}
#02
{
  #options(Seurat.object.assay.version = "v5")
  tmp1 <- Seurat::Read10X(paste0(path,"cellR_Iso/Ocon", "/outs/filtered_feature_bc_matrix"))
  tmp2 <- Seurat::Read10X(paste0(path,"cellR_Iso/Ohucp","/outs/filtered_feature_bc_matrix"))
  tmp3 <- Seurat::Read10X(paste0(path,"cellR_Iso/Ohucs","/outs/filtered_feature_bc_matrix"))
  abj1 <- Seurat::CreateSeuratObject(counts = tmp1); dim(abj1)
  abj2 <- Seurat::CreateSeuratObject(counts = tmp2); dim(abj2)
  abj3 <- Seurat::CreateSeuratObject(counts = tmp3); dim(abj3)
  # abj1 <- Seurat::CreateSeuratObject(counts = tmp1, min.cells = 3, min.features = 200); dim(abj1)#
  # abj2 <- Seurat::CreateSeuratObject(counts = tmp2, min.cells = 3, min.features = 200); dim(abj2)#
  # abj3 <- Seurat::CreateSeuratObject(counts = tmp3, min.cells = 3, min.features = 200); dim(abj3)#
  abj1@meta.data$group <- "Ocon"
  abj2@meta.data$group <- "Ohucp"
  abj3@meta.data$group <- "Ohucs"
  objs <- merge(abj1,y = list(abj2,abj3), add.cell.ids = c("Ocon","Ohucp","Ohucs"))
  #
  objs_1 <- objs
  rmGG <- c("Malat1", "Gm42418", "Gm26917", "AY036118") #Genes representing ribosomal contamination
  objs_1 = subset(objs_1, features = setdiff(rownames(objs_1),rmGG))
  rmMI <- rownames(objs_1)[grep("^mt",rownames(objs_1))]
  objs_1[["percMito"]] = PercentageFeatureSet(objs_1, features = rmMI)
  objs_1 = subset(objs_1, subset = percMito < 25)
  #
  options(future.globals.maxSize = 32 * 1e9)
  objs_2 <- objs_1
  objs_2 <- Seurat::NormalizeData(objs_2)
  objs_2 <- Seurat::FindVariableFeatures(objs_2, selection.method="vst", nfeatures=2000)
  objs_2 <- Seurat::ScaleData(objs_2)
  objs_2 <- JoinLayers(objs_2)
  #objs_2 <- Seurat::SCTransform(objs_2)
  #
  objs_3 <- objs_2
  objs_3 <- Seurat::RunPCA(objs_3)
  objs_3 <- Seurat::FindNeighbors(objs_3)
  objs_3 <- Seurat::FindClusters(objs_3, resolution = 0.5)
  data.used <- Seurat::Stdev(objs_3,reduction = "pca")
  vaue <- 0
  for (i in seq_along(data.used)) {
    vaue <- vaue + data.used[i]
    if(vaue > sum(data.used)*0.85) {
      print(i)
      vaue <- i
      break
    }
  }
  #
  objs_4 <- objs_3
  objs_4 <- Seurat::RunUMAP(objs_4, dims = 1:vaue)
  ###objs_4 <- Seurat::RunTSNE(objs_4, dims = 1:10)
  ###objs_4@meta.data$group <- factor(objs_4@meta.data$group, levels = c("Ocon","Ohucs","Ohucp"))
  ###objs_4 <- subset(objs_4, cells = WhichCells(objs_4, expression = group == "Ohucp"))
  p___01 <- Seurat::DimPlot(objs_4, reduction = "umap",
                            group.by = c("seurat_clusters"),alpha = 0.2,
                            cols = ggsci::pal_d3(palette = "category20")(20),
                            #label = T, label.size = 5, label.color = "black",
                            pt.size = .5)
  p___01
  objs_4@meta.data$group <- factor(objs_4@meta.data$group, levels = c("Ocon","Ohucs","Ohucp"))
  p___02 <- Seurat::DimPlot(objs_4, reduction = "umap",
                            group.by = c("group"),alpha = 0.2,
                            cols = ggsci::pal_d3(palette = "category20")(3),
                            #label = T, label.size = 5, label.color = "black",
                            pt.size = .5)
  p___02
  #
  #
  #
  objs_5 <- subset(objs_4, cells = WhichCells(objs_4, expression = group != "Ohucp"))
  p___03 <- Seurat::DimPlot(objs_5, reduction = "umap",
                            group.by = c("group"),alpha = 0.3,
                            cols = c("steelblue","indianred"),
                            #label = T, label.size = 5, label.color = "black",
                            pt.size = .5)
  p___03
  objs_6 <- subset(objs_4, cells = WhichCells(objs_4, expression = group != "Ohucs"))
  p___04 <- Seurat::DimPlot(objs_6, reduction = "umap",
                            group.by = c("group"),alpha = 0.3,
                            cols = c("steelblue","indianred"),
                            #label = T, label.size = 5, label.color = "black",
                            pt.size = .5)
  p___04
  objs_7 <- subset(objs_4, cells = WhichCells(objs_4, expression = group == "Ocon"))
  p___05 <- Seurat::DimPlot(objs_7, reduction = "umap",
                            group.by = c("group"),alpha = 0.3,
                            cols = c("steelblue"),
                            #label = T, label.size = 5, label.color = "black",
                            pt.size = .5)
  p___05
}
################################################################################
#
#
#
#
#----------------------------------分析与可视化---------------------------------
#
#
#
#
################################################################################
#01、分群
{
  ###objs_4@meta.data$group <- factor(objs_4@meta.data$group, levels = prex)
  ###objs_4 <- subset(objs_4, cells = WhichCells(objs_4, expression = group == prex[3]))
  p___01 <- Seurat::DimPlot(objs_4, reduction = "umap",
                            group.by = c("seurat_clusters"),alpha = 0.2,
                            cols = c(ggsci::pal_d3(palette = "category20")(20),"grey"),
                            label = T, label.size = 5, label.color = "black",
                            pt.size = .5)
  p___01
  # test <- subset(objs_4, cells = WhichCells(objs_4, expression = seurat_clusters %in% c(5,16,20)))
  # Seurat::DimPlot(test, reduction = "umap",
  #                 group.by = c("seurat_clusters"),alpha = 0.2,
  #                 cols = c(ggsci::pal_d3(palette = "category20")(20),"grey"),
  #                 label = T, label.size = 5, label.color = "black",
  #                 pt.size = .5)
}
################################################################################
#02、细胞占比
{
  #每个cluster中, 分属不同样本的细胞的数量占比
  p___02 <- ggplot(data = objs_4@meta.data) + 
    geom_bar(aes(x = seurat_clusters, fill = group), position = "fill") +
    scale_fill_jama()
  p___02
  #
  #
  #
  #每个cluster中, 分属不同样本的细胞的数量, 在相应样本的总细胞里的占比
  red1 <- objs_4@meta.data[,c("seurat_clusters","group")]
  red1 <- red1 %>% 
    dplyr::group_by(seurat_clusters,group) %>%
    dplyr::summarise(cou1 = n(), .groups = "drop") %>%
    dplyr::mutate(cou2 = sum(cou1), cou3 = cou1/cou2, .by = group)
  options(scipen=200)
  red1$cou3 <- round(red1$cou3,3)
  p___03 <- ggplot(data = red1) +
    geom_bar(aes(x = seurat_clusters, y = cou3), stat = "identity") +
    facet_wrap(vars(group),nrow = 3) +
    scale_y_continuous(limits = c(0,0.5)) +
    geom_text(aes(x = seurat_clusters, y = cou3+0.06, label = cou3),
              size=4,color="indianred",angle=75,vjust =.5,hjust=.5)
  p___03
  ggsave(plot = p___03,units = "cm",width = 20, height = 15,
         filename = paste0(path,"PLOTs/Iso/",samp,"_cellRatio-v1.pdf"))
}
# 跑完发现了某些类群的分属于不同组的细胞的比例发生变化
{
  test <- subset(objs_4, cells = WhichCells(objs_4, expression = seurat_clusters %in% c(0,1,4,6,8)))
  Seurat::DimPlot(test, reduction = "umap",
                  group.by = c("seurat_clusters"),alpha = 0.2,
                  cols = c(ggsci::pal_d3(palette = "category20")(20),"grey"),
                  label = T, label.size = 5, label.color = "black",
                  pt.size = .5)
}
################################################################################
#03、各个样本细胞的分布
{
  objs_4@meta.data$group <- factor(objs_4@meta.data$group, levels = prex)
  p___04 <- Seurat::DimPlot(objs_4, reduction = "umap",
                            group.by = c("group"),alpha = 0.2,
                            cols = ggsci::pal_d3(palette = "category20")(3),
                            #label = T, label.size = 5, label.color = "black",
                            pt.size = .5)
  p___04
  #
  #
  #
  objs_5 <- subset(objs_4, cells = WhichCells(objs_4, expression = group == prex[1]))
  p___05 <- Seurat::DimPlot(objs_5, reduction = "umap",
                            group.by = c("group"),alpha = 0.3,
                            cols = c("steelblue"),
                            #label = T, label.size = 5, label.color = "black",
                            pt.size = .5)
  p___05
  #
  #
  #
  objs_6 <- subset(objs_4, cells = WhichCells(objs_4, expression = group != prex[3]))
  p___06 <- Seurat::DimPlot(objs_6, reduction = "umap",
                            group.by = c("group"),alpha = 0.3,
                            cols = c("steelblue","indianred"),
                            #label = T, label.size = 5, label.color = "black",
                            pt.size = .5)
  p___06
  #
  #
  #
  objs_7 <- subset(objs_4, cells = WhichCells(objs_4, expression = group != prex[2]))
  p___07 <- Seurat::DimPlot(objs_7, reduction = "umap",
                            group.by = c("group"),alpha = 0.3,
                            cols = c("steelblue","indianred"),
                            #label = T, label.size = 5, label.color = "black",
                            pt.size = .5)
  p___07
}
################################################################################
#04、可视化feature的表达分布
{
  aims <- c("MT2-Mm","MERVL-int")
  RidgePlot(objs_4, features = aims, ncol = 2, group.by ="group")
  VlnPlot(objs_4, features = aims, ncol = 2, group.by ="group", pt.size = 0)
  FeaturePlot(objs_4, features = aims)
  DotPlot(objs_4, features = aims, group.by = "group") + RotatedAxis()
  DotPlot(objs_4, features = aims, group.by = "seurat_clusters") + RotatedAxis()
  DoHeatmap(objs_4, features = VariableFeatures(objs_4)[1:10], size = 3)
}
{
  aims <- c("TP53","CDKN2A","CDKN1A","IFNB1","IFNG","IFNAR1","EGFR")
  aims <- babelgene::orthologs(genes = aims, species = "mouse", human = T)$symbol
  #
  aims <- c("Brca1", "Brca2", "Myc", "Kras", "Apc", "Akt2","Eif5a2", "Rab25", "Pik3ca", "Egfr")
  FeaturePlot(objs_4, features = aims)
}
################################################################################
#05、细胞类型注释 [手动注释]
{
  objs_8 <- objs_4
  aims <- c("Col1a1","Bgn","Ogn","Dcn","Lum","Notch3","Cd68","Foxl2",
            "Inha","Cyp19a1","Star","Cyp11a1",
            "Epcam","Krt8","Muc1",#"Amh","Gpm6a","Ptgfr","Fshr","Lhcgr","Pgr","Hsd3b",
            "Srd5a1","C1qa","Cd34","Cd3g","Upk1b","Zp3","Cd79a")
  for (mm in aims) {
    print(mm)
    p___01 <- VlnPlot(objs_8, features = mm, ncol = 1, group.by ="seurat_clusters", pt.size = 0)
    ggsave(plot = p___01,
           height = 16, width = 16, units = "cm",
           file = paste0(path,"PLOTs/Iso/OvaryMarker/",mm,".Vln.pdf"))
    p___02 <- FeaturePlot(objs_8, features = mm)
    ggsave(plot = p___02,
           height = 16, width = 16, units = "cm",
           file = paste0(path,"PLOTs/Iso/OvaryMarker/",mm,".Dis.pdf"))
  }
  objs_8$cellType <- "Unknown"
  objs_8$cellType[which(objs_8$seurat_clusters %in% c(0,15))] <- "Stroma A"
  objs_8$cellType[which(objs_8$seurat_clusters %in% c(3,17))] <- "Stroma B"
  objs_8$cellType[which(objs_8$seurat_clusters %in% c(8))] <- "Stroma C"
  objs_8$cellType[which(objs_8$seurat_clusters %in% c(6,20))] <- "GCs"
  objs_8$cellType[which(objs_8$seurat_clusters %in% c(4,10))] <- "TCs"
  objs_8$cellType[which(objs_8$seurat_clusters %in% c(5,14))] <- "Phagocytes"
  objs_8$cellType[which(objs_8$seurat_clusters %in% c(2,18,19))] <- "Endothelial cells"
  objs_8$cellType[which(objs_8$seurat_clusters %in% c(11,17))] <- "T lymphocytes"
  objs_8$cellType[which(objs_8$seurat_clusters %in% c(7))] <- "Epithelial cells"
  objs_8$cellType[which(objs_8$seurat_clusters %in% c(9,12))] <- "Luteal cells"
  p___03 <- Seurat::DimPlot(objs_8, reduction = "umap",
                            group.by = c("cellType"),alpha = 0.2,
                            cols = c(ggsci::pal_d3(palette = "category20")(10),"grey60"),
                            label = T, label.size = 5, label.color = "black",
                            pt.size = .5)
  p___03
}
################################################################################
#06、细胞类型注释 [自动注释]
{
  refc <- celldex::MouseRNAseqData()
  objs_A <- GetAssayData(objs_8,layer="data")
  objs_B <- SingleR(test = objs_A, ref = refc, labels = refc$label.main)
  objs_8$typeCell <- objs_B$pruned.labels
  objs_8 <- objs_8[,which(!(is.na(objs_8$typeCell)))]
  p___03 <- Seurat::DimPlot(objs_8, reduction = "umap",
                            group.by = c("typeCell"),alpha = 0.2,
                            cols = c(ggsci::pal_d3(palette = "category20b")(17)),
                            label = T, label.size = 4, label.color = "black",
                            pt.size = .2)
  p___03
}
################################################################################
#07、差异分析---- [FindMarkers] [簇间]
{
  #目标：[5,19 __vs__ 1,12,13,14]
  #预处理
  objs_9 <- objs_4
  conditon <- c(5,19); controls <- c(1,12,13,14)
  clu__1 <- Idents(objs_9)
  clu__2 <- ifelse(clu__1 %in% conditon, "A", ifelse(clu__1 %in% controls,"B", "other"))
  Idents(objs_9) <- clu__2
  #
  #
  #基因的差异
  red1 <- read.table(file = "/Reference/aaSHY/BED/special/allGene-2.bed")
  red1 <- unique(red1$V4)
  aims__g <- intersect(rownames(objs_9), red1)
  objs_a <- subset(x = objs_9, features = aims__g)
  mak__a <- Seurat::FindMarkers(object = objs_a, min.pct = 0.20,
                                ident.1 = "A", ident.2 = "B", assay = "RNA", slot="data")
  volcan(f_Inputs = mak__a); p___01
  #
  #
  #重复序列的差异
  red2 <- read.csv(file = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/TE.Relation.csv")
  red2 <- gsub("_","-",red2$gene_id)
  aims__t <- intersect(rownames(objs_9), red2)
  objs_b <- subset(x = objs_9, features = aims__t)
  mak__b <- Seurat::FindMarkers(object = objs_b, min.pct = 0.20,
                                ident.1 = "A", ident.2 = "B", assay = "RNA", slot="data")
  volcan(f_Inputs = mak__b); p___01
  #
  #
  #
  #
  #
  #
  test <- mak__b[which(mak__b$p_val_adj <0.05 & mak__b$avg_log2FC > 0.58),]
  test <- test[order(test$avg_log2FC, decreasing = T),]
  rownames(test[1:41,])
  rownames(test[1:20,])
}
################################################################################
#08、差异分析---- [FindMarkers] [簇内]
{
  #目标：[在cluster中, Ohucs __vs__ Ocon]
  #目标：[在cluster中, Ohucp __vs__ Ocon]
  #预处理
  xx = 3; print(prex[xx]) #选择哪个处理组
  yy = 1 #选择哪个cluster
  objs_9 <- subset(objs_4, cells = WhichCells(objs_4, expression = group %in% c(prex[1],prex[xx])))
  objs_9 <- subset(objs_9, subset= seurat_clusters == yy)
  Idents(objs_9) <- objs_9$group
  #
  #
  #基因的差异
  red1 <- read.table(file = "/Reference/aaSHY/BED/special/allGene-2.bed")
  red1 <- unique(red1$V4)
  aims__g <- intersect(rownames(objs_9), red1)
  objs_a <- subset(x = objs_9, features = aims__g)
  mak__a <- Seurat::FindMarkers(object = objs_a, min.pct = 0.20,
                                assay = "RNA", slot="data",
                                ident.1 = prex[xx], ident.2 = prex[1])
  volcan(f_Inputs = mak__a); p___01
  ggsave(plot = p___01,
         height = 16, width = 16, units = "cm",
         file = paste0(path,"PLOTs/Iso/",prex[xx],"_vs_",prex[1],"-clu",yy,"-gene.pdf"))
  #
  #
  #重复序列的差异
  red2 <- read.csv(file = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/TE.Relation.csv")
  red2 <- gsub("_","-",red2$gene_id)
  aims__t <- intersect(rownames(objs_9), red2)
  objs_b <- subset(x = objs_9, features = aims__t)
  mak__b <- Seurat::FindMarkers(object = objs_b, min.pct = 0.20,
                                assay = "RNA", slot="data",
                                ident.1 = prex[xx], ident.2 = prex[1])
  volcan(f_Inputs = mak__b); p___01
  ggsave(plot = p___01,
         height = 16, width = 16, units = "cm",
         file = paste0(path,"PLOTs/Iso/",prex[xx],"_vs_",prex[1],"-clu",yy,"-TE.pdf"))
  #
  #
  #
  #
  #
  #
  test <- mak__b[which(mak__b$p_val_adj <0.05 & mak__b$avg_log2FC > 0.58),]
  test <- test[order(test$avg_log2FC, decreasing = T),]
  rownames(test[1:41,])
  rownames(test[1:20,])
}
################################################################################
#09、差异分析---- [FindAllMarkers]
{
  future::plan("multisession", workers=2)
  objs_C <- objs_8
  Idents(objs_C) <- objs_C$cellType
  #
  #
  #基因的差异
  red1 <- read.table(file = "/Reference/aaSHY/BED/special/allGene-2.bed")
  red1 <- unique(red1$V4)
  aims__g <- intersect(rownames(objs_C), red1)
  objs_C1 <- subset(x = objs_C, features = aims__g)
  mak__1 <- Seurat::FindAllMarkers(object = objs_C1, min.pct = 0.2, logfc.threshold = 0.58)
  scRNAtoolVis::jjVolcano(diffData = mak__1,
                          adjustP.cutoff = 0.05, 
                          tile.col = jjAnno::useMyCol("paired",n = 11))
  scRNAtoolVis::markerVolcano(markers = mak__1,topn = 2,
                              log2FC = 0.58,
                              labelCol = jjAnno::useMyCol("paired",n = 11))
  tmp_gg <- mak__1 %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(order_by = avg_log2FC, n = 2)
  scRNAtoolVis::jjDotPlot(object = objs_C1,
                          gene = tmp_gg$gene,
                          dot.col = ggsci::pal_npg(palette = c("nrc"))(10)[2:4],
                          plot.margin = c(1,0,0,1),
                          id = "celltype",split.by = "group",
                          xtree = F)
  #
  #
  #重复序列的差异
  red2 <- read.csv(file = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/TE.Relation.csv")
  red2 <- gsub("_","-",red2$gene_id)
  aims__t <- intersect(rownames(objs_C), red2)
  objs_C2 <- subset(x = objs_C, features = aims__t)
  mak__2 <- Seurat::FindAllMarkers(object = objs_C2, min.pct = 0.2, logfc.threshold = 0.58)
  scRNAtoolVis::jjVolcano(diffData = mak__2,
                          adjustP.cutoff = 0.05, 
                          tile.col = jjAnno::useMyCol("paired",n = 12),
                          polar = F, flip = F)
  scRNAtoolVis::markerVolcano(markers = mak__2,topn = 2,
                              log2FC = 0.58,
                              labelCol = jjAnno::useMyCol("paired",n = 11))
  tmp_tt <- mak__2 %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(order_by = avg_log2FC, n = 2)
  scRNAtoolVis::jjDotPlot(object = objs_C2,
                          gene = tmp_tt$gene,
                          dot.col = ggsci::pal_npg(palette = c("nrc"))(10)[2:4],
                          plot.margin = c(1,0,0,1),
                          id = "celltype",split.by = "group",
                          xtree = F)
}
#10、高变feature的热图
{
  future::plan("multisession", workers=2)
  objs_H <- objs_8
  Idents(objs_H) <- objs_8$cellType
  mak__3 <- Seurat::FindAllMarkers(object = objs_H, min.pct = 0.2, logfc.threshold = 0.58)
  red2 <- read.csv(file = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/TE.Relation.csv")
  red2 <- gsub("_","-",red2$gene_id)
  mak__3 <- mak__3[which(rownames(mak__3) %in% red2),]
  tttttt <- mak__3[which(mak__3$p_val_adj < 0.05),]
  #DoHeatmap(objs_H,slot = "counts",
  #          features = as.character(unique(tttttt$gene)),
  #          group.by = "cellType", assay = "RNA",size = 3)
  dittoHeatmap(object = objs_H, 
               scale = "row",treeheight_row = 10,
               scaled.to.max = F, show_rownames=T,
               breaks = seq(-2,2,by=0.01),
               heatmap.colors = colorRampPalette(c('#1A5592','white',"#B83D3D"))(400),
               genes = as.character(unique(tttttt$gene)),
               annot.by = c("cellType"))#,highlight.features = c("RLTR13A")
  
}
{
  ###future::plan("multisession", workers=2)
  ###objs_D <- objs_8
  ###red1 <- read.table(file = "/Reference/aaSHY/BED/special/allGene-2.bed")
  ###red1 <- unique(red1$V4)
  ###aims__g <- intersect(rownames(objs_D), red1)
  ###objs_D1 <- subset(x = objs_D, features = aims__g)
  ###mak__3 <- Seurat::FindAllMarkers(object = objs_D1, min.pct = 0.2, logfc.threshold = 0.58)
  ###gggggg <- mak__3[which(abs(mak__3$avg_log2FC) >5 & mak__3$p_val_adj < 0.01),]
  ###DoHeatmap(objs_D1,
  ###          features = as.character(unique(gggggg$gene)),
  ###          group.by = "seurat_clusters", assay = "RNA")
  ###DoHeatmap(objs_D1, 
  ###          features = VariableFeatures(objs_D1))
  #
  #
  #
  future::plan("multisession", workers=2)
  objs_D <- objs_8
  Idents(objs_D) <- objs_8$cellType
  red2 <- read.csv(file = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/TE.Relation.csv")
  red2 <- gsub("_","-",red2$gene_id)
  aims__t <- intersect(rownames(objs_D), red2)
  objs_D2 <- subset(x = objs_D, features = aims__t)
  mak__4 <- Seurat::FindAllMarkers(object = objs_D2, min.pct = 0.01, logfc.threshold = 0.58)
  tttttt <- mak__4#[which(mak__4$p_val_adj < 0.05),]
  tttttt[rownames(tttttt)[2],]
  View(tttttt)
  
  DoHeatmap(objs_D2,size = 3,
            features = unique(tttttt$gene)[1:886],
            disp.min = -4, disp.max = 4,
            group.by = "cellType", assay = "RNA") +
    scale_fill_gradientn(colors = c("navy", "white", "firebrick3"))
  
  intersect(Features(objs_D2),as.character(unique(tttttt$gene)))
  
  
  tttttt <- mak__4[which(abs(mak__4$avg_log2FC) >2 & mak__4$p_val_adj < 0.05),]
  dittoHeatmap(object = objs_D2, 
               scale = "row", treeheight_row = 10,
               scaled.to.max = F, show_rownames=F,
               breaks = seq(-2,2,by=0.01),
               heatmap.colors = colorRampPalette(c('navy','white',"firebrick3"))(400),
               genes = as.character(unique(tttttt$gene)),
               annot.by = c("cellType"))#,highlight.features = c("RLTR13A")
}




