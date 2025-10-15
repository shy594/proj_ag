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
#01、UMAP聚类 (基于Cell Ranger的结果)
{
  #options(Seurat.object.assay.version = "v5")
  tmp1 <- Seurat::Read10X(paste0(path,"cellR_Iso/Ocon", "/outs/filtered_feature_bc_matrix"))
  abj1 <- Seurat::CreateSeuratObject(counts = tmp1); dim(abj1)
  #tmp1 <- Seurat::Read10X(paste0(path,"cellR_Iso/Ocon", "/outs/raw_feature_bc_matrix"))
  #abj1 <- Seurat::CreateSeuratObject(counts = tmp1, 
  #                                   min.cells = 3, min.features = 300); dim(abj1)#4，4
  
  abj1@meta.data$group <- "Ocon"
  objs_1 <- abj1
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
  p___01 <- Seurat::DimPlot(objs_4, reduction = "umap",
                            group.by = c("seurat_clusters"), alpha = 0.5,
                            cols = ggsci::pal_d3(palette = "category20")(20),
                            #label = T, label.size = 5, label.color = "black",
                            pt.size = .5)
  p___01
  p___02 <- Seurat::DimPlot(objs_4, reduction = "tsne",
                            group.by = c("seurat_clusters"),alpha = 0.2,
                            cols = ggsci::pal_d3(palette = "category20")(20),
                            #label = T, label.size = 5, label.color = "black",
                            pt.size = .5)
  p___02
}
################################################################################
#
#
#
#
#
#
#----------------------------------分析与可视化---------------------------------
#
#
#
#
#
#
################################################################################
#01、查看具体CLUs的细胞分布
{
  test <- subset(objs_4, cells = WhichCells(objs_4, expression = seurat_clusters %in% c(0,1)))
  p___03 <- Seurat::DimPlot(test, reduction = "umap",
                            group.by = c("seurat_clusters"),alpha = 0.2,
                            cols = ggsci::pal_d3(palette = "category20")(20),
                            label = T, label.size = 5, label.color = "black",
                            pt.size = .5)
  p___03
}
################################################################################
#04、可视化feature的表达分布
{
  aims <- c("Bgn","Dcn")
  RidgePlot(objs_4, features = aims, ncol = 2, group.by ="seurat_clusters")
  VlnPlot(objs_4, features = aims, ncol = 2, group.by ="seurat_clusters", pt.size = 0)
  FeaturePlot(objs_4, features = aims)
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
  
  # 补充一些marker
  for (mm in aims) {
    print(mm)
    p___04 <- VlnPlot(objs_8, features = mm, ncol = 1, group.by ="seurat_clusters", pt.size = 0)
    ggsave(plot = p___04,
           height = 16, width = 16, units = "cm",
           file = paste0(path,"PLOTs/Iso/OvaryMarker/",mm,".Vln.pdf"))
    
    p___05 <- FeaturePlot(objs_8, features = mm); p___05
    ggsave(plot = p___05,
           height = 16, width = 16, units = "cm",
           file = paste0(path,"PLOTs/Iso/OvaryMarker/",mm,".Dis.pdf"))
  }
  objs_8 <- objs_4
  objs_8$cellType <- "Unknown"
  objs_8$cellType[which(objs_8$seurat_clusters %in% c(11))] <- "Epithelial"#Epcam, Krt8
  objs_8$cellType[which(objs_8$seurat_clusters %in% c(3,7))] <- "Luteal"#Cyp11a1
  objs_8$cellType[which(objs_8$seurat_clusters %in% c(0,1,5))] <- "Stroma"#Dcn, Ogn
  objs_8$cellType[which(objs_8$seurat_clusters %in% c(4,9))] <- "Granulosa"#Inha
  objs_8$cellType[which(objs_8$seurat_clusters %in% c(6,10))] <- "Lymphocytes"#Ptprc
  objs_8$cellType[which(objs_8$seurat_clusters %in% c(2))] <- "Endothelial" #Cd34
  objs_8$cellType[which(objs_8$seurat_clusters %in% c(8))] <- "Vascular smooth muscle" #Myh11,Pde3a
  # 缺少巨噬细胞、鞘细胞、T淋巴细胞
  p___04 <- Seurat::DimPlot(objs_8, reduction = "umap",
                            group.by = c("cellType"),alpha = 0.2,
                            #cols = ggsci::pal_npg()(10),
                            cols = ggthemes::few_pal()(8),
                            label = F, label.size = 4, label.color = "black",
                            pt.size = .5)
  p___04
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
  mak__3 <- Seurat::FindAllMarkers(object = objs_H, 
                                   min.pct = 0.2, logfc.threshold = 0.58)
  mak__3 <- mak__3[which(mak__3$p_val_adj < 0.05),]
  # mak__4 <- mak__3[which(mak__3$cluster %in% c(8)),]
  
  
  ervs_01 <- read.table("/ChIP_seq_1/aaSHY/pacbio/assemIso4/dumpROOM/isoform_gene_ervs-info.txt", 
                        header = T)
  ervs_02 <- unique(ervs_01$gene_name)
  
  DoHeatmap(objs_H,
            slot = "counts",
            features = as.character(ervs_02),
            group.by = "cellType", assay = "RNA", size = 3)
  dittoHeatmap(object = objs_H, 
               scale = "row", treeheight_row = 10,
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
################################################################################
{
  future::plan("multisession", workers=2)
  objs_H <- objs_8
  mak__3 <- Seurat::FindAllMarkers(object = objs_H, 
                                   min.pct = 0.2, logfc.threshold = 0.58)
  mak__3 <- mak__3[which(mak__3$p_val_adj < 0.05),]
  # mak__4 <- mak__3[which(mak__3$cluster %in% c(8)),]
  
  
  ervs_01 <- read.table("/ChIP_seq_1/aaSHY/pacbio/assemIso4/dumpROOM/isoform_gene_ervs-info.txt", 
                        header = T)
  ervs_01 <- ervs_01[which(ervs_01$subfamily_id %in% c("MERVL-int","MT2_Mm")),]
  ervs_02 <- unique(ervs_01$gene_name)
  # 1. 筛选在数据中存在的基因（排除未检测到的基因）
  valid_genes <- intersect(ervs_02, rownames(objs_8))
  cat("有效基因数量：", length(valid_genes), "\n")  # 确保保留了足够的基因
  
  # 2. 提取这些基因的表达矩阵（默认使用标准化后的data slot）
  expr_matrix <- GetAssayData(objs_8, slot = "data")[valid_genes, ]
  
  # 3. 计算每个细胞中这些基因的平均表达量（作为“虚拟基因”）
  # 按列（细胞）取均值，结果是一个细胞×1的向量
  combined_expr <- colMeans(expr_matrix)
  
  # 4. 将这个“虚拟基因”添加到Seurat对象的meta.data中
  objs_8$Combined_Genes <- combined_expr
  
  # 5. 使用FeaturePlot可视化这个“虚拟基因”的表达分布
  FeaturePlot(
    object = objs_8,
    features = "Combined_Genes",  # 使用meta.data中的列名
    reduction = "umap",  # 或"tsne"，根据你的降维方法
    pt.size = 1,
    cols = c("lightgrey", "red")  # 颜色从低表达到高表达
  )
}

















