#write by Sun Haiayng at 2025.01.15
###分析来自Mouse TBLCs的PSC
#01、搭环境
{
  rm(list=ls())
  for (i in c("ReactomePA","ggplot2","Seurat","SingleR","celldex",
              "SeuratObject","stringr","ChIPseeker","DESeq2","clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db","TxDb.Mmusculus.UCSC.mm10.knownGene")) {library(i, character.only = T)};rm(i)
}
################################################################################
#02、列索引
{
  path = "/ChIP_seq_2/aaSHY/ghyuu/ask/scPSC/"
  prex = "scPSC"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/crREf_Run/chrNo"
  inx4_iso = "/Reference/aaSHY/AppTools/cellR_assemIso4/chrNo"
  inx5 = "/Reference/aaSHY/AppTools/scTE/Data/shy.exclusive.idx"
}
################################################################################
#03、跑命令
{
  shel = paste0(path, "src/run_a.sh")
  cat("#!/bin/bash\n", append = F, file=shel)
  #跑cellranger
  cmd_01 <- paste0("#parallel-fastq-dump ...","\n")
  cmd_02 <- paste0("#cellranger count --id=samp_1 --transcriptome=", inx4," ",
                   "--fastqs ",path,"fq --sample ",prex," ",
                   "--create-bam=true --localmem=80 --localcores=4 ",
                   "--output-dir ",path,"cellR_Run/",prex,"\n")
  cmd_03 <- paste0("cellranger count --id=samp_1 --transcriptome=", inx4_iso," ",
                   "--fastqs ",path,"fq --sample ",prex," ",
                   "--create-bam=false --localmem=50 --localcores=4 ",
                   "--output-dir ",path,"cellR_Iso/",prex,"\n")
  cat(cmd_01, append = T, file = shel)
  cat(cmd_02, append = T, file = shel)
  cat(cmd_03, append = T, file = shel)
  print(paste0("nohup bash ",shel, " >",path,"Log/",basename(shel),".log"," 2>&1 &"))
  #
  #
  #
  #
  #
  #跑scTE (在scTE环境下)
  for (jjjj in prex) {
    cmd_03 <- paste0("samtools view -h ",path,
                     "cellR_Run/",jjjj,"/outs/possorted_genome_bam.bam |grep ",
                     "-P '^@|CB:' |grep -P '^@|UB:' |samtools view -bS >",path,
                     "cellR_Run/",jjjj,"/outs/possorted_genome_bam.filt.bam","\n")
    cmd_04 <- paste0("samtools index ",path,
                     "cellR_Run/",jjjj,"/outs/possorted_genome_bam.filt.bam","\n")
    cmd_05 <- paste0("mkdir -p ",path,"scTE_Run/",jjjj,"\n")
    cmd_06 <- paste0("cd ",path,"scTE_Run/",jjjj,"\n")
    cmd_07 <- paste0("scTE -p 4 -i ",path,"cellR_Run/",jjjj,
                     "/outs/possorted_genome_bam.filt.bam -x ",inx5," -CB CB ",
                     "-UMI UB -o ",jjjj,"\n")
  }
}
################################################################################
#04、Seurat
{
  red1 <- as.data.frame(fread(paste0(path,"scTE_Run/scPSC/scPSC.csv"), sep=",", header = T))
  rownames(red1) <- red1[,1]
  red1 <- red1[,-1]
  red1 <- t(red1)
  red1 <- as(as.matrix(red1), "dgCMatrix")
  red1 <- Seurat::CreateSeuratObject(counts=red1, min.cells=3)
  red1@meta.data$group <- prex[1]
  objs_1 <- red1
  #
  #
  #
  ###rmGG <- c("Malat1", "Gm42418", "Gm26917", "AY036118")
  ###objs_1 = subset(objs_1, features = setdiff(rownames(objs_1),rmGG))
  objs_1 <- subset(objs_1, subset = nCount_RNA > 4000 & nCount_RNA < 40000)
  rmMI <- rownames(objs_1)[grep("^mt",rownames(objs_1))]
  objs_1[["percMito"]] = PercentageFeatureSet(objs_1, features = rmMI)
  objs_1 = subset(objs_1, subset = percMito < 5)
  #
  options(future.globals.maxSize = 32 * 1e9)
  objs_2 <- objs_1
  objs_2 <- Seurat::NormalizeData(objs_2)
  objs_2 <- Seurat::FindVariableFeatures(objs_2, selection.method="vst", nfeatures=3000)
  objs_2 <- Seurat::ScaleData(objs_2)
  #objs_2 <- JoinLayers(objs_2)
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
                            group.by = c("seurat_clusters"),alpha = 0.8,
                            cols = c(ggsci::pal_d3(palette = "category20")(8)),
                            label = T, label.size = 5, label.color = "black",
                            pt.size = .5)
  p___01
}
################################################################################
#04、可视化feature的表达分布
{
  aims <- c("MM","MM-int")
  RidgePlot(objs_4, features = aims, ncol = 2, group.by ="seurat_clusters")
  VlnPlot(objs_4, features = aims, ncol = 2, group.by = "seurat_clusters", pt.size = 0)
  FeaturePlot(objs_4, features = aims)
  #
  #
  #
  FeaturePlot(objs_4, features = "ORR1A3-int")
  #
  #
  #
  DotPlot(objs_4, features = aims, group.by = "seurat_clusters") + RotatedAxis()
  DoHeatmap(objs_4, features = VariableFeatures(objs_4)[1:10], size = 3)
}
{
  aims <- c("TP53","CDKN2A","CDKN1A","IFNB1","IFNG","IFNAR1")
  aims <- babelgene::orthologs(genes = aims, species = "mouse", human = T)$symbol
  FeaturePlot(objs_4, features = aims)
}




