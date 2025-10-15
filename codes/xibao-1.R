#write by Sun Haiayng at 2024.12.18
#01、搭环境
{
  rm(list=ls())
  for (i in c("ReactomePA","ggplot2","Seurat","SingleR","celldex",
              "SeuratObject","stringr","ChIPseeker","DESeq2","clusterProfiler",
              "topGO","Rgraphviz","pathview","org.Mm.eg.db","TxDb.Mmusculus.UCSC.mm10.knownGene")) {library(i, character.only = T)};rm(i)
}
#02、给索引
{
  path = "/disk5/aaSHY/scAging/"
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/fa/GRCm38.dna.primary_assembly.fa"
  inx2 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx4 = "/Reference/aaSHY/AppTools/cellR_assemIso4/chrNo"
  #inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/cellR_Run/chrNo"
  inx5 = "/Reference/aaSHY/AppTools/scTE/Data/shy.exclusive.idx"
  inxA = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MERVL-int::ERVL::LTR.bed"
  inxB = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MT2_Mm::ERVL::LTR.bed"
  inxC = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/IAPEY_LTR::ERVK::LTR.bed"
  inxD = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/RLTR1B-int::ERV1::LTR.bed"
}
################################################################################
#03、建索引、文件名、核验
{
  shell <- paste0(path,"src/run_a.sh")#; shell <- paste0(path,"src/run_b.sh")
  cat("#!/bin/bash\n", file = shell)
  cmd_01 <- paste0("cd /Reference/aaSHY/AppTools/cellR_assemIso4/; ",
                   "cellranger mkref --genome=GRCm38 ",
                   "--fasta=GRCm38.dna.primary_assembly.fa ",
                   "--genes=pacbFlt-correct-2.gtf --nthreads=10 --memgb=64","\n")
  for (i in cmd_01) {cat(i, file = shell, append = T)}
  #cmd_02 <- paste0("#scTE_build -te /Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38.TE.bed -gene /Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf -g mm10 -o /Reference/aaSHY/AppTools/scTE/Data/shy -m inclusive")
  #cmd_03 <- paste0("#scTE_build -te /Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38.TE.bed -gene /Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf -g mm10 -o /Reference/aaSHY/AppTools/scTE/Data/shy -m exclusive")
  #for (i in cmd_02) {cat(i, file = shell, append = T)}
  #for (i in cmd_03) {cat(i, file = shell, append = T)}
  #red1 <- read.table("/disk5/aaSHY/scAging/fq/sampleName_clientId.txt",sep = "\t",header = T)
  #file.rename(from = paste0(path,"fq/",red1$sampleName,"_good_1.fq.gz"),to = paste0(path,"fq/",red1$clientId,"_S",seq(9),"_L001_R1_001.fastq.gz"))
  #file.rename(from = paste0(path,"fq/",red1$sampleName,"_good_2.fq.gz"),to = paste0(path,"fq/",red1$clientId,"_S",seq(9),"_L001_R2_001.fastq.gz"))
  #ff01 <- list.files(paste0(path,"fq"),"fastq.gz",full.names = T)
  #cmd_04 <- paste0(">",path,"fq/shy_md5.txt","\n")
  #cmd_05 <- paste0("md5sum ",ff01," >>",path,"fq/shy_md5.txt","\n")
  #for (i in cmd_04) {cat(i, file = shell, append = T)}
  #for (i in cmd_05) {cat(i, file = shell, append = T)}
  print(paste0("nohup bash ",shell, " >",paste0(path,"Log/"),basename(shell),".log"," 2>&1 &"))
}
################################################################################
#04、跑命令 [cellranger]
{
  #注意索引的变化
  #conda activate base
  shell <- paste0(path,"src/run_1.sh")
  cat("#!/bin/bash\n", file = shell)
  prex <- list.files(paste0(path,"fq"),pattern = "_R1_001.fastq.gz")
  prex <- sapply(str_split(prex,"_"),"[",1)
  cmd_01 <- paste0("mkdir -p ",path,"cellR_Iso/",prex,"\n")
  cmd_02 <- paste0("cellranger count --id=samp_1 --transcriptome=", inx4," ",
                   "--fastqs ",path,"fq --sample ",prex," ",
                   "--create-bam=true --localmem=80 --localcores=4 ",
                   "--output-dir ",path,"cellR_Iso/",prex,"\n")
  for (i in cmd_01) {cat(i, file = shell, append = T)}
  for (i in cmd_02) {cat(i, file = shell, append = T)}
  print(paste0("nohup bash ",shell, " >",paste0(path,"Log/"),basename(shell),".log"," 2>&1 &"))
}
################################################################################
#05、跑命令 [scTE]
{
  #conda activate scTE
  shell <- paste0(path,"src/run_c.sh")
  cat("#!/bin/bash\n", file = shell)
  prex <- list.files(paste0(path,"fq"),pattern = "_R1_001.fastq.gz")
  prex <- sapply(str_split(prex,"_"),"[",1)
  prex <- prex
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
                     "-UMI UB --min_genes 200 -o ",jjjj,"\n")
    for (i in cmd_03) {cat(i, file = shell, append = T)}
    for (i in cmd_04) {cat(i, file = shell, append = T)}
    for (i in cmd_05) {cat(i, file = shell, append = T)}
    for (i in cmd_06) {cat(i, file = shell, append = T)}
    for (i in cmd_07) {cat(i, file = shell, append = T)}
  }
  print(paste0("nohup bash ",shell, " >",paste0(path,"Log/"),basename(shell),".log"," 2>&1 &"))
}
################################################################################
#
#
#
#
#
#
#
#
#
#做测试--------
#做测试--------
#做测试--------
################################################################################
#06、UMAP [01] [使用的是cellranger的结果] [Integrate整合样本]
{
  #options(Seurat.object.assay.version = "v5")
  tmp1 <- Seurat::Read10X(paste0(path,"cellR_copy/Bcon"))
  tmp2 <- Seurat::Read10X(paste0(path,"cellR_copy/Bhucp"))
  tmp3 <- Seurat::Read10X(paste0(path,"cellR_copy/Bhucs"))
  abj1 <- Seurat::CreateSeuratObject(counts = tmp1, min.cells = 3,min.features = 200); dim(abj1)#
  abj2 <- Seurat::CreateSeuratObject(counts = tmp2, min.cells = 3,min.features = 200); dim(abj2)#
  abj3 <- Seurat::CreateSeuratObject(counts = tmp3, min.cells = 3,min.features = 200); dim(abj3)#
  objs <- list(abj1,abj2,abj3)
  rmGG <- c("Malat1", "Gm42418", "Gm26917", "AY036118")
  objs <- lapply(objs, function(x) {
    x = subset(x, features = setdiff(rownames(x),rmGG))
    rmMI <- rownames(x)[grep("^mt",rownames(x))]
    x[["percMito"]] = PercentageFeatureSet(x, features = rmMI)
    x = subset(x, subset = percMito < 25)
  })
  #
  options(future.globals.maxSize = 32 * 1e9)
  objs <- lapply(objs, function(x) {
    x <- Seurat::NormalizeData(x)
    x <- Seurat::FindVariableFeatures(x, selection.method="vst", nfeatures=2000)
    x <- Seurat::ScaleData(x)
    x <- Seurat::SCTransform(x)
  })
  #
  objs[[1]]@meta.data$group <- "Bcon"
  objs[[2]]@meta.data$group <- "Bhucp"
  objs[[3]]@meta.data$group <- "Bhucs"
  #
  scFe <- Seurat::SelectIntegrationFeatures(object.list = objs, nfeatures = 3000)
  objs <- Seurat::PrepSCTIntegration(object.list = objs, anchor.features = scFe)
  future::plan("multisession", workers=4)
  anch <- Seurat::FindIntegrationAnchors(object.list = objs, 
                                         normalization.method = "SCT",
                                         anchor.features = scFe)
  anch_I <- Seurat::IntegrateData(anchorset=anch, normalization.method = "SCT")
  future::plan("sequential")
  #
  oooo <- Seurat::RunPCA(anch_I)
  oooo <- Seurat::FindNeighbors(oooo)
  oooo <- Seurat::FindClusters(oooo, resolution = 0.2)
  #
  umap <- Seurat::RunUMAP(oooo, dims = 1:20)
  #umap@meta.data$group <- factor(umap@meta.data$group, levels = c("Bcon","Bhucs","Bhucp"))
  umap_1 <- subset(umap, cells = WhichCells(umap, expression = group == "Bhucp"))
  uuuu <- Seurat::DimPlot(umap_1, reduction = "umap",
                          group.by = "group",alpha = 0.2,
                          #cols = ggsci::pal_d3(palette = "category20")(3),
                          #label = T, label.size = 5, label.color = "black",
                          pt.size = .5)
  uuuu
  ###tsne <- Seurat::RunTSNE(objs, dims = 1:20)
  ###tttt <- Seurat::DimPlot(tsne, reduction = "tsne", group.by = "orig.ident",
  ###                        cols = pal_d3(palette = "category20")(2),
  ###                        #label = T, label.size = 5, label.color = "black",
  ###                        pt.size = 1)
  ###tttt
}
################################################################################
#06、UMAP [02] [使用的是cellranger的结果] [merge----整合样本]
{
  #options(Seurat.object.assay.version = "v5")
  tmp1 <- Seurat::Read10X(paste0(path,"cellR_Copy/Ocon"))
  tmp2 <- Seurat::Read10X(paste0(path,"cellR_Copy/Ohucp"))
  tmp3 <- Seurat::Read10X(paste0(path,"cellR_Copy/Ohucs"))
  abj1 <- Seurat::CreateSeuratObject(counts = tmp1, min.cells = 3,min.features = 200); dim(abj1)#
  abj2 <- Seurat::CreateSeuratObject(counts = tmp2, min.cells = 3,min.features = 200); dim(abj2)#
  abj3 <- Seurat::CreateSeuratObject(counts = tmp3, min.cells = 3,min.features = 200); dim(abj3)#
  abj1@meta.data$group <- "Bcon"
  abj2@meta.data$group <- "Bhucp"
  abj3@meta.data$group <- "Bhucs"
  objs <- merge(abj1,y = list(abj2,abj3), add.cell.ids = c("Bcon","Bhucp","Bhucs"))
  View(objs_2)
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
  ###objs_4@meta.data$group <- factor(objs_4@meta.data$group, levels = c("Bcon","Bhucs","Bhucp"))
  ###objs_4 <- subset(objs_4, cells = WhichCells(objs_4, expression = group == "Bhucp"))
  p___01 <- Seurat::DimPlot(objs_4, reduction = "umap",
                            group.by = c("seurat_clusters"),alpha = 0.2,
                            cols = ggsci::pal_d3(palette = "category20")(18),
                            #label = T, label.size = 5, label.color = "black",
                            pt.size = .5)
  p___01
  objs_4@meta.data$group <- factor(objs_4@meta.data$group, levels = c("Bcon","Bhucs","Bhucp"))
  p___02 <- Seurat::DimPlot(objs_4, reduction = "umap",
                            group.by = c("group"),alpha = 0.2,
                            cols = ggsci::pal_d3(palette = "category20")(3),
                            #label = T, label.size = 5, label.color = "black",
                            pt.size = .5)
  p___02
  #
  #
  #
  objs_5 <- subset(objs_4, cells = WhichCells(objs_4, expression = group != "Bhucp"))
  p___03 <- Seurat::DimPlot(objs_5, reduction = "umap",
                            group.by = c("group"),alpha = 0.3,
                            cols = c("steelblue","indianred"),
                            #label = T, label.size = 5, label.color = "black",
                            pt.size = .5)
  p___03
  objs_6 <- subset(objs_4, cells = WhichCells(objs_4, expression = group != "Bhucs"))
  p___04 <- Seurat::DimPlot(objs_6, reduction = "umap",
                            group.by = c("group"),alpha = 0.3,
                            cols = c("steelblue","indianred"),
                            #label = T, label.size = 5, label.color = "black",
                            pt.size = .5)
  p___04
  objs_7 <- subset(objs_4, cells = WhichCells(objs_4, expression = group == "Bcon"))
  p___05 <- Seurat::DimPlot(objs_7, reduction = "umap",
                            group.by = c("group"),alpha = 0.3,
                            cols = c("steelblue"),
                            #label = T, label.size = 5, label.color = "black",
                            pt.size = .5)
  p___05
}
################################################################################
#07、UMAP [03] [使用的是scTE的结果------] [merge----整合样本]
{
  #options(Seurat.object.assay.version = "v5")
  prex <- c("Bcon","Bhucs","Bhucp")
  prex <- c("Mcon","Mhucs","Mhucp")
  prex <- c("Ocon","Ohucs","Ohucp")
  for (op in seq_along(prex)) {
    print(op)
    dirr <- paste0(path,"scTE/",prex[op],"/",prex[op],".csv")
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
  ###objs_4 <- Seurat::RunTSNE(objs_4, dims = 1:10)
  ###objs_4@meta.data$group <- factor(objs_4@meta.data$group, levels = prex)
  ###objs_4 <- subset(objs_4, cells = WhichCells(objs_4, expression = group == prex[3]))
  p___01 <- Seurat::DimPlot(objs_4, reduction = "umap",
                            group.by = c("seurat_clusters"),alpha = 0.2,
                            cols = c(ggsci::pal_d3(palette = "category20")(20),"grey"),
                            label = T, label.size = 5, label.color = "black",
                            pt.size = .5)
  p___01
  #
  #
  #
  p___02 <- ggplot(data = objs_4@meta.data) + 
    geom_bar(aes(x = seurat_clusters, fill = group), position = "fill") +
    scale_fill_jama()
  p___02
  #
  #
  #
  {
    #每个样本的总细胞数量里有多少比例的细胞分布在cluster里
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
  }
  #
  #
  #
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
  #
  #
  #
  aims <- c("MT2-Mm","MERVL-int")
  RidgePlot(objs_4, features = aims, ncol = 2, group.by ="group")
  VlnPlot(objs_4, features = aims, ncol = 2, group.by ="group", pt.size = 0)
  FeaturePlot(objs_4, features = aims)
  DotPlot(objs_4, features = aims, group.by = "group") + RotatedAxis()
  DotPlot(objs_4, features = aims, group.by = "seurat_clusters") + RotatedAxis()
  DoHeatmap(objs_4, features = VariableFeatures(objs_4)[1:10], size = 3)
}







