my_DESeq2 <- function(counts,
                      group,
                      Rdata.name,
                      csv.name,
                      vst.name,
                      org = "org.Rn.eg.db"){
  # meta.data
  meta.data <- group
  
  ### DESeq2
  ### 要有数据countData，这里就是exprSet
  ### 要有分组信息，在colData中，这里是metadata
  ### 第一列如果是基因名称，需要自动处理，设置参数tidy=TRUE
  ### design部分是分组信息，格式是~group
  dds <-DESeqDataSetFromMatrix(countData = counts, 
                               colData = meta.data, 
                               design = ~ group,
                               tidy = F)
  nrow(dds)
  rownames(dds)
  
  ### 筛选样本，counts函数提取表达量数据
  dds <- dds[rowSums(counts(dds))>1,]
  nrow(dds)
  
  #########################################################
  ### 4.数据质量判断
  ### vst标准化处理
  vsd <- vst(dds, blind = FALSE)
  ### 内置函数plotPCA进行主成分分析画图
  plotPCA(vsd, "group")
  
  ### 用内置函数plotCounts来进行快速，简易作图
  ### 找到阳性基因,此处ESR1
  ### dds来自于上一步
  ### gene 输入的是ensemble ID
  ### intgroup 输入的是metadata中分组信息
  #plotCounts(dds, gene = "Piezo1", intgroup=c("group"))
  
  ### 导出标准化后的表达数据
  ### assay函数提取vst标准化后的数据，保存数据用于热图
  exprSet_vst <- as.data.frame(assay(vsd))
  write.csv(exprSet_vst, file = vst.name)
  
  ### 5.正式运行DESeq主程序
  dds <- DESeq(dds)
  ### 依次是，1.分组信息(metadata中的列) 2.处理组，3.对照组
  contrast=c("group", 
             levels(meta.data$group)[1], 
             levels(meta.data$group)[2])
  ### results函数获取差异分析的结果
  dd1 <- results(dds, contrast=contrast, alpha = 0.05)
  ### 内置函数plotMA作图
  plotMA(dd1, ylim=c(-5,5))
  ### logFC矫正
  dd2 <- lfcShrink(dds,contrast=contrast, res=dd1,type="ashr")
  plotMA(dd2, ylim=c(-5,5))
  ### 7.导出差异分析的结果
  library(dplyr)
  library(tibble)
  library(tidyr)
  res <- dd2 %>% 
    as.data.frame() %>% 
    rownames_to_column("gene_id") 
  
  ### 8.基因注释
  ### 当前基因名称是ENSEMBL
  library(AnnotationDbi)
  library(org.Rn.eg.db)
  ### 增加基因名称SYMBOL
  res$entrz <- mapIds(org.Rn.eg.db,
                      keys=res$gene_id,
                      column="ENTREZID",
                      keytype="SYMBOL",
                      multiVals="first")
  ### 修改logFC，p值的名称，为的是跟火山图的代码匹配
  colnames(res) <- c("symbol","baseMean","logFC","lfcSE","P.Value","adj.P.Val","entrez")
  
  ### 9.保存数据
  save(res, file = Rdata.name)
  write.csv(res,file = csv.name)
}
