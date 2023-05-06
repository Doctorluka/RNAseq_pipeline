#########################################################
###################### step1: 准备数据 ###################
#########################################################

# 以R project为起点
# R project下建立文件夹：
# data：储存Rdata/RDS等文件
# figure：储存图片
# result：储存结果文件，如csv
# scripts：储存脚本代码

#常规操作
rm(list = ls())
options(stringsAsFactors = F)
gc()

#读入数据
library(openxlsx) #excel文件读入
exp <- read.xlsx("data/ref_trans_full_table.xlsx", sheet = 2) #原始数据放在data文件夹里
str(exp) #查看数据是否正确
# 'data.frame':	29316 obs. of  7 variables:
# $ gene_name: chr  "Gad1" "Alx4" "Tmco5b" "Cbln1" ...
# $ A71_count: chr  "0" "0" "0" "0" ...
# $ A72_count: chr  "13" "3" "0" "9" ...
# $ A73_count: chr  "0" "4" "0" "4" ...
# $ I85_count: num  31 0 0 4 11 ...
# $ I93_count: num  17 0 0 0 7 ...
# $ I99_count: num  11 0 0 2 8 ...

#构建分组数据
group <- data.frame(sample = colnames(exp)[2:7],
                    group = c(rep("Ctrl",3), rep("treat",3)))
head(group)
#      sample group
# 1 A71_count  Ctrl
# 2 A72_count  Ctrl
# 3 A73_count  Ctrl
# 4 I85_count treat
# 5 I93_count treat
# 6 I99_count treat

save(exp,group, file = "data/1.rawcounts_group.Rdata") #保存在合适的文件夹



#########################################################
###################### step2: 过滤基因 ###################
#########################################################

rm(list = ls())
gc()

library(dplyr)
library(tibble)
library(tidyr)

load("data/1.rawcounts_group.Rdata")

# 过滤重复基因和无表达量基因
filter_counts <- exp %>% 
  dplyr::mutate(rowmean = rowMeans(.[,2:ncol(.)])) %>%  # 增加一列计算行平均值(注意：这里第一列是基因名)
  dplyr::arrange(desc(rowmean)) %>%  # 按rowmeans倒序排列
  dplyr::distinct(gene_name, .keep_all = T) %>%  # 去除重复的基因行，保留表达量最高的那个
  dplyr::filter(rowmean > 0) %>%  # 去除表达量为0的基因(可以根据自己要求调整)
  dplyr::select(-rowmean) %>% # 去除rowmeans列
  na.omit() # 去除NA值
dim(filter_counts)
# [1] 17625     7    # 最终应获得约16000-18000左右的基因比较常见

filter_counts[1:4,1:4] #最后矩阵长这样
#   gene_name A71_count A72_count A73_count
# 1      Lyz2    624681    505091    536698
# 2    Mt-co1    341124    433895    251043
# 3     Sftpc    348630    311728    329555
# 4    Sftpa1    162027    161295    140403

group #查看样品名与group是否对应
#      sample group
# 1 A71_count  Ctrl
# 2 A72_count  Ctrl
# 3 A73_count  Ctrl
# 4 I85_count treat
# 5 I93_count treat
# 6 I99_count treat

#########################################################
###################### step3: DESeq2 ####################
#########################################################

library(DESeq2)
dds <-DESeqDataSetFromMatrix(countData=filter_counts, 
                             colData=group_CM, 
                             design=~group,
                             tidy=TRUE) #TRUE是指filter_counts第一列是基因名，如果基因名在行名上应该选FALSE
nrow(dds)
rownames(dds)


### 筛选样本，counts函数提取表达量数据
dds <- dds[rowSums(counts(dds))>1,]
nrow(dds)

### 4.数据质量判断
### vst标准化处理
vsd <- vst(dds, blind = FALSE)
### 内置函数plotPCA进行主成分分析画图，可以用ggplot2语法修改
plotPCA(vsd, "group")

### 用内置函数plotCounts来进行快速，简易作图
### 找到阳性基因,此处ESR1
### dds来自于上一步
### gene 输入的是ensemble ID
### intgroup 输入的是metadata中分组信息
plotCounts(dds, gene = "Piezo1", intgroup=c("group"))

### 导出标准化后的表达数据
### assay函数提取vst标准化后的数据，保存数据用于热图
exprSet_vst <- as.data.frame(assay(vsd))
test <- exprSet_vst %>% head() #查看数据

### 保存数据,用于表达量作图，比如差异分析，热图
save(exprSet_vst,file = "data/2.exprSet_vst.Rdata")


### 5.正式运行DESeq主程序
dds <- DESeq(dds)

### 6.logFC矫正，RNAseq很重要的一步
### contrast参数设置
### 依次是，1.分组信息(metadata中的列) 2.处理组，3.对照组
contrast=c("group", "MCT", "Ctrl")

### results函数获取差异分析的结果
dd1 <- results(dds, contrast=contrast, alpha = 0.05)

### 内置函数plotMA作图
plotMA(dd1, ylim=c(-5,5))

### logFC矫正
dd2 <- lfcShrink(dds,contrast=contrast, res=dd1,type="ashr")
plotMA(dd2, ylim=c(-5,5))


### 7.导出差异分析的结果
res <- dd2 %>% 
  as.data.frame() %>% 
  rownames_to_column("gene_id") 

#### 增加ENTREZID,用于KEGG分析
# 根据数据的实际情况调整，这里的目的就是获得symbol和entrezid两列方便后续分析
library(clusterProfiler)
library(org.Rn.eg.db)
entrez <- bitr(geneID = res$gene_id,
               fromType = "SYMBOL",
               toType = "ENTREZID",
               OrgDb = "org.Rn.eg.db") #这里是大鼠的注释信息,人的"org.Hs.eg.db",小鼠的"org.Mm.eg.db"
res <- res %>% inner_join(entrez, ., by =c("SYMBOL" = "gene_id"))
dim(res)

### 修改logFC，p值的名称，为的是跟火山图的代码匹配
colnames(res) <- c("gene_id","entrez","baseMean","logFC","lfcSE","P.Value","adj.P.Val")

head(res)
#   gene_id entrez baseMean       logFC      lfcSE      P.Value    adj.P.Val
# 1    Lyz2  25211 566390.1  0.01583649 0.06252924 6.078288e-01 0.8905920100
# 2   Sftpc  50683 322508.1 -0.01971806 0.06582212 5.395627e-01 0.8589446657
# 3  Sftpa1  24773 168197.6  0.15176507 0.11144402 1.656068e-02 0.1414188917
# 4  Eef1a1 171361 152705.7 -0.07807301 0.09461995 9.130075e-02 0.3987862937
# 5     Gsn 296654 153036.0 -0.41776936 0.23682222 1.775595e-03 0.0291449051
# 6   Epas1  29452 136472.1 -0.84452321 0.21459858 2.079157e-06 0.0001354173

# 保存差异分析结果
save(res,file = "data/2.res.Rdata")


#########################################################
################# step3: DESeq2偷懒法 ####################
#########################################################
# 重要前提！！
# 必须走完一遍上面的流程，熟悉各个数据的格式，熟悉之后可以按下面流程走

# 首先按照自建函数方法保存 my_DESeq2.R 文件到scripts文件夹
source("scripts/my_DESeq2.R")

# 不要运行，这是示例的默认参数
my_DESeq2 <- function(counts, # 格式见下方
                      group, # 格式见下方
                      Rdata.name, # DESeq2分析结果
                      csv.name,  # DESeq2分析结果
                      vst.name, # 标准化表达矩阵
                      org = "org.Rn.eg.db") # 这里是大鼠的注释信息,人的"org.Hs.eg.db",小鼠的"org.Mm.eg.db"
# 各输入参数 参考格式
counts
#   gene_name A71_count A72_count A73_count
# 1      Lyz2    624681    505091    536698
# 2    Mt-co1    341124    433895    251043
# 3     Sftpc    348630    311728    329555
# 4    Sftpa1    162027    161295    140403

group #查看样品名与group是否对应
#      sample group
# 1 A71_count  Ctrl
# 2 A72_count  Ctrl
# 3 A73_count  Ctrl
# 4 I85_count treat
# 5 I93_count treat
# 6 I99_count treat

# 运行
my_DESeq2(counts = exp, group = group,
          Rdata.name = "results/3.DESeq2.Rdata",
          csv.name = "results/3.DESeq2_res.csv",
          vst.name = "results/3.exprset_vst.csv")
# 加载运行结果，与上面常规分析结果一致
load("results/3.DESeq2.Rdata")
exprSet_vst <- read.csv("results/3.exprset_vst.csv", row.names = 1)



#########################################################
################# step4: volcano plot ###################
#########################################################

# 方法一：自己敲代码画
# 仅提供示例代码，需要自定义的自己上网找

ggplot(res, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(alpha=0.5, size=3.5, 
             aes(color = regulated))+
  ylab("-log10(Pvalue)") +
  scale_color_manual(values=c("blue", "grey","#D31638")) +
  geom_vline(xintercept=c(-logFC_cutoff,logFC_cutoff),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(pvalue_cutoff),lty=4,col="black",lwd=0.8) +
  xlim(-5, 5)+
  theme_bw()+
  theme(legend.justification = c(1,0), 
        legend.position = c(1,0),
        legend.background = element_rect(fill = 'white', colour = 'black'))+
  coord_flip()

# 方法二：偷懒法1   #PS:长什么样子见 README.md
# 首先保存 my_DESeq2.R 文件到scripts文件夹
source(file = "scripts/function.volcanoplot.1.R")

# 示例参数
my_volcanoplot.1 <- function(diff,  # DESeq2分析结果，并按照上述方法改了列名的
                             logFC_cutoff = 0.58, # logFC筛选阈值,默认0.58, 即FC=1
                             pvalue_cutoff = 0.05, # p值筛选阈值
                             label = FALSE, # 是否添加基因名
                             label.n = 7, # 添加最显著的基因名，默认上下调各前7个，label = TRUE时有效
                             label.color = "#A58533", # 基因名颜色，label = TRUE时有效
                             color = c("blue", "grey","#D31638"), # 颜色，依次是下调、无差异、上调基因
                             save = FALSE, # 是否保存图片
                             filename = NULL, # 保存图片的名称，save = TRUE时有效
                             save.width = 6, # 保存图片的宽和长，save = TRUE时有效
                             save.height = 6)

head(res)
#   gene_id entrez baseMean       logFC      lfcSE      P.Value    adj.P.Val
# 1    Lyz2  25211 566390.1  0.01583649 0.06252924 6.078288e-01 0.8905920100
# 2   Sftpc  50683 322508.1 -0.01971806 0.06582212 5.395627e-01 0.8589446657
# 3  Sftpa1  24773 168197.6  0.15176507 0.11144402 1.656068e-02 0.1414188917
# 4  Eef1a1 171361 152705.7 -0.07807301 0.09461995 9.130075e-02 0.3987862937
# 5     Gsn 296654 153036.0 -0.41776936 0.23682222 1.775595e-03 0.0291449051
# 6   Epas1  29452 136472.1 -0.84452321 0.21459858 2.079157e-06 0.0001354173

my_volcanoplot(data, label = T, label.n = 10,
               save = T, filename = "figure/2.volcano.pdf")

# 方法三：偷懒法2
# 自建函数-2
source(file = "scripts/function.volcanoplot.2.R")

# 示例参数
volcanoPlotter <- function(data, label.gene) #只提供标记基因

# 默认标记前30，可自定义
volcanoPlotter(res,label.gene = NULL)

# labels defined by yourself
label.gene.inf <- c("Nfkbia","Bcr","Cdkn1a","Cdk1","Cdkn3","Foxs1",
                    "Foxc2","Mmp14","Mmp2","Mmp23","Mmp12","Melk",
                    "Mmp8","Mmp9","Timp1","Pdgfc","Pdgfrl","Ets1",
                    "Angptl7","Tgfbr3","Tgfbr1","Tgfb3","Tgfbr2",
                    "Tgfb1","Ctgf","Pcna")
volcanoPlotter(res, label.gene = label.gene.inf)



#########################################################
################# step4: heatmap ########################
#########################################################

# heatmap
library(pheatmap)

# 自己选择基因，或者直接用差异基因
label.gene.inf <- c("Nfkbia","Bcr","Cdkn1a","Cdk1","Cdkn3","Foxs1",
                    "Foxc2","Mmp14","Mmp2","Mmp23","Mmp12","Melk",
                    "Mmp8","Mmp9","Timp1","Pdgfc","Pdgfrl","Ets1",
                    "Angptl7","Tgfbr3","Tgfbr1","Tgfb3","Tgfbr2",
                    "Tgfb1","Ctgf","Pcna")

heatdata <- exprSet_vst[label.gene.inf,] # 用前面标准化的矩阵
head(heatdata) # 这里的示例数据跟前面的有所不同，不用在意
#             Ctrl      Ctrl      Ctrl     treat     treat     treat
# Nfkbia 11.958522 12.229809 12.024335 13.074980 12.281264 12.892269
# Bcr    11.403141 11.563972 11.421307 11.336955 10.803347 11.158113
# Cdkn1a 11.572151 11.856884 11.214111 10.924862 10.867515 10.959334
# Cdk1   10.108383  9.818906 10.103609  9.279739  9.320314  9.182204
# Cdkn3   8.976968  9.183561  9.198010  8.673351  8.568731  8.646242
# Foxs1   9.641579  9.513576  9.413421  9.018259  9.087003  8.967537

#注释信息
annotation_col <- data.frame(group = c(rep("Ctrl",3), rep("treat",3)))
rownames(annotation_col) <- colnames(heatdata) #要跟表达矩阵的列名(sample)的分组对应
annotation_col 
#   group
# 1  Ctrl
# 2  Ctrl
# 3  Ctrl
# 4 treat
# 5 treat
# 6 treat

pheatmap(heatdata, #热图的数据
         cluster_rows = TRUE,#行聚类
         cluster_cols = F,#列聚类，可以看出样本之间的区分度
         annotation_col =annotation_col, #标注样本分类
         annotation_names_col = F, #不显示注释名称
         annotation_legend=TRUE, # 显示注释
         show_rownames = T,# 显示行名
         show_colnames = F,# 显示列名
         scale = "row", #以行来标准化，这个功能很不错
         color =colorRampPalette(c("#68A5EE", "white","#EC5B22"))(100),#调色
         #filename = "heatmap_F.pdf",#是否保存
         cellwidth = 25, cellheight = 12,# 格子比例
         fontsize = 10) #字体大小
library(export)
graph2pdf(file = "figure/3heatmap_select.pdf")


#########################################################
################# step4: GSEA ###########################
#########################################################

# 个人喜欢GSEA
# 个人现在几乎不用差异基因进行富集分析了，所以没有常规富集分析代码

library(msigdbr)
library(clusterProfiler)
library(stringr)
library(ggplot2)
library(export)

### 1.获取基因logFC
geneList <- res$logFC
### 2.命名
names(geneList) = res$symbol
## 3.排序很重要
geneList = sort(geneList, decreasing = TRUE)
head(geneList)

#GSEA GO gene set
rat_df = msigdbr(species = "Rattus norvegicus") #这里是大鼠的学名，小鼠是 Mus musculus，人是 Homo sapiens(函数默认是人)
colnames(rat_df)
rat_df %>% 
  distinct(gs_cat, gs_subcat) %>% 
  arrange(gs_cat, gs_subcat) %>% 
  print(n = 23)

# 提取KEGG数据库信息,可以自定义，详情查询GSEA官网
KEGG_df <- msigdbr(species = "Rattus norvegicus", 
                 category = "C2", 
                 subcategory = "CP:KEGG")
KEGG_df <- KEGG_df %>% 
  dplyr::select(gs_name, gene_symbol)
head(KEGG_df)

set.seed(111) # 每次运行结果有所差别，设置随机种子
KEGG.res <- GSEA(geneList,TERM2GENE = KEGG_df)

# 画图，图就不展示了，很常规
# 自定义美化自己查代码
# GSEA结果用不了barplot()
dotplot(KEGG.res,
        color = "pvalue",
        showCategory=10,
        label_format = 60,
        split=".sign",
        font.size = 10)+facet_grid(~.sign)
# 保存信息
yd <- as.data.frame(KEGG.res)
write.csv(yd, file = "results/4.KEGG_Res.csv")


# 提取Wiki pathway数据库信息
WP_df <- msigdbr(species = "Rattus norvegicus", 
                 category = "C2", 
                 subcategory = "CP:WIKIPATHWAYS")
WP_df <- WP_df %>% 
  dplyr::select(gs_name, gene_symbol)
head(WP_df)

set.seed(111)
go.res <- GSEA(geneList,TERM2GENE = WP_df)
