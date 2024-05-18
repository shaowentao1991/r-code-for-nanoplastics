#载入seurat包，并看一下umap图
install.packages("Seurat")
library(Seurat)
DimPlot(project,reduction = "umap")

#载入自动调颜色的包
install.packages("ggsci")
library(ggsci)
library(ggplot2)
#载入umap图小坐标轴的包
install.packages("tidydr")
library(tidydr)
###总览####
umapplot2<-umapplot1+scale_color_npg()  #nature 色系
plot(umapplot2)  #显示图片
#自己定义细胞类型，然后将细胞类型加入到project中，用nature色系作图
clustercelltype1 <- c("1"= "Epithelial cell", 
                      "2"= "Epithelial cell", 
                      "3"= "Immune cell", Neutrophils
                      "4"= "Endothelial cell", 
                      "5"= "Fibroblast cell",
                      "6"= "Fibroblast cell", 
                      "7"= "Fibroblast cell", 
                      "8"= "Epithelial cell",
                      "9"= "Immune cell", NKT 
                      "10"= "Immune cell", Macrophages 
                      "11"= "Fibroblast cell",
                      "12"= "Smooth muscle cell",
                      "13"= "Fibroblast cell",
                      "14"= "Immune cell", DC 
                      "15"= "Epithelial cell",
                      "16"= "Epithelial cell",
                      "17"= "Immune cell", ILC 
                      "18"= "Fibroblast cell",
                      "19"= "Mesothelial cell",
                      "20"= "Smooth muscle cell",
                      "21"= "Immune cell",Macrophages
                      "22"= "Epithelial cell",
                      "23"= "Immune cell",ILC 
                      "24"= "Endothelial cell",
                      "25"= "Immune cell",Macrophages
                      "26"= "Smooth muscle cell",
                      "27"= "Immune cell",
                      "28"= "Immune cell", B cells
                      "29"= "Fibroblast cell",
                      "30"= "Immune cell",DC
                      "31"= "Epithelial cell",
                      "32"= "Epithelial cell"
                      )
project[["cell_type"]]=unname(clustercelltype1[project@meta.data$integrated_snn_res.0.9])

#设置图片中细胞类型的展示顺序
project@meta.data$cell_type<-factor(project@meta.data$cell_type,
                                    levels = c("Epithelial cell","Mesothelial cell",
                                               "Fibroblast cell","Endothelial cell",
                                               "Smooth muscle cell","Immune cell"))
table(project@meta.data$Cluster,project@meta.data$cell_type)
#储存图片时设置为1500*1500
DimPlot(project, reduction = 'umap', group.by = 'cell_type',
        label = FALSE, label.size=8,pt.size = 2)
  scale_color_npg()+
  theme_dr()+ #小坐标轴
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#设置三个样本的展示顺序
samplename <- c("mGBctrl"= "Vehicle",
                "mgb3w"= "3weeks",
                "mgb6w"= "6wkeeks")
project[["samplename1"]]=unname(samplename[project@meta.data$Sample])
project@meta.data$samplename1<-factor(project@meta.data$samplename1,
                                      levels = c("Vehicle","3weeks","6wkeeks"))
#输出像素为2400*800
DimPlot(project, reduction = 'umap', split.by = "samplename1",group.by = 'cell_type',
        label.size=8,pt.size = 2)+
  scale_color_npg() +
  theme_dr()+ #小坐标轴
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#做柱状图,先建立一个数据框
sample_table <- as.data.frame(table(project@meta.data$samplename1,project@meta.data$cell_type))
names(sample_table) <- c("samplename1","cell_type","CellNumber")

#做堆积柱状图,输出800*800
ggplot(sample_table,aes(x=samplename1,weight=CellNumber,fill=cell_type))+  
  geom_bar(position="fill")+ #做堆积柱状图
  scale_fill_npg()+ #配色方案
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
  )+labs(y="Percentage")+RotatedAxis()

#做绝对数量柱状图,输出800*800
bar_stack<-ggplot(sample_table,aes(x=samplename1,weight=CellNumber,fill=cell_type))+  
  geom_bar(position="stack")+ #做绝对数量柱状图
  scale_fill_npg()+ #配色方案
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
  )+labs(y="Percentage")+RotatedAxis()
bar_stack

#设置marker基因
markergenes=c("Epcam","Slc10a2","Krt8","Muc1","Krt18","Cdh1","Krt19","Msln","Upk3b","Slc39a8",
              "Col1a2","Col1a1","Dcn","Col6a1","Col3a1","Col5a1","Mmp2","Cdh5",
              "Pecam1","Flt4","Cldn5","Kdr","Acta2","Tagln","Myl9","Mylk",
              "Cnn1","Myh11","Notch3","Cd52","Tyrobp","Csf1r","Csf3r","S100a8",
              "Guca1a")
markergenes2=c("Epcam","Slc10a2","Krt8","Muc1","Krt18","Cdh1","Krt19","Msln","Upk3b","Slc39a8",
              "Col1a2","Col1a1","Dcn","Col6a1","Col3a1","Col5a1","Mmp2","Cdh5",
              "Pecam1","Kdr","Acta2","Tagln","Myl9","Mylk",
              "Cnn1","Myh11","Notch3","Cd52","Tyrobp")
#按照细胞类型画气泡图，输出1000*400
DotPlot(project,
        features = unique(markergenes),
        group.by = 'Cluster',
        assay="RNA",)+RotatedAxis()
#给cluster按照细胞类型设置顺序
project@meta.data$Cluster<-factor(project@meta.data$Cluster,
                                    levels = c("1","2","8","15","16","22","31",
                                               "32","19","5","6","7","11","13","18",
                                               "29","4","24","12","20","26","3",
                                               "9","10","14","17","21","23","25","27",
                                               "28","30"))
project@meta.data$Cluster<-factor(project@meta.data$Cluster,
                                  levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32"))
#做所有cluster的气泡图，输出1000*1000
DotPlot(project,
        features = as.character(unique(markergenes)),
        group.by = "Cluster",
        assay="RNA")+RotatedAxis()
#按照细胞类型画热图，输出800*800
DoHeatmap(project,features =markergenes,
          group.by = 'cell_type',
          assay = "RNA",
          group.colors = c("#E64B35B2","#4DBBD5B2","#00A087B2","#3C5488B2","#F39B7FB2")
          )+
  scale_fill_gradientn(colors = c("navy","white","firebrick3"))#设置热图颜色 
#按照cluster画热图，输出800*800
DoHeatmap(project,features =markergenes,
          group.by = 'Cluster',
          assay = "RNA")+
  scale_fill_gradientn(colors = c("navy","white","firebrick3"))#设置热图颜色 

#制作单基因umap图，输出1500*3000
FeaturePlot(project, features = markergenes,cols = c('lightgrey', 'blue','seagreen2'), pt.size = 0.5)+  
  theme_dr()+ #小坐标轴
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
#制作单基因umap图，输出800*800
FeaturePlot(project, features = "Tyrobp",cols = c('lightgrey', 'blue','seagreen2'), pt.size = 0.5)+  
  theme_dr()+ #小坐标轴
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
##绘制小提琴图
install.packages("remotes")  
remotes::install_github("lyc-1995/MySeuratWrappers")#通过链接安装包  
library(MySeuratWrappers)  
library(remotes)  
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')#颜色设置  
#以cluster,输出1000*700
VlnPlot(project, features = markergenes,  
        group.by = 'Cluster',
        stacked=T,pt.size=0,
        cols = my36colors,#颜色 
        direction = "horizontal",     #水平作图 
        x.lab = '', y.lab = '')+    #横纵轴不标记任何东西 
  theme(axis.text.x = element_blank(),   
                axis.ticks.x = element_blank())#不显示坐标刻度
#以celltype,输出1000*400
VlnPlot(project, features = markergenes2,  
        group.by = 'cell_type',
        stacked=T,pt.size=0,
        cols = my36colors,#颜色 
        direction = "horizontal",     #水平作图 
        x.lab = '', y.lab = '')+    #横纵轴不标记任何东西 
  theme(axis.text.x = element_blank(),   
        axis.ticks.x = element_blank())#不显示坐标刻度


#提取epithlia cell,fibroblast做亚群分析####
table(Idents(project))
Idents(project)<-project$cell_type
epifibmescell<-subset(project,idents=c("Epithelial cell","Fibroblast cell","Mesothelial cell"))
fibrocell<-subset(project,cell_type=="Fibroblast cell") #提取亚群
epicell2<-subset(project,cell_type=="Epithelial cell") #提取亚群
epifibmescell<-subset(project,cell_type=c("Epithelial cell","Fibroblast cell","Mesothelial cell"))
epicell2 <- FindVariableFeatures(epicell2, selection.method = 'vst', nfeatures = 2000) #寻找高表达变化基因
epicell2<- RunPCA(epicell2, features = VariableFeatures(object = epicell2)) #重新PCA
epicell2 <- JackStraw(epicell2, num.replicate = 100)
epicell2 <- ScoreJackStraw(epicell2, dims = 1:20)
ElbowPlot(epicell2) #做碎石图
epicell2 <- FindNeighbors(epicell2, dims = 1:10)#选择10个维度
epicell2 <- FindClusters(epicell2, resolution = c(seq(.1,2,.2))) #按分辨率从0.1-2作分类，间隔0.2
table(epicell2$seurat_clusters)
head(Idents(epicell2), 5)
colnames((epicell2@meta.data))#看分类后的表头

##clustree看不同分辨率下的分群
install.packages("clustree")
library(clustree)
clustree(epicell2@meta.data, prefix = "RNA_snn_res.")##输出1000*800
epicell2 <- RunUMAP(epicell2, dims = 1:10)#做umap图数据
DimPlot(epicell2,reduction = "umap")#看umap图
P1=DimPlot(epicell2,reduction = "umap",group.by = "RNA_snn_res.0.1",label = TRUE)
P2=DimPlot(epicell2,reduction = "umap",group.by = "RNA_snn_res.0.5",label = TRUE)
P1+P2

#选择分辨率为0.9重新run umap
epicell2 <- FindClusters(epicell2, resolution = 0.9) #按分辨率0.5
colnames((epicell2@meta.data))#看分类后的表头
DimPlot(epicell2, reduction = 'umap',split.by = "samplename1",
        label = TRUE, label.size=8,pt.size = 1)+
  theme_dr()+ #小坐标轴
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
#寻找所有marker基因
epicellmarkers <- FindAllMarkers(object = epicell2, only.pos = TRUE)
head(epicellmarkers,n=5)
#top10markers
epicelltop10<-epicellmarkers %>%
  group_by(cluster) %>%
  top_n(n=10,wt=avg_log2FC)
DoHeatmap(epicell2,features = epicelltop10$gene)

##epi--寻找不同组别间差异基因,vehicle与3周######
Idents(epicell2)<-"samplename1" #将组别设置为idents
vevs3fiffgenes<- FindMarkers(object = epicell2, only.pos = FALSE,assay = 'RNA',
                             ident.1 = "3weeks",ident.2 = "Vehicle",min.pct = 0)#vehicle和3weeks之间差异基因，log2大于0表示在ident1上调，小于0表示在ident1下调
sig_vevs3diffgenes <- subset(vevs3fiffgenes, p_val_adj<0.01)

VlnPlot(epicell2, features = c("Cxcl5","Cxcl12","Aqp1","S100a8","S100a9","Avir","Trpm5","Chil4","Dcn","Fn1","Col1a1","Lcn2","Ccl11","Apoe","Sult1a1","Sult1c2","Il25","Dclk1","Chat","Alox5","Pou2f3"),  
        group.by = 'seurat_clusters',
        stacked=T,pt.size=0,
        cols = my36colors,#颜色 
        direction = "horizontal",     #水平作图 
        x.lab = '', y.lab = '')+    #横纵轴不标记任何东西 
  theme(axis.text.x = element_blank(),   
        axis.ticks.x = element_blank())#不显示坐标刻度

##epi--寻找不同组别间差异基因,vehicle与6周######
Idents(epicell2)<-"samplename1" #将组别设置为idents
vevs6fiffgenes<- FindMarkers(object = epicell2, only.pos = FALSE,assay = 'RNA',
                             ident.1 = "6wkeeks",ident.2 = "Vehicle",min.pct = 0)#vehicle和3weeks之间差异基因，log2大于0表示在ident1上调，小于0表示在ident1下调
sig_vevs6diffgenes <- subset(vevs6fiffgenes, p_val_adj<0.01)


##3周与6周
w3vsw6fiffgenes<- FindMarkers(object = epicell2, only.pos = FALSE,assay = 'RNA',
                             ident.1 = "6wkeeks",ident.2 = "3weeks",min.pct = 0)#6weeks和3weeks之间差异基因，log2大于0表示在ident1上调，小于0表示在ident1下调
sig_w3vsw6diffgenes <- subset(w3vsw6fiffgenes, p_val_adj<0.01)
w3vsw6fiffgenes_up=rownames(sig_w3vsw6diffgenes[sig_w3vsw6diffgenes$avg_log2FC > 0,])#获取6weeks上调基因
w3vsw6fiffgenes_down=rownames(sig_w3vsw6diffgenes[sig_w3vsw6diffgenes$avg_log2FC < 0,])#获取6weeks下调基因


##上调取交集，画图
install.packages("VennDiagram")
library(VennDiagram)
venn.diagram(
  list(week3 = vevs3fiffgenes_up, week6 =vevs6fiffgenes_up), 
  cat.default.pos='text',imagetype = "png", alpha=c(0.8,0.8),height = 4000, width = 4000,
  fill=c("#53A85F","#C1E6F3"), 
  filename = "Venn_up2.png",scaled=T) #画图

venn.diagram(
  list(week3 = vevs3fiffgenes_up, week6 =vevs6fiffgenes_up, w3_vs_w6=w3vsw6fiffgenes_up), 
  cat.default.pos='text',imagetype = "png", alpha=c(0.8,0.8,0.8),height = 4000, width = 4000,
  fill=c("#53A85F","#C1E6F3","#E95C59"), 
  filename = "Venn_up3.png",scaled=T) #画图
up_gene <- intersect(vevs3fiffgenes_up, vevs6fiffgenes_up)
up_gene <- intersect(up_gene,w3vsw6fiffgenes_up)

VlnPlot(epicell2, features = c("Sult1a1","Slc6a6","Slc4a4","Akr1c19","Gas6","Fam213a","Ptprk","Cyp2f2","Gstm1","Gstm2","Egln2","Prss23","Abcb1a","Slc26a3","Aqp1","Slc2a2","Muc3","Cldn4"),  
        group.by = 'samplename1',
        stacked=T,pt.size=0,
        cols = my36colors,#颜色 
        direction = "horizontal",     #水平作图 
        x.lab = '', y.lab = '')+    #横纵轴不标记任何东西 
  theme(axis.text.x = element_blank(),   
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(color='black'),
        axis.line.x = element_line(color='black'))#不显示坐标刻度
up_genego_CC <- enrichGO(gene= up_gene,
                           OrgDb         = 'org.Mm.eg.db',
                           keyType       = 'SYMBOL',
                           ont           = "CC",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.01,
                           qvalueCutoff  = 0.05)
up_genego_MF <- enrichGO(gene     = up_gene,
                           OrgDb         = 'org.Mm.eg.db',
                           keyType       = 'SYMBOL',
                           ont           = "MF",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.01,
                           qvalueCutoff  = 0.05)
up_genego_BP <- enrichGO(gene       = up_gene,
                           OrgDb         = 'org.Mm.eg.db',
                           keyType       = 'SYMBOL',
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.01,
                           qvalueCutoff  = 0.05)   
up_genego_bp<-data.frame(up_genego_BP)
write.csv(up_genego_bp,"up_genego_bp.csv")#输出结果
up_gene_BP <- barplot(up_genego_BP,showCategory = 10) + ggtitle("up_gene for Biological process")
up_gene_CC <- barplot(up_genego_CC,showCategory = 10) + ggtitle("up_gene for Cellular component")
up_gene_MF <- barplot(up_genego_MF,showCategory = 10) + ggtitle("up_gene for Molecular function")
up_gene_BP/up_gene_CC/up_gene_MF

##下调取交集，画图
venn.diagram(
  list(week3 = vevs3fiffgenes_down, week6 =vevs6fiffgenes_down), 
  imagetype = "png", alpha=c(0.8,0.8),height = 4000, width = 4000,
  fill=c("#53A85F","#C1E6F3"), 
  filename = "Venn_down2-2.png",scaled=T) #画图

venn.diagram(
  list(week3 = vevs3fiffgenes_down, week6 =vevs6fiffgenes_down, w3_vs_w6=w3vsw6fiffgenes_down), 
  cat.default.pos='text',imagetype = "png", alpha=c(0.8,0.8,0.8),height = 4000, width = 4000,
  fill=c("#53A85F","#C1E6F3","#E95C59"), 
  filename = "Venn_down3.png",scaled=T) #画图
down_gene <- intersect(vevs3fiffgenes_down, vevs6fiffgenes_down)
down_gene <- intersect(down_gene,w3vsw6fiffgenes_down)
VlnPlot(epicell2, features = down_gene,  
        group.by = 'samplename1',
        stacked=T,pt.size=0,
        cols = my36colors,#颜色 
        direction = "horizontal",     #水平作图 
        x.lab = '', y.lab = '')+    #横纵轴不标记任何东西 
  theme(axis.text.x = element_blank(),   
        axis.ticks.x = element_blank())#不显示坐标刻度
up_genego_CC <- enrichGO(gene= up_gene,
                         OrgDb         = 'org.Mm.eg.db',
                         keyType       = 'SYMBOL',
                         ont           = "CC",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05)
up_genego_MF <- enrichGO(gene     = up_gene,
                         OrgDb         = 'org.Mm.eg.db',
                         keyType       = 'SYMBOL',
                         ont           = "MF",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05)
up_genego_BP <- enrichGO(gene       = up_gene,
                         OrgDb         = 'org.Mm.eg.db',
                         keyType       = 'SYMBOL',
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05)   
up_genego_bp<-data.frame(up_genego_BP)
write.csv(up_genego_bp,"up_genego_bp.csv")#输出结果
up_gene_BP <- barplot(up_genego_BP,showCategory = 10) + ggtitle("up_gene for Biological process")
up_gene_CC <- barplot(up_genego_CC,showCategory = 10) + ggtitle("up_gene for Cellular component")
up_gene_MF <- barplot(up_genego_MF,showCategory = 10) + ggtitle("up_gene for Molecular function")
up_gene_BP/up_gene_CC/up_gene_MF




#############细胞类型注释，SingleR包安装和使用####
if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("SingleR",force = TRUE)  #在官网https://bioconductor.org/packages/release/bioc/html/SingleR.html 找到下载代码
library(SingleR)
BiocManager::install("celldex")
library(celldex)
#下载参考marker数据集
ImmGen.se<-ImmGenData() #(鼠)
Mouse.se<-MouseRNAseqData() #(鼠)

##用markers做epi的小提琴图
epimarkers4 = c("Epcam","Col1a1","Acta2","Muc3","Cxcl5","Lcn2","Dclk1","Trpm5","Mki67")

VlnPlot(epicell2, features = epimarkers4,  
          group.by = 'RNA_snn_res.0.9',
        stacked=T,pt.size=0,
        cols = my36colors,#颜色 
        direction = "horizontal",     #水平作图 
        x.lab = '', y.lab = '')+    #横纵轴不标记任何东西 
  theme(axis.text.x = element_blank(),   
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(color='black'),
        axis.line.x = element_line(color='black'))#不显示坐标刻度
#细胞注释
epictype <- c("0"="Muc3 high",
              "1"="Other cells",
              "2"="Lcn2high/Cxcl15high",
              "3"="Lcn2high/Cxcl15high",
              "4"="Lcn2high/Cxcl15high",
              "5"="Other cells",
              "6"="Other cells",
              "7"="Muc3 high",
              "8"="Muc3 high",
              "9"="Other cells",
              "10"="Muc3 high",
              "11"="Lcn2high/Cxcl15high",
              "12"="Fibroblast-like cell",
              "13"="Fibroblast-like cell",
              "14"="Other cells",
              "15"="Muc3 high",
              "16"="Mki67 high",
              "17"="Tuft cell")
epicell2@meta.data$seurat_clusters<-factor(epicell2@meta.data$seurat_clusters,
                                           levels = c("12","13","2","3","4","11","0","7","8","10","15","5","6","17","16","1","9","14"))

                                           
epicell2[["epicell_type"]]=unname(epictype[epicell2@meta.data$RNA_snn_res.0.9])
epicell2[["epicell_type1"]]=unname(epictype1[epicell2@meta.data$RNA_snn_res.0.9])
epicell2[["epicell_type2"]]=unname(epictype2[epicell2@meta.data$seurat_clusters])
epicell2[["epicell_type3"]]=unname(epictype3[epicell2@meta.data$seurat_clusters])
epicell2[["epicell_type4"]]=unname(epictype4[epicell2@meta.data$seurat_clusters])

epicell2@meta.data$epicell_type1<-factor(epicell2@meta.data$epicell_type1,
                                    levels = c("Stem-like cell","Other cell"))
epicell2@meta.data$epicell_type2<-factor(epicell2@meta.data$epicell_type2,
                                         levels = c("EMT","Other cell"))
epicell2@meta.data$epicell_type3<-factor(epicell2@meta.data$epicell_type3,
                                         levels = c("Muc3+ cell","Other cell"))
epicell2@meta.data$epicell_type4<-factor(epicell2@meta.data$epicell_type4,
                                         levels = c("Fibroblast-like cell","Other cell"))
epicell2@meta.data$epicell_type<-factor(epicell2@meta.data$epicell_type,
                                         levels = c("Fibroblast-like cell","Muc3 high","Lcn2high/Cxcl15high","Mki67 high","Tuft cell","Other cells"))

#按照cluster画热图，输出800*800
DoHeatmap(epicell2,features = epicelltop10$gene,
          group.by = 'epicell_type',
          assay = "RNA")+
  scale_fill_gradientn(colors = c("navy","white","firebrick3")) #设置热图颜色
#输出原始样本的umap图,1000*1000
DimPlot(epicell2, reduction = 'umap',group.by = 'samplename1',
        label.size=8,pt.size = 2)+
  scale_color_npg() +
  theme_dr()+ #小坐标轴
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
#输出像素为2400*1600
DimPlot(epicell2, reduction = 'umap',group.by = 'epicell_type',
        label.size=8,pt.size = 2)+
  scale_color_npg() +
  theme_dr()+ #小坐标轴
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#做柱状图,先建立一个数据框
epic_table<- as.data.frame(table(epicell2@meta.data$samplename1,epicell2@meta.data$epicell_type))
names(epic_table) <- c("samplename1","epicell_type","CellNumber")
epictable<-subset(epic_table,epicell_type!= "Fibroblast cell")#删除Fibroblast cell
#做堆积柱状图,输出800*1000
ggplot(epic_table,aes(x=samplename1,weight=CellNumber,fill=epicell_type))+  
  geom_bar(position="fill")+ #做堆积柱状图
  scale_fill_npg()+ #配色方案
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
  )+labs(y="Percentage")+RotatedAxis()

#细胞类型小提琴图，300*300
VlnPlot(epicell2, features = c("Epcam","Col1a1","Acta2","Fn1","Cxcl5","Olfm4","Hmgcs2","Muc3"),  
        group.by = 'seurat_clusters',
        stacked=T,pt.size=0,
        cols = my36colors,#颜色 
        direction = "horizontal",     #水平作图 
        x.lab = '', y.lab = '')+    #横纵轴不标记任何东西 
  theme(axis.text.x = element_blank(),   
        axis.ticks.x = element_blank())#不显示坐标刻度

#制作单基因umap图，输出800*800
FeaturePlot(epicell2, features = "Acta2",cols = c('lightgrey', 'blue','seagreen2'), pt.size = 0.5)+  
  theme_dr()+ #小坐标轴
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


#################################拟时序轨迹分析##################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("monocle")
library(monocle)
#先把cluster顺序回归默认
epicount<-as(as.matrix(epicell2@assays$RNA@counts),'sparseMatrix') #seurat中提取counts数据
epipd <- new('AnnotatedDataFrame', data = epicell2@meta.data)
fData <- data.frame(gene_short_name = row.names(epicount), row.names = row.names(epicount))
epifd <- new('AnnotatedDataFrame', data = fData)

fibcount<-as(as.matrix(fibrocell@assays$RNA@counts),'sparseMatrix') #seurat中提取counts数据
fibpd <- new('AnnotatedDataFrame', data = fibrocell@meta.data)
fibfData <- data.frame(gene_short_name = row.names(fibcount), row.names = row.names(fibcount))
fibfd <- new('AnnotatedDataFrame', data = fibfData)

efmcount<-as(as.matrix(epifibmescell@assays$RNA@counts),'sparseMatrix') #seurat中提取counts数据
efmpd <- new('AnnotatedDataFrame', data = epifibmescell@meta.data)
efmfData <- data.frame(gene_short_name = row.names(efmcount), row.names = row.names(efmcount))
efmfd <- new('AnnotatedDataFrame', data = efmfData)
#构建S4对象，CellDataSet
epitime <- newCellDataSet(epicount,
                       phenoData = epipd,
                       featureData = epifd,
                       lowerDetectionLimit = 0.5,
                       expressionFamily = negbinomial.size()) 
fibtime <- newCellDataSet(fibcount,
                          phenoData = fibpd,
                          featureData = fibfd,
                          lowerDetectionLimit = 0.5,
                          expressionFamily = negbinomial.size()) 
efmtime <- newCellDataSet(efmcount,
                          phenoData = efmpd,
                          featureData = efmfd,
                          lowerDetectionLimit = 0.5,
                          expressionFamily = negbinomial.size()) 
#提取的是单细胞稀疏矩阵，用negbinomial.size()，如果是UMI的话不要标准化；FPKM值用tobit()；logFPKM值用gaussianff()
##估计size factors 和 dispersions
epitime <- estimateSizeFactors(epitime)
epitime <- estimateDispersions(epitime)

efmtime <- estimateSizeFactors(efmtime)
efmtime <- estimateDispersions(efmtime)

fibtime <- estimateSizeFactors(fibtime)
fibtime <- estimateDispersions(fibtime)

saveRDS(efmtime,file="efmtime.rds")
saveRDS(epitime,file="epitime.rds")
#过滤细胞
epitime <- detectGenes(epitime, min_expr = 2 )
efmtime <- detectGenes(efmtime, min_expr = 2 )
fibtime <- detectGenes(fibtime, min_expr = 2 )
print(head(fData(epitime )))
##过滤低质量细胞
epiexpressed_genes <- row.names(subset(fData(epitime ),
                                    num_cells_expressed >= 10))
length(epiexpressed_genes)
head(pData(epitime))
##根据染毒时间选择差异基因
epidiff_test_res <- differentialGeneTest(epitime[epiexpressed_genes,],
                                         fullModelFormulaStr = "~Sample")



####提取免疫细胞亚群分析####
immunecell<-subset(project,cell_type=="Immune cell") #提取亚群
immunecell <- FindVariableFeatures(immunecell, selection.method = 'vst', nfeatures = 2000) #寻找高表达变化基因
immunecell<- RunPCA(immunecell, features = VariableFeatures(object = immunecell)) #重新PCA
immunecell <- JackStraw(immunecell, num.replicate = 100)
immunecell <- ScoreJackStraw(immunecell, dims = 1:20)
ElbowPlot(immunecell) #做碎石图
immunecell <- FindNeighbors(immunecell, dims = 1:10)#选择10个维度
immunecell <- FindClusters(immunecell, resolution = c(seq(.1,2,.2))) #按分辨率从0.1-2作分类，间隔0.2
table(immunecell$seurat_clusters)
head(Idents(immunecell), 5)
colnames((immunecell@meta.data))#看分类后的表头
##clustree看不同分辨率下的分群
install.packages("clustree")
library(clustree)
clustree(immunecell@meta.data, prefix = "RNA_snn_res.")##输出1000*800
immunecell <- RunUMAP(immunecell, dims = 1:10)#做umap图数据
DimPlot(immunecell,reduction = "umap")#看umap图
P1=DimPlot(immunecell,reduction = "umap",group.by = "RNA_snn_res.0.1",label = TRUE)
P2=DimPlot(immunecell,reduction = "umap",group.by = "RNA_snn_res.0.5",label = TRUE)
P1+P2
#选择分辨率为0.5重新run umap
immunecell <- FindClusters(immunecell, resolution = 0.5) #按分辨率0.5
colnames((immunecell@meta.data))#看分类后的表头
DimPlot(immunecell, reduction = 'umap',split.by = "samplename1",
        label = TRUE, label.size=8,pt.size = 1)+
  theme_dr()+ #小坐标轴
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
#寻找所有marker基因
immunecellmarkers <- FindAllMarkers(object = immunecell, only.pos = TRUE)
head(immunecellmarkers,n=5)
#top10markers
library(dplyr)
immunecelltop10<-immunecellmarkers %>%
  group_by(cluster) %>%
  top_n(n=10,wt=avg_log2FC)
DoHeatmap(immunecell,features = immunecelltop10$gene)
immunemarkers1=c("Csf1r","Lyz1","Cd300e","Mafb","Batf3","Krt79","Msr1",
                 "Siglech","Ccr9","Bst2","Pacsin1","Tcf4",
                 "S100a9","S100a8","Mmp9","Csf3r","Mmp8","Il1rn","Cxcr2",
                 "Il6","Gata2","Cpa3","Ms4a2","Fcer1a",
                 "Cd79a","Fcmr","Cd79b","Cd19","Fcer2a","Pax5","Cd22",
                 "Cd3d","Cd3e","Cd3g","Cd28",
                 "Gzma","Klrb1c","Ncr1","Klra4","Eomes","Gzmb","Fasl","Klrk1")
immunemarkers2=c("Adgre1","Itgam","Fcgr1","Cd68","C1qa","C1qb","Ms4a7","C1qc",
                 "Lyz2","Ly6c2","Plac8","Ccr2",
                 "Mrc1","Cd163","Folr2","Arg1","Chil3",
                 "Epcam","H2-Ab1","H2-Eb1","Cd207","Cd24a")
immunemarkers3=c("Mgl2","Retnla","Arg1","Chil3","Fizz1","Ym1","Ccl24","Irf4")
immunemarkers2=c("Vcan","Ly6c1","Ly6c2","Cd300e","Itgal","Ace","Fcgr4","","","","","","","","","","","","","",)
immunemarkers4=c("Cd79b","Cd79a",
                 "Bst2","Tcf4",
                 "Adgre1","Cd68","C1qa","C1qb","C1qc",
                 "S100a9","S100a8","Guca1a","Csf3r",
                 "Cd3d","Cd3e","Cd3g",
                 "Ccr2")
immunecell@meta.data$Cluster<-factor(immunecell@meta.data$Cluster,
                                  levels = c("28","14","30","10","21","25","3","27","9","17","23"))
VlnPlot(immunecell, features = immunemarkers4,  
        group.by = 'Cluster',
        stacked=T,pt.size=0,
        cols = my36colors,#颜色 
        direction = "horizontal",     #水平作图 
        x.lab = '', y.lab = '')+    #横纵轴不标记任何东西 
  theme(axis.text.x = element_blank(),   
        axis.ticks.x = element_blank())#不显示坐标刻度
####使用singler注释
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SingleR",force = TRUE)  #在官网https://bioconductor.org/packages/release/bioc/html/SingleR.html 找到下载代码
library(SingleR)
BiocManager::install("celldex")
library(celldex)
#下载参考marker数据集
ImmGen.se<-ImmGenData() #(鼠)
Mouse.se<-MouseRNAseqData() #(鼠)

immenecell_for_SingleR <- GetAssayData(immunecell, slot="data")
immunecell.ImmGen <- SingleR(test = immenecell_for_SingleR, ref = ImmGen.se, labels = ImmGen.se$label.main)
immunecell.ImmGen_fine <- SingleR(test = immenecell_for_SingleR, ref = ImmGen.se, labels = ImmGen.se$label.fine)
immunecell.ImmGen
table(immunecell.ImmGen$labels, immunecell$seurat_clusters)
table(immunecell.ImmGen_fine$labels, immunecell$seurat_clusters)
immunecell@meta.data$labelmain <-immunecell.ImmGen$labels
immunecell@meta.data$labelfine <-immunecell.ImmGen_fine$labels
##免疫细胞注释
immunetype <- c("0"="Neutrophils",
               "1"="Neutrophils",
               "2"="NKT",
               "3"="NKT",
               "4"="DC",
               "5"="M2+Macrophages",
               "6"="Macrophages",
               "7"="Macrophages",
               "8"="innate lymphoid cells",
               "9"="Macrophages",
               "10"="DC",
               "11"="Other cells",
               "12"="DC",
               "13"="B cells")
immunecell[["immune_type"]]=unname(immunetype[immunecell@meta.data$seurat_clusters])

#immunetype <- c("28"= "B cells", 
                "14"= "DC",
                "30"= "DC",
                "10"= "Macrophages",
                "21"= "Macrophages",
                "25"= "Macrophages",
                "3"= "Neutrophils",
                "27"= "Neutrophils",
                "9"= "T", 
                "17"= "T", 
                "23"= "Monocyt" )
immunecell[["immune_type"]]=unname(immunetype[immunecell@meta.data$Cluster])#
#输出原始样本的umap图,1000*1000
DimPlot(immunecell, reduction = 'umap',group.by = 'immune_type',
        label = TRUE,label.size=6,pt.size = 2)+
  scale_color_npg() +
  theme_dr()+ #小坐标轴
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



#####成纤维细胞亚群分析#####
fibrocell<-subset(project,cell_type=="Fibroblast cell") #提取亚群
fibrocell <- FindVariableFeatures(fibrocell, selection.method = 'vst', nfeatures = 2000) #寻找高表达变化基因
fibrocell<- RunPCA(fibrocell, features = VariableFeatures(object = fibrocell)) #重新PCA
fibrocell <- JackStraw(fibrocell, num.replicate = 100)
fibrocell <- ScoreJackStraw(fibrocell, dims = 1:20)
ElbowPlot(fibrocell) #做碎石图
fibrocell <- FindNeighbors(fibrocell, dims = 1:7)#选择7个维度
fibrocell <- FindClusters(fibrocell, resolution = c(seq(.1,2,.2))) #按分辨率从0.1-2作分类，间隔0.2
table(fibrocell$seurat_clusters)
head(Idents(fibrocell), 5)
colnames((fibrocell@meta.data))#看分类后的表头
##clustree看不同分辨率下的分群
install.packages("clustree")
library(clustree)
clustree(fibrocell@meta.data, prefix = "RNA_snn_res.")##输出1000*800
fibrocell <- RunUMAP(fibrocell, dims = 1:10)#做umap图数据
DimPlot(fibrocell,reduction = "umap")#看umap图
P1=DimPlot(fibrocell,reduction = "umap",group.by = "RNA_snn_res.0.5",label = TRUE)
P2=DimPlot(fibrocell,reduction = "umap",group.by = "RNA_snn_res.0.7",label = TRUE)
P1+P2
#选择分辨率为0.3重新run umap
fibrocell <- FindClusters(fibrocell, resolution = 0.3) #按分辨率0.3
colnames((fibrocell@meta.data))#看分类后的表头
DimPlot(fibrocell, reduction = 'umap',group.by = "samplename1",
        label = FALSE, label.size=8,pt.size = 1)+
  theme_dr()+ #小坐标轴
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#寻找所有marker基因
fibrocellmarkers <- FindAllMarkers(object = fibrocell, only.pos = TRUE)
head(fibrocellmarkers,n=5)
#top10markers
library(dplyr)
fibrocelltop10<-fibrocellmarkers %>%
  group_by(cluster) %>%
  top_n(n=10,wt=avg_log2FC)
DoHeatmap(fibrocell,features = fibrocelltop10$gene)



fibromarkers=c("Col1a1","Col3a1","Dcn","Col12a1","Acta2","Ccl2","Ccl7","Il6","Cxcl12")
VlnPlot(fibrocell, features = fibromarkers,  
        group.by = 'seurat_clusters',
        stacked=T,pt.size=0,
        cols = my36colors,#颜色 
        direction = "horizontal",     #水平作图 
        x.lab = '', y.lab = '')+    #横纵轴不标记任何东西 
  theme(axis.text.x = element_blank(),   
        axis.ticks.x = element_blank())#不显示坐标刻度
DotPlot(fibrocell,
        features = unique(fibromarkers),
        group.by = 'seurat_clusters',
        assay="RNA",)+RotatedAxis()
#单基因umap图800*800
FeaturePlot(fibrocell, features = "Acta2",cols = c('lightgrey', 'blue','seagreen2'), pt.size = 0.5)+  
  theme_dr()+ #小坐标轴
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

##亚群细胞注释
fibrotype <- c("0"="myofibroblast-Cxcl12high",
                "1"="myofibroblast-Cxcl12high",
                "2"="fibroblasts",
                "3"="myofibroblast-Cxcl12high",
                "4"="fibroblasts-Col12a1high",
                "5"="myofibroblast-Ccl2+Ccl7+",
                "6"="inflammation fibrobalsts",
                "7"="fibroblasts")
fibrocell[["fibrotype"]]=unname(fibrotype[fibrocell@meta.data$seurat_clusters])

#输出原始样本的umap图,1000*1000
DimPlot(fibrocell, reduction = 'umap',group.by = 'fibrotype',split.by = "samplename1",
        label = FALSE,label.size=6,pt.size = 2)+
  scale_color_npg() +
  theme_dr()+ #小坐标轴
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

fibro_table2<- as.data.frame(table(fibrocell@meta.data$samplename1,fibrocell@meta.data$fibrotype))
names(fibro_table2) <- c("samplename1","fibrotype","CellNumber")
ggplot(fibro_table2,aes(x=samplename1,weight=CellNumber,fill=fibrotype))+  
  geom_bar(position="fill")+ #做堆积柱状图
  scale_fill_npg()+ #配色方案
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
  )+labs(y="Percentage")+RotatedAxis()

##亚群差异基因通路富集KEGG GO#####
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
install.packages("corrplot")
install.packages("tidyverse")
BiocManager::install("clusterProfiler")
library(corrplot)
library(tidyverse)
library(RColorBrewer)
library(clusterProfiler)
myo013markers<- FindMarkers(object = fibrocell, only.pos = FALSE,
                              ident.1 = c(0,1,3),min.pct = 0.25)#013和其他之间差异基因
myo013gene_up=rownames(myo013markers[myo013markers$avg_log2FC > 0,])#获取上调基因
myo013gene_down=rownames(myo013markers[myo013markers$avg_log2FC < 0,])#获取下调基因
## 把SYMBOL改为ENTREZID
myo013gene_upENTREZID=as.character(na.omit(AnnotationDbi::select(org.Mm.eg.db,
                                                   keys = myo013gene_up,
                                                   columns = 'ENTREZID',
                                                   keytype = 'SYMBOL')[,2]))
myo013gene_downENTREZID=as.character(na.omit(AnnotationDbi::select(org.Mm.eg.db,
                                                     keys = myo013gene_down,
                                                     columns = 'ENTREZID',
                                                     keytype = 'SYMBOL')[,2]))
##KEGG通路
myo013gene_upENTREZID <- unique(myo013gene_upENTREZID)
myokk.up <- enrichKEGG(gene = myo013gene_upENTREZID,
                    organism = "mmu",
                    pvalueCutoff = 0.9,
                    qvalueCutoff = 0.9)
myo013gene_downENTREZID <- unique(myo013gene_downENTREZID)
myokk.down <- enrichKEGG(gene = myo013gene_downENTREZID,
                       organism = "mmu",
                       pvalueCutoff = 0.9,
                       qvalueCutoff = 0.9)
dotplot(myokk.up)
dotplot(myokk.down )
##GO通路
myogo_CC <- enrichGO(gene= myo013gene_up,
                       OrgDb         = 'org.Mm.eg.db',
                       keyType       = 'SYMBOL',
                       ont           = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05)
myogo_MF <- enrichGO(gene     = myo013gene_up,
                       OrgDb         = 'org.Mm.eg.db',
                       keyType       = 'SYMBOL',
                       ont           = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05)
myogo_BP <- enrichGO(gene       = myo013gene_up,
                       OrgDb         = 'org.Mm.eg.db',
                       keyType       = 'SYMBOL',
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05)    
myogo_CC@result$Description <- substring(myogo_CC@result$Description,1,70)
myogo_MF@result$Description <- substring(myogo_MF@result$Description,1,70)
myogo_BP@result$Description <- substring(myogo_BP@result$Description,1,70)
myop_BP <- barplot(myogo_BP,showCategory = 10) + ggtitle("myobarplot for Biological process")
myop_CC <- barplot(myogo_CC,showCategory = 10) + ggtitle("myobarplot for Cellular component")
myop_MF <- barplot(myogo_MF,showCategory = 10) + ggtitle("myobarplot for Molecular function")
myoplotc <- myop_BP/myop_CC/myop_MF
ggsave('myoenrichGO.png',myoplotc, width = 12,height = 20)

myogo_CCdown <- enrichGO(gene= myo013gene_down,
                     OrgDb         = 'org.Mm.eg.db',
                     keyType       = 'SYMBOL',
                     ont           = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05)
myogo_MFdown <- enrichGO(gene     = myo013gene_down,
                     OrgDb         = 'org.Mm.eg.db',
                     keyType       = 'SYMBOL',
                     ont           = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05)
myogo_BPdown <- enrichGO(gene       = myo013gene_down,
                     OrgDb         = 'org.Mm.eg.db',
                     keyType       = 'SYMBOL',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05)    
myogo_CCdown@result$Description <- substring(myogo_CCdown@result$Description,1,70)
myogo_MFdown@result$Description <- substring(myogo_MFdown@result$Description,1,70)
myogo_BPdown@result$Description <- substring(myogo_BPdown@result$Description,1,70)
myop_BPdown <- barplot(myogo_BPdown,showCategory = 10) + ggtitle("myodownbarplot for Biological process")
myop_CCdown <- barplot(myogo_CCdown,showCategory = 10) + ggtitle("myodownbarplot for Cellular component")
myop_MFdown <- barplot(myogo_MFdown,showCategory = 10) + ggtitle("myodownbarplot for Molecular function")
myoplotcdown <- myop_BPdown/myop_CCdown/myop_MFdown
myoplotcdown 

##细胞通讯分析#####
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
suppressMessages(if(!require(CellChat))devtools::install_github("sqjin/CellChat"))##发现ComplexHeatmap版本不对无法安装，所以先用bioconduct安装ComplexHeatmap
library(patchwork)
library(ggalluvial)
library(ggplot2)
library(igraph)
library(dplyr)
library(CellChat)
table(project@meta.data$integrated_snn_res.0.9,project@meta.data$cell_type) #检查标签对不对
#安装不同版本的包 remotes::install_version
#对免疫细胞进行细分，根据前期分类和singleR一起参考
project_for_SingleR <- GetAssayData(project, slot="data")
project.ImmGen <- SingleR(test = project_for_SingleR, ref = ImmGen.se, labels = ImmGen.se$label.main)
table(project.ImmGen$labels, project$Cluster)
clustercelltype3 <- c("1"= "Epithelial cell", 
                      "2"= "Epithelial cell", 
                      "3"= "Neutrophils", 
                      "4"= "Endothelial cell", 
                      "5"= "Fibroblast cell",
                      "6"= "Fibroblast cell", 
                      "7"= "Fibroblast cell", 
                      "8"= "Epithelial cell",
                      "9"= "T",  
                      "10"= "Macrophages",  
                      "11"= "Fibroblast cell",
                      "12"= "Smooth muscle cell",
                      "13"= "Fibroblast cell",
                      "14"= "DC",  
                      "15"= "Epithelial cell",
                      "16"= "Epithelial cell",
                      "17"= "T",  
                      "18"= "Fibroblast cell",
                      "19"= "Mesothelial cell",
                      "20"= "Smooth muscle cell",
                      "21"= "Macrophages",
                      "22"= "Epithelial cell",
                      "23"= "Monocyt",
                      "24"= "Endothelial cell",
                      "25"= "Macrophages",
                      "26"= "Smooth muscle cell",
                      "27"= "Neutrophils",
                      "28"= "B cells", 
                      "29"= "Fibroblast cell",
                      "30"= "DC",
                      "31"= "Epithelial cell",
                      "32"= "Epithelial cell")
project[["cell_type3"]]=unname(clustercelltype3[project@meta.data$seurat_clusters])
#从seurat中提取用于cellchat的数据,也可以使用data.input2<-project[["RNA"]]@data
data.input <- GetAssayData(project, assay = "RNA", slot = "data")
meta_data <- project@meta.data #得到含有细胞名、样本名的meta data 数据
cell_Vehicle <- rownames(meta_data)[meta_data$samplename1 == "Vehicle"]
cell_3weeks <- rownames(meta_data)[meta_data$samplename1 == "3weeks"]
cell_6weeks <- rownames(meta_data)[meta_data$samplename1 == "6wkeeks"]#三组细胞的标签分别提取出来
#按细胞标签将三组input数据分别提取出来
data.input_Vehicle <- data.input[, cell_Vehicle]
data.input_3weeks <- data.input[, cell_3weeks]
data.input_6weeks <- data.input[, cell_6weeks]
#三组细胞的meta信息提取出来
meta_data_Vehicle = meta_data[cell_Vehicle, ]
meta_data_3weeks  = meta_data[cell_3weeks, ]
meta_data_6weeks = meta_data[cell_6weeks, ]
#创建每组的cellchat对象,group.by指定通讯间的对象，用meta中的注释作为分组依据
cellchat_Vehicle <- createCellChat(object = data.input_Vehicle, meta_data_Vehicle, group.by = "cell_type3")
cellchat_3weeks <- createCellChat(object = data.input_3weeks, meta_data_3weeks, group.by = "cell_type3")
cellchat_6weeks <- createCellChat(object = data.input_6weeks, meta_data_6weeks, group.by = "cell_type3")
#载入数据库，use CellChatDB.human if running on human data
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)#看细胞通讯数据库的类型
glimpse(CellChatDB$interaction)#dplyr包的函数展示互作记录
CellChatDB.secreted <- subsetDB(CellChatDB, search = "Secreted Signaling")#选择Secreted Signaling做后续分析
#CellChatDB.use <- CellChatDB#或者使用默认
#将数据库内容载入cellchat对象中
cellchat_Vehicle@DB <- CellChatDB.use
cellchat_3weeks@DB <- CellChatDB.use
cellchat_6weeks@DB <- CellChatDB.use
#表达量预处理
cellchat_Vehicle <- subsetData(cellchat_Vehicle,features = NULL)#取出表达数据
cellchat_3weeks <- subsetData(cellchat_3weeks,features = NULL)#取出表达数据
cellchat_6weeks <- subsetData(cellchat_6weeks,features = NULL)#取出表达数据
cellchat_Vehicle <- identifyOverExpressedGenes(cellchat_Vehicle)#寻找高表达的基因#
cellchat_3weeks  <- identifyOverExpressedGenes(cellchat_3weeks )#寻找高表达的基因#
cellchat_6weeks <- identifyOverExpressedGenes(cellchat_6weeks)#寻找高表达的基因#
###FindVariableFeatures()
cellchat_Vehicle <- identifyOverExpressedInteractions(cellchat_Vehicle)#寻找高表达的通路
cellchat_Vehicle <- projectData(cellchat_Vehicle, PPI.mouse)#投影到PPI
cellchat_Vehicle <- computeCommunProb(cellchat_Vehicle)#默认计算方式为type = "truncatedMean",#默认cutoff的值为20%，即表达比例在25%以下的基因会被认为是0， trim = 0.1可以调整比例阈值
cellchat_Vehicle <- filterCommunication(cellchat_Vehicle, min.cells = 10)#去掉通讯数量很少的细胞
##提取通路信息与配受体信息
df.net_Vehicle <- subsetCommunication(cellchat_Vehicle)#将细胞通讯预测结果以数据框的形式取出
df.net_3weeks <- subsetCommunication(cellchat_3weeks)
df.net_6weeks <- subsetCommunication(cellchat_6weeks)
#DT::datatable(df.net_Vehicle )
write.csv(df.net_Vehicle,'01.df.net_Vehicle.csv')
write.csv(df.net_3weeks,'01.df.net_3weeks.csv')
write.csv(df.net_6weeks,'01.df.net_6weeks.csv')
##单提取通路信息
df.netP_Vehicle <- subsetCommunication(cellchat_Vehicle,slot.name = "netP")
df.netP_3weeks <- subsetCommunication(cellchat_3weeks,slot.name = "netP")
df.netP_6weeks <- subsetCommunication(cellchat_6weeks,slot.name = "netP")
write.csv(df.netP_Vehicle,'01.df.netP_Vehicle.csv')
write.csv(df.netP_3weeks,'01.df.netP_3weeks.csv')
write.csv(df.netP_6weeks,'01.df.netP_6weeks.csv')

##这种方式只取通路，数据结构更简单
#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))
#source可以以细胞类型的名称定义，也可以按照细胞名称中的顺序以数值向量直接取
#指定输入与输出的细胞集群
#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))#指定通路提取
cellchat_Vehicle<- computeCommunProbPathway(cellchat_Vehicle)
#每对配受体的预测结果存在net中，每条通路的预测结果存在netp中
cellchat_Vehicle <- aggregateNet(cellchat_Vehicle)
cellchat_Vehicle <- netAnalysis_computeCentrality(cellchat_Vehicle, slot.name = "netP")#确定signaling角色(例如，主要的发送者，接收者)以及主要的贡献singnaling
#计算联路数与通讯概率，可用sources.use and targets.use指定来源与去向

###FindVariableFeatures()-3weeks
cellchat_3weeks <- identifyOverExpressedInteractions(cellchat_3weeks)#寻找高表达的通路
cellchat_3weeks <- projectData(cellchat_3weeks, PPI.mouse)#投影到PPI
cellchat_3weeks <- computeCommunProb(cellchat_3weeks)#默认计算方式为type = "truncatedMean",#默认cutoff的值为20%，即表达比例在25%以下的基因会被认为是0， trim = 0.1可以调整比例阈值
cellchat_3weeks <- filterCommunication(cellchat_3weeks, min.cells = 10)#去掉通讯数量很少的细胞
cellchat_3weeks<- computeCommunProbPathway(cellchat_3weeks)
cellchat_3weeks <- aggregateNet(cellchat_3weeks)
cellchat_3weeks <- netAnalysis_computeCentrality(cellchat_3weeks, slot.name = "netP")

###FindVariableFeatures()-6weeks
cellchat_6weeks <- identifyOverExpressedInteractions(cellchat_6weeks)#寻找高表达的通路
cellchat_6weeks <- projectData(cellchat_6weeks, PPI.mouse)#投影到PPI
cellchat_6weeks <- computeCommunProb(cellchat_6weeks)#默认计算方式为type = "truncatedMean",#默认cutoff的值为20%，即表达比例在25%以下的基因会被认为是0， trim = 0.1可以调整比例阈值
cellchat_6weeks <- filterCommunication(cellchat_6weeks, min.cells = 10)#去掉通讯数量很少的细胞
cellchat_6weeks<- computeCommunProbPathway(cellchat_6weeks)
cellchat_6weeks <- aggregateNet(cellchat_6weeks)
cellchat_6weeks <- netAnalysis_computeCentrality(cellchat_6weeks, slot.name = "netP")

#将三组数据合并后再做图
obejectlist<-list(Vehicle=cellchat_Vehicle, weeks3=cellchat_3weeks,weeks6=cellchat_6weeks)
cellchatmerge <- mergeCellChat(obejectlist, add.names = names(obejectlist))
#简单看看不同组别间互作数量和权重是否一样
gg1 <- compareInteractions(cellchatmerge, show.legend = F, group = c(1,2,3))
gg2 <- compareInteractions(cellchatmerge, show.legend = F, group = c(1,2,3), measure = "weight")
gg1 + gg2
#细胞通路在三组间的富集程度
gg1<-rankNet(cellchatmerge,  comparison = c(1, 2),mode = "comparison",stacked = T,do.stat = TRUE)
gg2<-rankNet(cellchatmerge,comparison = c(1, 2),mode = "comparison",stacked = F,do.stat = TRUE)
gg1 + gg2
##根据使用netVisual_circle显示所有celltype之间的通讯次数
weight.max <- getMaxWeight(obejectlist, attribute = c("idents","count"))
par(mfrow = c(1,3), xpd=TRUE)
for (i in 1:length(obejectlist)) {
  netVisual_circle(obejectlist[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 10, title.name = paste0("Number of interactions - ", names(obejectlist)[i]))
}
#根据使用netVisual_circle显示所有celltype之间的总通讯强度
weight.max <- getMaxWeight(obejectlist, attribute = c("idents","count"))
par(mfrow = c(1,3), xpd=TRUE)
for (i in 1:length(obejectlist)) {
  netVisual_circle(obejectlist[[i]]@net$weight, weight.scale = T, label.edge= F, edge.width.max = 10, title.name = paste0("Interaction weights/strength - ", names(obejectlist)[i]))
}

#根据使用netVisual_circle显示两组之间的总通讯数目和强度差别，红色表示在第二组上调，蓝色表示在第二组下调
par(mfrow = c(1,4), xpd=TRUE)
netVisual_diffInteraction(cellchatmerge,comparison = c(1, 2), edge.width.max = 4,weight.scale = T,title.name = paste0("Number of interactions - vehicle+3weeks "))
netVisual_diffInteraction(cellchatmerge,comparison = c(1, 2), edge.width.max = 4,weight.scale = T, measure = "weight",title.name = paste0("strength of interactions -vehicle+3weeks "))
netVisual_diffInteraction(cellchatmerge,comparison = c(1, 3), edge.width.max = 4,weight.scale = T,title.name = paste0("Number of interactions - vehicle+6weeks "))
netVisual_diffInteraction(cellchatmerge,comparison = c(1, 3), edge.width.max = 4,weight.scale = T, measure = "weight",title.name = paste0("strength of interactions -vehicle+6weeks "))
#0-3week单个细胞类型Epithelial cell的netVisual_circle显示两组之间的总通讯数目和强度差别
par(mfrow = c(1,4), xpd=TRUE)
netVisual_diffInteraction(cellchatmerge,comparison = c(1, 2), sources.use = c("Epithelial cell"),weight.scale = T,title.name = paste0("Number of interactions - vehicle+3weeks "))
netVisual_diffInteraction(cellchatmerge,comparison = c(1, 2), sources.use = c("Epithelial cell"),weight.scale = T, measure = "weight",title.name = paste0("strength of interactions -vehicle+3weeks "))
netVisual_diffInteraction(cellchatmerge,comparison = c(1, 2), targets.use = c("Epithelial cell"),weight.scale = T,title.name = paste0("Number of interactions - vehicle+3weeks "))
netVisual_diffInteraction(cellchatmerge,comparison = c(1, 2),targets.use = c("Epithelial cell"), weight.scale = T, measure = "weight",title.name = paste0("strength of interactions -vehicle+3weeks "))
#0-6week单个细胞类型Epithelial cell的netVisual_circle显示两组之间的总通讯数目和强度差别
par(mfrow = c(1,4), xpd=TRUE)
netVisual_diffInteraction(cellchatmerge,comparison = c(1, 3), sources.use = c("Epithelial cell"),weight.scale = T,title.name = paste0("Number of interactions - vehicle+6weeks "))
netVisual_diffInteraction(cellchatmerge,comparison = c(1, 3), sources.use = c("Epithelial cell"),weight.scale = T, measure = "weight",title.name = paste0("strength of interactions -vehicle+6weeks "))
netVisual_diffInteraction(cellchatmerge,comparison = c(1, 3), targets.use = c("Epithelial cell"),weight.scale = T,title.name = paste0("Number of interactions - vehicle+6weeks "))
netVisual_diffInteraction(cellchatmerge,comparison = c(1, 3),targets.use = c("Epithelial cell"), weight.scale = T, measure = "weight",title.name = paste0("strength of interactions -vehicle+6weeks "))
#0-3week单个细胞类型Fibroblast cell的netVisual_circle显示两组之间的总通讯数目和强度差别
par(mfrow = c(1,4), xpd=TRUE)
netVisual_diffInteraction(cellchatmerge,comparison = c(1, 2), sources.use = c("Fibroblast cell"),weight.scale = T,title.name = paste0("Number of interactions - vehicle+3weeks "))
netVisual_diffInteraction(cellchatmerge,comparison = c(1, 2), sources.use = c("Fibroblast cell"),weight.scale = T, measure = "weight",title.name = paste0("strength of interactions -vehicle+3weeks "))
netVisual_diffInteraction(cellchatmerge,comparison = c(1, 2), targets.use = c("Fibroblast cell"),weight.scale = T,title.name = paste0("Number of interactions - vehicle+3weeks "))
netVisual_diffInteraction(cellchatmerge,comparison = c(1, 2),targets.use = c("Fibroblast cell"), weight.scale = T, measure = "weight",title.name = paste0("strength of interactions -vehicle+3weeks "))
#0-6week单个细胞类型Fibroblast cell的netVisual_circle显示两组之间的总通讯数目和强度差别
par(mfrow = c(1,4), xpd=TRUE)
netVisual_diffInteraction(cellchatmerge,comparison = c(1, 3), sources.use = c("Fibroblast cell"),weight.scale = T,title.name = paste0("Number of interactions - vehicle+6weeks "))
netVisual_diffInteraction(cellchatmerge,comparison = c(1, 3), sources.use = c("Fibroblast cell"),weight.scale = T, measure = "weight",title.name = paste0("strength of interactions -vehicle+6weeks "))
netVisual_diffInteraction(cellchatmerge,comparison = c(1, 3), targets.use = c("Fibroblast cell"),weight.scale = T,title.name = paste0("Number of interactions - vehicle+6weeks "))
netVisual_diffInteraction(cellchatmerge,comparison = c(1, 3),targets.use = c("Fibroblast cell"), weight.scale = T, measure = "weight",title.name = paste0("strength of interactions -vehicle+6weeks "))
#0-3week单个细胞类型Macrophages的netVisual_circle显示两组之间的总通讯数目和强度差别
par(mfrow = c(1,4), xpd=TRUE)
netVisual_diffInteraction(cellchatmerge,comparison = c(1, 2), sources.use = c("Macrophages"),weight.scale = T,title.name = paste0("Number of interactions - vehicle+3weeks "))
netVisual_diffInteraction(cellchatmerge,comparison = c(1, 2), sources.use = c("Macrophages"),weight.scale = T, measure = "weight",title.name = paste0("strength of interactions -vehicle+3weeks "))
netVisual_diffInteraction(cellchatmerge,comparison = c(1, 2), targets.use = c("Macrophages"),weight.scale = T,title.name = paste0("Number of interactions - vehicle+3weeks "))
netVisual_diffInteraction(cellchatmerge,comparison = c(1, 2),targets.use = c("Macrophages"), weight.scale = T, measure = "weight",title.name = paste0("strength of interactions -vehicle+3weeks "))
#0-6week单个细胞类型Macrophages的netVisual_circle显示两组之间的总通讯数目和强度差别
par(mfrow = c(1,4), xpd=TRUE)
netVisual_diffInteraction(cellchatmerge,comparison = c(1, 3), sources.use = c("Macrophages"),weight.scale = T,title.name = paste0("Number of interactions - vehicle+6weeks "))
netVisual_diffInteraction(cellchatmerge,comparison = c(1, 3), sources.use = c("Macrophages"),weight.scale = T, measure = "weight",title.name = paste0("strength of interactions -vehicle+6weeks "))
netVisual_diffInteraction(cellchatmerge,comparison = c(1, 3), targets.use = c("Macrophages"),weight.scale = T,title.name = paste0("Number of interactions - vehicle+6weeks "))
netVisual_diffInteraction(cellchatmerge,comparison = c(1, 3),targets.use = c("Macrophages"), weight.scale = T, measure = "weight",title.name = paste0("strength of interactions -vehicle+6weeks "))

#使用热图展示不同组别间差别
names(obejectlist) #看出来vehicle是第一组，weeks3是第二组，weeks6是第三组
gg1 <- netVisual_heatmap(cellchatmerge,comparison = c(1, 2),title.name = "Number-vehicle/3weeks ") #第1，2组的通讯数目差别,红色表示在第二组上调，蓝色表示在第二组下调
gg2 <- netVisual_heatmap(cellchatmerge,comparison = c(1, 2), measure = "weight",title.name = "weight-vehicle/3weeks ")#第1，2组的通讯权重差别,红色表示在第二组上调，蓝色表示在第二组下调
gg3 <- netVisual_heatmap(cellchatmerge,comparison = c(1, 3),title.name = "Number-vehicle/6weeks ") #第1，3组的通讯数目差别,红色表示在第3组上调，蓝色表示在第3组下调
gg4 <- netVisual_heatmap(cellchatmerge,comparison = c(1, 3), measure = "weight",title.name = "weight-vehicle/6weeks ")#第1，3组的通讯权重差别,红色表示在第3组上调，蓝色表示在第3组下调
gg1+gg2 +gg3+gg4

#在2D空间中可视化主要的发送者(源)和接收者(目标)。
num.link <- sapply(obejectlist, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(obejectlist)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(obejectlist[[i]], title = names(obejectlist)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)

#比较与每个细胞群相关的传出（或传入）信号
#outgoing
i = 1
pathway.union1 <- union(obejectlist[[i]]@netP$pathways, obejectlist[[i+1]]@netP$pathways)
pathway.union <- union( pathway.union1,obejectlist[[i+2]]@netP$pathways)#合并三个数据库，只能两两合并
ht1 = netAnalysis_signalingRole_heatmap(obejectlist[[i]], pattern = "outgoing", signaling = pathway.union, title = names(obejectlist)[i], width = 5, height =25)
ht2 = netAnalysis_signalingRole_heatmap(obejectlist[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(obejectlist)[i+1],width = 5, height = 25)
ht3 = netAnalysis_signalingRole_heatmap(obejectlist[[i+2]], pattern = "outgoing", signaling = pathway.union, title = names(obejectlist)[i+2],width = 5, height = 25)
draw(ht1 + ht2+ht3, ht_gap = unit(0.5, "cm"))
#incoming
ht4 = netAnalysis_signalingRole_heatmap(obejectlist[[i]], pattern = "incoming", signaling = pathway.union, title = names(obejectlist)[i], width = 5, height = 25)
ht5 = netAnalysis_signalingRole_heatmap(obejectlist[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(obejectlist)[i+1],width = 5, height = 25)
ht6 = netAnalysis_signalingRole_heatmap(obejectlist[[i+2]], pattern = "incoming", signaling = pathway.union, title = names(obejectlist)[i+2],width = 5, height = 25)
draw(ht4 + ht5+ht6, ht_gap = unit(0.5, "cm"))
##分别找出各个组别的前20通路
pathway.vehicle<-obejectlist[[i]]@netP$pathways
pathway.3weeks<-obejectlist[[i+1]]@netP$pathways
pathway.6weeks<-obejectlist[[i+2]]@netP$pathways
ht7 = netAnalysis_signalingRole_heatmap(obejectlist[[i]], pattern = "incoming", signaling = pathway.vehicle, title = names(obejectlist)[i], width = 5, height = 25)
ht8 = netAnalysis_signalingRole_heatmap(obejectlist[[i+1]], pattern = "incoming", signaling = pathway.3weeks, title = names(obejectlist)[i+1],width = 5, height = 25)
ht9 = netAnalysis_signalingRole_heatmap(obejectlist[[i+2]], pattern = "incoming", signaling = pathway.6weeks, title = names(obejectlist)[i+2],width = 5, height = 25)
#取各组的前20通路
pathway.top20_1 <- union(pathway.vehicle[1:20],pathway.3weeks[1:20])
pathway.top20<-union(pathway.top20_1,pathway.6weeks[1:20])
ht10 = netAnalysis_signalingRole_heatmap(obejectlist[[i]], pattern = "incoming", signaling = pathway.top20, title = names(obejectlist)[i], width = 5, height = 8)
ht11 = netAnalysis_signalingRole_heatmap(obejectlist[[i+1]], pattern = "incoming", signaling = pathway.top20, title = names(obejectlist)[i+1],width = 5, height = 8)
ht12 = netAnalysis_signalingRole_heatmap(obejectlist[[i+2]], pattern = "incoming", signaling = pathway.top20, title = names(obejectlist)[i+2],width = 5, height = 8)
draw(ht10 + ht11+ht12, ht_gap = unit(0.5, "cm"))
ht13 = netAnalysis_signalingRole_heatmap(obejectlist[[i]], pattern = "outgoing", signaling = pathway.top20, title = names(obejectlist)[i], width = 5, height = 8)
ht14 = netAnalysis_signalingRole_heatmap(obejectlist[[i+1]], pattern = "outgoing", signaling = pathway.top20, title = names(obejectlist)[i+1],width = 5, height = 8)
ht15 = netAnalysis_signalingRole_heatmap(obejectlist[[i+2]], pattern = "outgoing", signaling = pathway.top20, title = names(obejectlist)[i+2],width = 5, height = 8)
draw(ht13 + ht14+ht15, ht_gap = unit(0.5, "cm"))
ht16 = netAnalysis_signalingRole_heatmap(obejectlist[[i]], pattern = "all", signaling = pathway.top20, title = names(obejectlist)[i], width = 5, height = 8)
ht17 = netAnalysis_signalingRole_heatmap(obejectlist[[i+1]], pattern = "all", signaling = pathway.top20, title = names(obejectlist)[i+1],width = 5, height = 8)
ht18 = netAnalysis_signalingRole_heatmap(obejectlist[[i+2]], pattern = "all", signaling = pathway.top20, title = names(obejectlist)[i+2],width = 5, height = 8)
draw(ht16 + ht17+ht18, ht_gap = unit(0.5, "cm"))

#识别上调和下调的信号配体-受体对
levels (cellchatmerge@idents$joint)#看cell group的名称和排序
#Epithelial cell发出到其他类型细胞，
netVisual_bubble(cellchatmerge, sources.use = 4, targets.use = c(1,2,4,5,6,7,8,9,11), 
                 comparison = c(1, 2,3),max.dataset = 1, title.name = "Increased signaling in vehicle",angle.x = 45)
netVisual_bubble(cellchatmerge, sources.use = 4, targets.use = c(1,2,4,5,6,7,8,9,11), 
                 comparison = c(1, 2,3),max.dataset = 3, title.name = "Increased signaling in weeks6",angle.x = 45)
#Fibroblast cell发出到其他类型细胞，
netVisual_bubble(cellchatmerge, sources.use = 5, targets.use = c(1,2,4,5,6,7,8,9,11), 
                 comparison = c(1, 2,3),max.dataset = 1, title.name = "Increased signaling in vehicle",angle.x = 45)
netVisual_bubble(cellchatmerge, sources.use = 5, targets.use = c(1,2,4,5,6,7,8,9,11), 
                 comparison = c(1, 2,3),max.dataset = 3, title.name = "Increased signaling in weeks6",angle.x = 45)

##根据功能相似性对信号进行识别
cellchatmerge <- computeNetSimilarityPairwise(cellchatmerge, type = "functional")
#> Compute signaling network similarity for datasets 1 2 3
install.packages("uwot")
library(uwot)
cellchatmerge <- netEmbedding(cellchatmerge, umap.method = 'uwot', type = "functional")#出现BUG，进行修改
#> Manifold learning of the signaling networks for datasets 1 2 3
cellchatmerge <- netClustering(cellchatmerge, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchatmerge, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2 3
#识别两组间通路的差异性大小，基于功能相似性，较大的距离意味着两个数据集之间的通信网络在功能或结构相似性方面的差异较大
rankSimilarity(cellchatmerge, comparison2 = c(1, 2), type = "functional")#对照和3周
rankSimilarity(cellchatmerge, comparison2 = c(1, 3), type = "functional")#对照和6周

##单个信号通路作图
#COLLAGE的Nnetvisual
pathways.COLLAGEN <- c("COLLAGEN") 
weight.max <- getMaxWeight(obejectlist, slot.name = c("netP"), attribute = pathways.COLLAGEN)
par(mfrow = c(1,3), xpd=TRUE)
for (i in 1:length(obejectlist)) {
  netVisual_aggregate(obejectlist[[i]], signaling = pathways.COLLAGEN,layout = "circle", edge.weight.max = weight.max[1],edge.width.max = 12,signaling.name = paste(pathways.COLLAGEN, names(obejectlist)[i]))
}
#其他的Nnetvisual
pathways.show <- c("CXCL") 
weight.max <- getMaxWeight(obejectlist, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,3), xpd=TRUE)
for (i in 1:length(obejectlist)) {
  netVisual_aggregate(obejectlist[[i]], signaling = pathways.show,layout = "circle", edge.weight.max = weight.max[1],edge.width.max = 8,signaling.name = paste(pathways.show, names(obejectlist)[i]))
}

pathways.show <- c("TGFb") 
weight.max <- getMaxWeight(obejectlist, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,3), xpd=TRUE)
for (i in 1:length(obejectlist)) {
  netVisual_aggregate(obejectlist[[i]], signaling = pathways.show,layout = "circle", edge.weight.max = weight.max[1],edge.width.max = 8,signaling.name = paste(pathways.show, names(obejectlist)[i]))
}
#单个通路上调和下调的信号配体-受体对
levels (cellchatmerge@idents$joint)#看cell group的名称和排序
pathways.show <- c("TGFb") 
#COLLAGEN pathway to epi
netVisual_bubble(cellchatmerge, signaling = pathways.show,sources.use =c(1,2,3,4,5,7,8,9,10,11), targets.use =4, 
                 comparison = c(1, 2,3), title.name = "COLLAGEN pathway to epi",angle.x = 45)
#LAMININ pathway from Fibroblast cell，
netVisual_bubble(cellchatmerge, signaling = pathways.show,sources.use =5 , targets.use =c(1,2,3,4,5,6,7,8,9,10,11) , 
                 comparison = c(1, 2,3), title.name = "LAMININ pathway FROM fibro",angle.x = 45)
#CXCL pathway from Fibroblast cell，
netVisual_bubble(cellchatmerge, signaling = pathways.show,sources.use =5 , targets.use =c(1,3,5,6,7,8,9,11) , 
                 comparison = c(1, 2,3), title.name = "CXCL pathway FROM fibro",angle.x = 45)
#CD45 pathway to Macrophages，
netVisual_bubble(cellchatmerge, signaling = pathways.show,sources.use =c(1,2,6,8,9,11) , targets.use =6, 
                 comparison = c(1, 2,3), title.name = "CD45 pathway to Macrophages",angle.x = 45)
#TGFb pathway to Fibroblast cell，
netVisual_bubble(cellchatmerge, signaling = pathways.show,sources.use =c(1,2,3,4,5,6,7,8,9,10,11) , targets.use =5, 
                 comparison = c(1, 2,3), title.name = "TGFb pathway to Fibroblast cell",angle.x = 45)

