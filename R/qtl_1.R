setwd("~/Downloads")
rm(list=ls())
pacman::p_load(pacman,qtl,tidyr,tibble,magrittr,data.table,reshape2,qtl2)

gt_dir="~/Downloads/qtl/gt"
phe_dir="~/Downloads/qtl/phe"
blue_raw=readRDS("~/res_phe6_BLUE.rds")
dir.create(gt_dir,recursive = T)
dir.create(phe_dir,recursive = T)

blue_phe=dcast(data=blue_raw,line ~ phe,value.var="BLUE") %>% 
  column_to_rownames("line") %>% as.matrix()

gt_raw_f="~/Downloads/gt_raw.csv"
gt_f=paste0(gt_dir,"/gt.csv")
gt_raw=read.csv(gt_raw_f)
cM_ref_f=paste0(gt_dir,"/all.map.order.csv")
cM_ref=read.csv(cM_ref_f) %>% column_to_rownames("Locus")

gt_raw2=gt_raw %>% column_to_rownames("Locus")
gt_raw2$Pos=cM_ref[gt_raw$Locus,"genetic_posi"]

gt=gt_raw2%>% t() %>% as.data.frame() %>%
  .[-c(3:5),] %>% .[c("Chr","Pos",rownames(blue_phe)),] %>% rownames_to_column("id")
gt[1,1]=""
gt[2,1]=""

write.csv(gt,file = gt_f,quote = F,row.names = F)

for(i in colnames(blue_phe)){
  print(i)
  sub_phe_f=paste0(phe_dir,"/",i,".csv")
  sub_phe=blue_phe[,i,drop=F] %>% as.data.frame() %>% rownames_to_column("id") %>%
    set_colnames(c("id","phe")) %>% .[,c("phe","id")]
  write.csv(sub_phe,sub_phe_f,quote = F,row.names = F)
  s=qtl::read.cross("csvs", "", gt_f, sub_phe_f, genotypes=c("a","h","b","d","c"), 
                    alleles=c("a","b"),)
  
  # for(i in 3:NROW(gt)){
  #   # print(gt[i,1])
  #   if(length(unique(unname(t(gt[i,]))))>4) print(gt[i,1])
  # }
# 
#   summary(s)
#   plotMap(s)
#   
#   plotPheno(s, pheno.col=1)
  
  ## 计算基因型概率
  s <- calc.genoprob(s, step=1)
  ## 使用默认方法进行single-QTL全基因组扫描
  # out.em <- scanone(s)
  # ## 查看扫描结果
  # summary(out.em)
  # ## 挑选LOD > 3的结果
  # summary(out.em, threshold=3)
  # ## 展示结果
  # plot(out.em)
  
  
  out.hk <- scanone(s, method="hk")
  mythreshold=0.2
  
  ## 进行1000次Permutation test
  operm <- scanone(s, method="hk", n.perm=1000)
  # operm <- scanone(s, n.perm=1000)
  # ## 获得显著性阈值
  # summary(operm, alpha=c(0.05, 0.2))
  summary(operm, alpha=c(mythreshold))
  summary(operm, 0.05)
  
  
  ## 从扫描结果中挑选显著的位点
  summary(out.hk, perms=operm, alpha=0.2, pvalues=TRUE)
  summary(out.hk, threshold=2.5,format="allpeaks")
  
  scanoneFiltGenedenovo=function(scanoneOutObj,lod_threshold){
    p_load(dplyr)
    # scanoneOutObj=out.hk
    # lod_threshold=4.51
    # lod_threshold=3.40051814931523
    M=scanoneOutObj %>% as.data.frame()
    filtM=M %>% mutate(chr=as.character(chr)) %>% dplyr::filter(lod>lod_threshold)
    chr_list=filtM$chr %>% unique()
    if(length(chr_list>0)){
      for(each_chr in chr_list){
        subFiltM=filtM %>% filter(chr==each_chr)
        
      }
      
    }
    
  }
}






qtl::read.cross()

?read.cross

## 安装R包
install.packages("qtl")
## 加载R包
library("qtl")
## 导入基因型和表型数据
sug <- read.cross("csvs", ".", "gen.csv", "phe.csv")
## 查看输入文件相关信息
summary(sug)

## 样本数
nind(sug)
## 染色体数
nchr(sug)
## 标记数
totmar(sug)
## 每个染色体上的标记数
nmar(sug)
## 表型数
nphe(sug)
# 除了文字信息外，我们还可以用图来展示这些信息。
plot(sug)


## 展示缺失基因型数据（黑色为缺失的基因型）
plotMissing(sug)

## 绘制遗传图谱
plotMap(sug)

## 绘制表型分布直方图
plotPheno(sug, pheno.col=1)


## 计算基因型概率
sug <- calc.genoprob(sug, step=1)
## 使用默认方法进行single-QTL全基因组扫描
out.em <- scanone(sug)
## 查看扫描结果
summary(out.em)
## 挑选LOD > 3的结果
out.em2=out.em %>% filter(chr=="CM014323.1")
summary(out.em, threshold=3)
## 展示结果
plot(out.em2)

?cim





