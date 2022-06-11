# qtl2 pipe

setwd("~/proj/advqtl/R")
rm(list=ls())
install.packages("BiocManager")
pacman::p_load(pacman,qtl,tidyr,tibble,magrittr,data.table,reshape2,qtl2,yaml)

#
phe=readRDS("~/proj/advqtl/data/res_phe6_BLUE.rds") %>%
  dcast(data=.,line ~ phe,value.var="BLUE")
colnames(phe)[1]="id"
write.csv(phe,"~/proj/advqtl/data/pheno.csv",quote = F,row.names = F)

gt_f="~/proj/advqtl/data/gt.csv"
gt_raw=read.csv(gt_f)
gt=gt_raw[-c(1:2),]
write.csv(gt,"~/proj/advqtl/data/geno.csv",quote = F,row.names = F)


gmap_f="~/proj/advqtl/data/all.map.order.csv"
gmap=read.csv(gmap_f) %>% dplyr::select(-c("phy_posi")) %>%
  set_colnames(c("marker","chr","pos"))
write.csv(gmap,"~/proj/advqtl/data/gmap.csv",quote = F,row.names = F)



# write_yaml(data.frame(a=1:10, b=letters[1:10], c=11:20), "~/Downloads/qtl/gt/test2.yaml")



s=qtl2::read_cross2("~/proj/advqtl/data/test.yaml")

map <- insert_pseudomarkers(s$gmap, step=1)
probs <- calc_genoprob(s, map, error_prob=0.002)

pheno <- s$pheno
# covar <- match(s$covar$sex, c("f", "m")) # make numeric
# names(covar) <- rownames(s$covar)
# Xcovar <- get_x_covar(s)
# perform genome scan
out <- scan1(probs, pheno)
kinship <- calc_kinship(probs)
out_pg <- scan1(probs, pheno, kinship)

kinship_loco <- calc_kinship(probs, "loco")
out_pg_loco <- scan1(probs, pheno, kinship_loco)
# leave-one-chromosome-out kinship matrices
# kinship <- calc_kinship(probs, "loco")

# genome scan with a linear mixed model
# out_lmm <- scan1(probs, pheno, kinship)

# find just the highest peak on each chromosome
# find_peaks(out, map, threshold=3)

# possibly multiple peaks per chromosome
# find_peaks(out, map, threshold=3, peakdrop=1)

# possibly multiple peaks, also getting 1-LOD support intervals
# find_peaks(out, map, threshold=3, peakdrop=1, drop=1)

# possibly multiple peaks, also getting 90% Bayes intervals
# find_peaks(out, map, threshold=3, peakdrop=1, prob=0.9)

# ?find_peaks
operm <- scan1perm(probs, s$pheno, n_perm=1000)
summary(operm,alpha=0.2)
# LOD thresholds (1000 permutations)
# LA   OA Oil   PA Prot   SA
# 0.05 3.46 3.48 3.4 3.46 3.41 3.33

find_peaks(out, map,lodcolumn=2, threshold=2.72, peakdrop=1, drop=1,expand2markers=FALSE)
find_peaks(out_pg, map, threshold=2.72, peakdrop=1, drop=1)
find_peaks(out_pg_loco, map, threshold=2.72, peakdrop=1, drop=1,expand2markers=FALSE)

lod_int(out, map, lodcolumn="LA", peakdrop=1, drop=1,expand2markers=FALSE)
bayes_int(out, map, lodcolumn="OA", prob=0.95, expand2markers=FALSE)

s=out %>% as.data.frame()


plot(out, map, lodcolumn=2, col="blue", main=colnames(s$pheno)[2])










blue_raw=readRDS("~/res_phe6_BLUE.rds")

gt_dir="~/Downloads/qtl/gt"
phe_dir="~/Downloads/qtl/phe"

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


