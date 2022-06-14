# qtl2 pipe
# only for crosstype: riself

setwd("~/proj/advqtl/R")
rm(list=ls())
# install.packages("BiocManager")
pacman::p_load(BiocManager,stringr,pacman,qtl,tidyr,tibble,magrittr,data.table,reshape2,qtl2,yaml)

# input
pheno_raw_f="~/proj/advqtl/data/res_phe6_BLUE.rds"
gt_raw_f="~/proj/advqtl/data/gt.csv"
gmap_raw_f="~/proj/advqtl/data/all.map.order.csv"
qtl2_dir="~/proj/advqtl/data"



if(!dir.exists(qtl2_dir)) dir.create(qtl2_dir)
phe_f=paste0(qtl2_dir,"/pheno.csv")
gt_f=paste0(qtl2_dir,"/geno.csv")
gmap_f=paste0(qtl2_dir,"/gmap.csv")

if(!(file.exists(phe_f))){
  phe=readRDS(pheno_raw_f) %>%
    dcast(data=.,line ~ phe,value.var="BLUE")
  colnames(phe)[1]="id"
  write.csv(phe,phe_f,quote = F,row.names = F)
}

if(!(file.exists(phe_f))){
  gt_raw=read.csv(gt_raw_f)
  gt=gt_raw[-c(1:2),]
  write.csv(gt,gt_f,quote = F,row.names = F)
}

if(!(file.exists(gmap_f))){
  gmap=read.csv(gmap_raw_f) %>% dplyr::select(-c("phy_posi")) %>%
    set_colnames(c("marker","chr","pos"))
  write.csv(gmap,gmap_f,quote = F,row.names = F)
}

date_str = function(full=F){
  if(full==T){
    return(paste(gsub("-","_",Sys.Date()),"_",paste(unlist(str_split(unlist(str_split(Sys.time()," "))[2],":")),collapse = "_"),sep = ""))
  }else{
    return(gsub("-","_",Sys.Date()))
  }
}

config_f=paste0(date_str(T),".yaml")
cat(paste0("see config yaml file in ",qtl2_dir,"/",config_f))
cat("\n")
setwd(qtl2_dir)
write("crosstype: riself",config_f)
write("geno: geno.csv",config_f,append = T)
write("pheno: pheno.csv",config_f,append = T)
write("gmap: gmap.csv",config_f,append = T)
write("alleles:",config_f,append = T)
write("- a",config_f,append = T)
write("- b",config_f,append = T)
write("genotypes:",config_f,append = T)
write("  a: 1",config_f,append = T)
write("  h: 2",config_f,append = T)
write("  b: 3",config_f,append = T)
write("na.strings:",config_f,append = T)
write("- '-'",config_f,append = T)
write("- NA",config_f,append = T)


s=qtl2::read_cross2(config_f)

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
find_peaks(out_pg_loco, map, threshold=2.71, peakdrop=1, drop=1,expand2markers=FALSE)

lod_int(out, map, lodcolumn="LA", peakdrop=1, drop=1,expand2markers=FALSE)
bayes_int(out, map, lodcolumn="OA", prob=0.95, expand2markers=FALSE)

s=out %>% as.data.frame()


plot(out, map, lodcolumn=2, col="blue", main=colnames(s$pheno)[2])






