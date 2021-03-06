#' call_IM_peaks use genotype file, loci mapping file and phenotype file (can be the result file of pheno2BLUX)
#'
#' @param gt_file three csv files with the structure below
#' @param gmap_file three csv files with the structure below
#' @param phe_file three csv files with the structure below
#' @param pval_threshold pvalue for the determination of threshold lod value, default 0.05
#' @param crosstype one from "bc","f2","risib","riself","4way"
#'
#'
#' @return example see data(result_example)
#' @export
#'
#' @examples
#'
#' gt_file (.csv) example:
#' Chr,Pos,p1(P),p2(M),Locus,Classification,og1,og10,og100
#' CM014315.1,1694837,0|0,1|1,CM014315.1_1694837,"(a,h,b)",b,a,b
#' CM014315.1,1932696,0|0,1|1,CM014315.1_1932696,"(a,h,b)",b,a,b
#' CM014315.1,2058342,0|0,1|1,CM014315.1_2058342,"(a,h,b)",b,a,b
#'
#' gmap_file (.csv) example:
#' Locus,LG,genetic_posi,phy_posi
#' CM014315.1_1694837,CM014315.1,0,1694837
#' CM014315.1_1932696,CM014315.1,0.25,1932696
#' CM014315.1_2058342,CM014315.1,0.5,2058342
#' CM014315.1_2400548,CM014315.1,0.75,2400548
#'
#' phe_file (.csv) example:
#' Line,OA,LA,Oil,Prot
#' og1,64.504,21.069,49.792,24.649
#' og10,70.415,17.5575,45.94,28.79
#' og100,50.626,34.521,51.175,27.577
#'
#'
#'




call_IM_peaks=function(gt_file, gmap_file, phe_file, pval_threshold=0.05, crosstype=""){

  crosslist=c("bc","f2","risib","riself","4way")
  if (!(crosstype %in% crosslist)) stop(paste0("crosstype canc only be one of the following: ",paste0(crosslist,collapse = " ")))
  # gt_file="~/proj/advqtl/data/test/gt_raw.csv"
  # gmap_file="~/proj/advqtl/data/test/all.map.order.csv"
  # phe_file="~/proj/advqtl/data/test/BLUE.csv"
  if(!(grepl(".csv", gt_file, fixed = TRUE))) stop("gt_file should be a .csv file!")
  if(!(grepl(".csv", gmap_file, fixed = TRUE))) stop("gmap_file should be a .csv file!")
  if(!(grepl(".csv", phe_file, fixed = TRUE))) stop("phe_file should be a .csv file!")

  # load required packages
  wqtl_pkg_url="https://cran.r-project.org/src/contrib/Archive/eqtl/eqtl_1.1-7.tar.gz"
  if(!require(eqtl)) install.packages(wqtl_pkg_url, repos = NULL, type = "source")
  if(!require(pacman)) install.packages("pacman")
  pacman::p_load(pacman,BiocManager,stringr,pacman,qtl,tidyr,tibble,magrittr,data.table,reshape2,eqtl)

  input_phe_file=paste0(tempdir(),"/input_phe.csv")
  input_gt_file=paste0(tempdir(),"/input_gt.csv")

  phe_raw=read.csv(phe_file)
  colnames(phe_raw)[1]="id"
  phe=cbind(phe_raw[,-1],phe_raw[,1,drop=F])

  # merge gt and gmap file
  gt_raw=read.csv(gt_file)
  cM_ref=read.csv(gmap_file) %>% column_to_rownames("Locus")
  gt_raw2=gt_raw %>% column_to_rownames("Locus")
  gt_raw2$Pos=cM_ref[gt_raw$Locus,"genetic_posi"]
  gt=gt_raw2%>% t() %>% as.data.frame() %>%
    .[-c(3:5),] %>% .[c("Chr","Pos",phe$id),] %>% rownames_to_column("id")
  gt[1,1]=""
  gt[2,1]=""

  write.csv(phe,file = input_phe_file,quote = F,row.names = F)
  write.csv(gt,file = input_gt_file,quote = F,row.names = F)

  s=qtl::read.cross("csvs", "",
                    input_gt_file, input_phe_file,
                    genotypes=c("a","h","b","d","c"),
                    alleles=c("a","b"),crosstype="riself")
  s <- calc.genoprob(s, step=1)
  out.em <- scanone(s, pheno.col=1:(NCOL(phe)-1), model='normal', method='hk')

  res_raw=""

  for (i in names(s$pheno)[-c(length(s$pheno))]){
    print(i)
    operm <- scanone(s, pheno.col= i,method="hk", n.perm=1000)
    lod_threshold=summary(operm, alpha=c(pval_threshold)) %>% as.numeric()
    print(paste0(i," get a lod threshold of ",format(lod_threshold,digit=4),
                 " when the pvalue was set to ",pval_threshold))
    each_peaklist=eqtl::define.peak(out.em,lodcolumn=i,th=lod_threshold,si=1.5,graph=FALSE,window.size=10)
    for(j in names(each_peaklist[[i]])){
      if(is.data.frame(each_peaklist[[i]][[j]])){
        each_chr_peak=each_peaklist[[i]][[j]]
        each_chr_peak2=cbind(trait=i,each_chr_peak)
        res_raw=rbind(res_raw,each_chr_peak2)
      }
    }
  }

  res=res_raw[-1,] %>% as.data.frame()
  return(res)
}

# s=call_IM_peaks("~/proj/advqtl/data/test/gt_raw.csv",
#                 "~/proj/advqtl/data/test/all.map.order.csv",
#                 "~/proj/advqtl/data/test/phe_raw_official.csv",
#                 0.05,"riself")
#
# write.csv(s,"~/proj/advqtl/data/test/res_reproduce.csv")

