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
  pacman::p_load(pacman,ggplotify,ggplot2,BiocManager,stringr,pacman,qtl,tidyr,tibble,magrittr,data.table,reshape2,eqtl)

  input_phe_file=paste0(tempdir(),"/input_phe.csv")
  input_gt_file=paste0(tempdir(),"/input_gt.csv")

  phe_raw=read.csv(phe_file)
  colnames(phe_raw)[1]="id"

  gt_raw=read.csv(gt_file)
  id_inter=intersect(phe_raw$id,colnames(gt_raw)[7:NCOL(gt_raw)])

  print(paste0("A total of ",length(id_inter)," samples kept after checking for intersect."))

  phe_raw2=phe_raw %>% column_to_rownames("id") %>% .[id_inter,] %>% rownames_to_column("id")
  gt_raw1=cbind(gt_raw[,1:6],gt_raw[,id_inter])

  phe=cbind(phe_raw2[,-1],phe_raw2[,1,drop=F])

  # merge gt and gmap file
  cM_ref=read.csv(gmap_file) %>% column_to_rownames("Locus")
  gt_raw2=gt_raw1 %>% column_to_rownames("Locus")
  gt_raw2$Pos=cM_ref[gt_raw1$Locus,"genetic_posi"]
  gt=gt_raw2%>% t() %>% as.data.frame() %>%
    .[-c(3:5),] %>% .[c("Chr","Pos",phe$id),] %>% rownames_to_column("id")
  gt[1,1]=""
  gt[2,1]=""

  write.csv(phe,file = input_phe_file,quote = F,row.names = F)
  write.csv(gt,file = input_gt_file,quote = F,row.names = F)

  s=qtl::read.cross("csvs", "",
                    input_gt_file, input_phe_file,
                    genotypes=c("a","h","b","d","c"),
                    alleles=c("a","b"),crosstype=crosstype)
  s <- calc.genoprob(s, step=1)
  out.em <- scanone(s, pheno.col=1:(NCOL(phe)-1), model='normal', method='hk')

  res_raw=""
  pic_list=list()

  for (i in names(s$pheno)[-c(length(s$pheno))]){
    print(i)
    operm <- scanone(s, pheno.col= i,method="hk", n.perm=1000)

    out_i=scanone(s, pheno.col=i, model='normal', method='hk')
    out_i$chr=as.character(out_i$chr)
    p = ggplot(out_i,aes(pos,lod))+
      geom_line()+
      facet_grid(. ~chr, scales = "free", space = "free")+
      theme_bw()

    pic_list[[i]]=p

    lod_threshold=summary(operm, alpha=c(pval_threshold)) %>% as.numeric()
    print(paste0(i," get a lod threshold of ",format(lod_threshold,digit=4),
                 " when the pvalue was set to ",pval_threshold))
    each_peaklist=eqtl::define.peak(out.em,lodcolumn=i,th=lod_threshold,si=1.5,graph=FALSE,window.size=10)
    subRes_raw=""
    for(j in names(each_peaklist[[i]])){
      if(is.data.frame(each_peaklist[[i]][[j]])){
        each_chr_peak=each_peaklist[[i]][[j]]
        each_chr_peak2=cbind(trait=i,chr=j,each_chr_peak)
        subRes_raw=rbind(subRes_raw,each_chr_peak2)
      }
    }
  #### get pve values for each trait
    resGetPve=function(res,crossObj,pheno.col){
      sigQtlNumber=res %>% as.data.frame() %>% remove_rownames() %>% NROW()
      myformula=as.formula(paste0("y ~ ",paste0("Q",seq(sigQtlNumber),collapse = " + ")))
      mqtl=makeqtl(crossObj, chr=res[,2], pos=as.numeric(res[,5]), what=c("prob"))
      fqtl<-fitqtl(crossObj,pheno.col=pheno.col,dropone=T,get.ests=T,
                   model="normal", qtl=mqtl,method="hk",formula=myformula)
      if(sigQtlNumber>1){
        pveRes=summary(fqtl)$result.drop[,c(4),drop=F] %>% as.data.frame() %>% set_colnames(c("pve"))
      }else{
        pveRes=summary(fqtl)$result.full[1,c(5),drop=F] %>% as.data.frame() %>% set_colnames(c("pve"))
      }
      resOut1=cbind(res[,1:5],pveRes[,1,drop=F])
      resOut2=cbind(resOut1,res[,6:NCOL(res)])
      return(resOut2)
    }

    if(!(is.data.frame(subRes_raw))) next
    subRes_raw2=subRes_raw[-1,] %>% as.data.frame()
    subRes=resGetPve(subRes_raw2,s,i)
    res_raw=rbind(res_raw,subRes)
  }
  if(!(is.data.frame(res_raw))) stop("no significant loci for any of the traits !!!")
  res=res_raw[-1,] %>% as.data.frame()
  return(list(sig_res=res,pic=pic_list))
}



#
# ## wrf input
#
# s=call_IM_peaks("~/proj/advqtl/data/test/gt.csv",
#                 "~/proj/advqtl/data/test/gmap.csv",
#                 "~/proj/advqtl/data/test/phe.csv",
#                 0.05,"riself")
#
# write.csv(s,"~/proj/advqtl/data/test/res.csv")
#
# gt_file="~/proj/advqtl/data/test/gt.csv"
# gmap_file="~/proj/advqtl/data/test/gmap.csv"
# phe_file="~/proj/advqtl/data/test/phe.csv"
# pval_threshold=0.05
# crosstype="riself"
#
# ## gda input
#
# s=call_IM_peaks("~/proj/advqtl/data/test2/gt.csv",
#                 "~/proj/advqtl/data/test2/all.map.order.csv",
#                 "~/proj/advqtl/data/test2/phe_GdaImpute.csv",
#                 0.05,"riself")
#
# write.csv(s,"~/proj/advqtl/data/test2/res_reproduce.csv")
#
# gt_file="~/proj/advqtl/data/test2/gt.csv"
# gmap_file="~/proj/advqtl/data/test2/all.map.order.csv"
# phe_file="~/proj/advqtl/data/test2/phe_GdaImpute.csv"
# pval_threshold=0.05
# crosstype="riself"
#
#

