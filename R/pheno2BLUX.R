pheno2BLUX=function(phe_raw_f,phe_BLU_f,BLUP_or_BLUE="BLUP",impute_phe_f=NA){

  # xlsx example:↓

  # id	OA2016GZ	OA2017GZ	OA2018GZ	OA2019GZ	LA2016GZ	LA2017GZ	LA2018GZ	LA2019GZ
  # og1	NA	NA	NA	77.11	NA	NA	NA	15.68
  # og10	50.06	66.17	78.49	86.94	34.63	19.92	8.87	6.81
  # og100	38.73	NA	56.21	64.12	39.81	NA	29.54	30.13
  # og101	42.87	39.3	48.63	59.84	34.8	39.85	34.84	30.76
  # og102	NA	NA	NA	65.49	NA	NA	NA	25.72
  # og103	59.36	51.14	85.48	92.72	19.19	29.31	2.9	0.19
  # og104	NA	45.8	59.73	65.5	NA	35.57	25.51	25.67


  # make sure data in the first sheet
  options (warn = - 1)
  if(!require("pacman")) install.packages("pacman")
  pacman::p_load(openxlsx,pacman,magrittr,xts,quantmod,zoo,abind,ROCR,tidyr,tibble,
                 dplyr,reshape2,purrr,emmeans,lme4,pbkrtest,lmerTest,PerformanceAnalytics)
  DMwR_pkg_url="https://mran.microsoft.com/snapshot/2016-05-02/src/contrib/DMwR_0.4.1.tar.gz"
  if(!require("DMwR")) install.packages(DMwR_pkg_url, repos = NULL, type = "source")
  p_load(DMwR)

  # phe_raw_f="~/proj/advqtl/data/test/phe_raw.fdas"
  if(!(grepl(".csv", phe_raw_f, fixed = TRUE))) stop("phe_raw_f should be a .csv file!")
  phe_raw=read.csv(phe_raw_f)

    sub_id_list=colnames(phe_raw) %>% grep("id",.,value=T,invert = T)
  sub_phe=phe_raw[,c(sub_id_list)]
  rownames(sub_phe)=phe_raw$id

  for (i in sub_id_list) sub_phe[,i]=as.numeric(sub_phe[,i,drop=T])

  if(any(is.na(sub_phe))){
    print("detect NAs in input data, use knn imputation...")
    knnImput_phe=knnImputation(sub_phe, k = 10, scale = T, meth = "weighAvg", distData = NULL)
  }else{
    knnImput_phe=sub_phe
  }

  # out_phe_f="~/Downloads/knnImput_phe.csv"
  if(!(is.na(impute_phe_f)))  write.csv(knnImput_phe,impute_phe_f)

  phe2_df=colnames(knnImput_phe) %>% as.data.frame() %>%
    separate(.,1,c("phe","year"),"20",remove = F) %>%
    as.data.frame() %>%
    set_colnames(c("id","phe","year")) %>%
    mutate(id=as.character(id))
  phe2_df$year=paste0("20",phe2_df$year) %>% gsub("GZ$","",.)

  if(BLUP_or_BLUE=="BLUE"){
    res_raw=""
    for(i in unique(phe2_df$phe)){
      print(i)
      sub_phe2=phe2_df %>% filter(phe==i)
      subid=sub_phe2$id
      sub_knnImput_phe=knnImput_phe %>% rownames_to_column("line") %>% .[,c("line",subid)]
      sub_knnImput_phe_long=melt(sub_knnImput_phe,id=c("line"))
      sub_df=merge(phe2_df,sub_knnImput_phe_long,by.y="variable",by.x="id")
      cols=3:4
      sub_df[,cols]=sub_df%>% select(all_of(cols)) %>% map_df(as.factor)
      # print(dim(sub_df))
      # colnames(sub_df)[5]="OA"
      # m1=lmer(value ~ line + (1|year) + (1|line:year), data=sub_df)
      m1=lmer(value ~ line + (1|year), data=sub_df)
      summary(m1)
      re1=emmeans(m1,"line") %>% as.data.frame()
      re2=re1 %>% .[,1:2] %>% mutate(phe=i)
      res_raw=suppressWarnings(rbind(res_raw,re2))
    }
  }else{
    res_raw=""
    for(i in unique(phe2_df$phe)){
      print(i)
      sub_phe2=phe2_df %>% filter(phe==i)
      subid=sub_phe2$id
      sub_knnImput_phe=knnImput_phe %>% rownames_to_column("line") %>% .[,c("line",subid)]
      sub_knnImput_phe_long=melt(sub_knnImput_phe,id=c("line"))
      sub_df=merge(phe2_df,sub_knnImput_phe_long,by.y="variable",by.x="id")
      cols=3:4
      sub_df[,cols]=sub_df%>% select(all_of(cols)) %>% map_df(as.factor)
      print(dim(sub_df))
      # colnames(sub_df)[5]="OA"
      # m1=lmer(value ~ line + (1|year) + (1|line:year), data=sub_df)
      m1=lmer(value ~ (1|line) + (1|year), data=sub_df)
      summary(m1)
      re1=ranef(m1)$line %>% as.data.frame() %>% rownames_to_column("line") %>% set_colnames(c("line","BLUP"))
      re2=re1 %>% .[,1:2] %>% mutate(phe=i)
      res_raw=suppressWarnings(rbind(res_raw,re2))
    }
  }

  # not sure select is BLUE or BLUP, use BLU instead
  res=res_raw[-1,] %>% as.data.frame() %>%
    set_colnames(c("line","BLU","phe")) %>%
    mutate(BLU=as.numeric(BLU)) %>% as.data.frame() %>%
    mutate(line=as.character(line))

  # saveRDS(res,"~/res_phe6_BLUE.rds")
  # res=readRDS("~/res_phe6_BLUE.rds")
  res2=dcast(data=res,line ~ phe,value.var="BLU") %>% column_to_rownames("line") %>% as.matrix()
  write.csv(res2,phe_BLU_f)

  # chart.Correlation(res2,histogram = T,pch=19)
  # return(res2)
}


# phe_raw_f="~/proj/advqtl/data/phe_raw.xlsx"
# impute_phe_f="~/impute.csv"
# phe_BLU_f="~/BLUE.csv"
# BLUP_or_BLUE="BLUE"

# s=pheno2BLUX("~/proj/advqtl/data/phe_raw.xlsx","~/imput.csv","~/BLUE.csv","BLUE")
