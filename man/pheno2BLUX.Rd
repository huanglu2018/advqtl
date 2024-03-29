% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pheno2BLUX.R
\name{pheno2BLUX}
\alias{pheno2BLUX}
\title{pheno2BLUX calculates BLUP and BLUE values for data with multiple trials}
\usage{
pheno2BLUX(phe_raw_f, phe_BLU_f, BLUP_or_BLUE = "BLUP", impute_phe_f = NA)
}
\arguments{
\item{phe_raw_f}{csv file with the structure below}

\item{phe_BLU_f}{BLUP or BLUE value data in .csv format see data(BLUE_example)}

\item{BLUP_or_BLUE}{to get BLUP or BLUE values, default "BLUP"}

\item{impute_phe_f}{where to write down the knn-immputed phe data, default NA}
}
\value{
no value return, phe_BLU_f is what you need
}
\description{
pheno2BLUX calculates BLUP and BLUE values for data with multiple trials
}
\examples{

phe_raw_f (.csv) example:↓
id	OA2016GZ	OA2017GZ	OA2018GZ	OA2019GZ	LA2016GZ	LA2017GZ	LA2018GZ	LA2019GZ
og1	NA	NA	NA	77.11	NA	NA	NA	15.68
og10	50.06	66.17	78.49	86.94	34.63	19.92	8.87	6.81
og100	38.73	NA	56.21	64.12	39.81	NA	29.54	30.13
og101	42.87	39.3	48.63	59.84	34.8	39.85	34.84	30.76
og102	NA	NA	NA	65.49	NA	NA	NA	25.72
og103	59.36	51.14	85.48	92.72	19.19	29.31	2.9	0.19
og104	NA	45.8	59.73	65.5	NA	35.57	25.51	25.67
}
