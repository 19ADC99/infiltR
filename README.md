# infiltR

Andrea Del Cortona

# Introduction

infiltR is a toolkit for evaluating infiltrating immune cells in tumor RNA-seq samples. It includes three common bulk RNA-seq deconvolution tools for quantification of tumor-infiltrating immune cells: [CIBERSORT](https://doi.org/10.1038/nmeth.3337), [MCP-counter](https://doi.org/10.1186/s13059-016-1070-5) and [quanTIseq](https://doi.org/10.1186/s13073-019-0638-6). MCP-counter is used to obtain absolute quantification of infiltrating immune cells, while a more detailed profiling of relative immune cells subtypes is performed with CIBERSORT and quanTIseq.

infiltR requires a matrix of gene counts. Lowly expressed genes should be removed (e.g.: with `r edgeR::filterByExpr()`) and counts can be linearly transformed, e.g.: `r TMM`-normalized. As benchmarked by [Avila Cobos et al.](https://doi.org/10.1038/s41467-020-19015-1) and [Zhong et al.](https://doi.org/10.1038/nmeth.1830), logarithmic transformation led to worse results when performing computational deconvolution than the linear (un-transformed) data. Next to the count matrix, infiltR takes a metadata table, where row names corresponds to the column names of the count matrix, with specified a sample_groups column that indicates the groups infiltR will use to compare the relative abundances of infiltrating immune cells across samples. An example of sample_groups could be: "treatment1", "treatment2", "untreated", or "responsive", "non-responsive".

infiltR first runs each deconvolution tool separately. Then, it generates publication-ready figures where extensive comparisons of absolute and relative abundance of infiltrating-immune cells between the sample groups are performed. Lastly, it reports on the concordance between scaled relative abundance of common infiltrating immune cells shared between the three methods (CIBERSORT, MCP-counter, quanTIseq).


# Installation

## 1. Install the package from github

You can install the development version of infiltR from [github](https://github.com/19ADC99/infiltR) with:

``` r
install.packages("devtools")
devtools::install_github("19ADC99/infiltR")
```

## 2. Load the package into R session

Load infitR to your R session:

``` r
library("infiltR")
```


# Quick start

``` r
# load library
library(infiltR)

# run infiltR
infiltR(
  counts_table,
  metadata,
  "cancer_type",
  my_palette = "default",
  save_plots = TRUE,
  outdir = "./",
)

```



# Acknowledgments



# License

infiltR is provided as GLP-3 package.

According to the [CIBERSORTx website](https://cibersortx.stanford.edu/#myModalagree):

> The Board of Trustees of the Leland Stanford Junior University (“Stanford”) provides CIBERSORTx website features and services (“Service”) free of charge for non-commercial use only. Use of the Service by any commercial entity for any purpose, including research, is prohibited.

Given the statement above on CIBERSORT license, I have retroengineered the original nu-SVM algorithm described by [Newman et al., 2015](https://doi.org/10.1038/nmeth.3337), inspired by the work of [Moonerss](https://github.com/Moonerss/). The algorithm is provided as-is and I cannot guarantee that it performs as the original CIBERSORT tool since I do not have access to it.


# Session Info


``` r
utils::sessionInfo()
```

```
R version 4.4.2 (2024-10-31)
Platform: x86_64-suse-linux-gnu
Running under: openSUSE Tumbleweed

Matrix products: default
BLAS:   /usr/lib64/R/lib/libRblas.so 
LAPACK: /usr/lib64/R/lib/libRlapack.so;  LAPACK version 3.12.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=it_IT.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: Europe/Brussels
tzcode source: system (glibc)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] infiltR_0.1.0

loaded via a namespace (and not attached):
 [1] SummarizedExperiment_1.34.0 gtable_0.3.6                ggplot2_3.5.1              
 [4] rstatix_0.7.2               Biobase_2.64.0              lattice_0.22-6             
 [7] quadprog_1.5-8              vctrs_0.6.5                 tools_4.4.2                
[10] doSNOW_1.0.20               generics_0.1.3              stats4_4.4.2               
[13] parallel_4.4.2              tibble_3.2.1                proxy_0.4-27               
[16] fansi_1.0.6                 pkgconfig_2.0.3             ggpp_0.5.8-1               
[19] Matrix_1.7-1                S4Vectors_0.42.1            lifecycle_1.0.4            
[22] GenomeInfoDbData_1.2.12     stringr_1.5.1               compiler_4.4.2             
[25] MatrixModels_0.5-3          munsell_0.5.1               codetools_0.2-20           
[28] carData_3.0-5               limSolve_1.5.7.1            SparseM_1.84-2             
[31] GenomeInfoDb_1.40.1         quantreg_5.99               snow_0.4-4                 
[34] class_7.3-22                Formula_1.2-5               preprocessCore_1.66.0      
[37] car_3.1-3                   tidyr_1.3.1                 pillar_1.9.0               
[40] crayon_1.5.3                MASS_7.3-61                 DelayedArray_0.30.1        
[43] iterators_1.0.14            abind_1.4-8                 foreach_1.5.2              
[46] quantiseqr_1.12.0           tidyselect_1.2.1            stringi_1.8.4              
[49] purrr_1.0.2                 dplyr_1.1.4                 splines_4.4.2              
[52] cowplot_1.1.3               grid_4.4.2                  SparseArray_1.4.8          
[55] colorspace_2.1-1            cli_3.6.3                   magrittr_2.0.3             
[58] S4Arrays_1.4.1              survival_3.7-0              utf8_1.2.4                 
[61] broom_1.0.7                 e1071_1.7-16                backports_1.5.0            
[64] ggpmisc_0.6.1               scales_1.3.0                UCSC.utils_1.0.0           
[67] lubridate_1.9.3             timechange_0.3.0            XVector_0.44.0             
[70] httr_1.4.7                  matrixStats_1.4.1           lpSolve_5.6.21             
[73] GenomicRanges_1.56.2        IRanges_2.38.1              rlang_1.1.4                
[76] glue_1.8.0                  polynom_1.4-1               MCPcounter_1.2.0           
[79] BiocGenerics_0.50.0         rstudioapi_0.17.1           jsonlite_1.8.9             
[82] R6_2.5.1                    MatrixGenerics_1.16.0       zlibbioc_1.50.0
```
