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

