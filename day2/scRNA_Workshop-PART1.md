---
title: "single cell 10x single-cell analysis - part1"
author: "UC Davis Bioinformatics Core"
date: "3/20/2018"
output:
  html_document:
    keep_md: true
---

[Seurat](http://satijalab.org/seurat/) is a popular R package that is designed for QC, analysis, and exploration of single cell RNA-seq data. Seurat aims to enable users to identify and interpret sources of heterogeneity from single cell transcriptomic measurements, and to integrate diverse types of single cell data. Further, the authors provide several [tutorials](http://satijalab.org/seurat/get_started.html), on their website.

Dowload and expand the expression_tables.tar.gz file to extract the single cell matrix files for the three samples. These are isolated mouse cells ran on the 10X genomics platform for single cell RNA sequencing, sequenced with UC Davis on 1 HiSeq 4000.

* sample1, UCD_VitE_Def
* sample2, UCD_Supp_VitE
* sample3, UCD_Adj_VitE

We start with loading needed libraries for R, at this time all we need is the package [Seurat](http://satijalab.org/seurat/).

```
## Loading required package: ggplot2
```

```
## Loading required package: cowplot
```

```
## 
## Attaching package: 'cowplot'
```

```
## The following object is masked from 'package:ggplot2':
## 
##     ggsave
```

```
## Loading required package: Matrix
```

## Load the Cell Ranger Matrix Data and create the base Seurat object.
Cell Ranger provides a function `cellranger aggr` that will combine multiple samples into a single matrix file. However, when processing data in R and Seurat this is unnecessary and we can aggregate them in R.

Seurat provides a function `Read10X` to read in 10X data folder. First we read in data from each individual sample folder. First, we initialize the Seurat object (`CreateSeuratObject`) with the raw (non-normalized data). Keep all genes expressed in >= 10 cells. Keep all cells with at least 200 detected genes. Also extracting sample names, calculating and adding in the metadata mitochondrial percentage of each cell. Adding in the metadata batchid. Finally, saving the raw Seurat object.

```r
## Dataset for analysis
dataset_loc <- "expression_tables"
ids <- c("sample1", "sample2", "sample3")

d10x.data <- sapply(ids, function(i){
  d10x <- Read10X(file.path(dataset_loc,i,"outs/filtered_gene_bc_matrices/mm10/"))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
  d10x
})

experiment.data <- do.call("cbind", d10x.data)

experiment.aggregate <- CreateSeuratObject(
  experiment.data,
  project = "scRNA workshop", 
  min.cells = 10,
  min.genes = 200,
  names.field = 2,
  names.delim = "\\-")
```
### Calc mitocondrial content
Calculate percent mitochondrial genes per cell. In mouse these genes can be identified as those that begin with 'mt', in human data they begin with MT.

```r
mito.genes <- grep("^mt-", rownames(experiment.aggregate@data), value = T)
percent.mito <- Matrix::colSums(experiment.aggregate@raw.data[mito.genes, ]) / Matrix::colSums(experiment.aggregate@raw.data)

# AddMetaData adds columns to object@data.info, and is a great place to stash QC stats
experiment.aggregate <- AddMetaData(
  object = experiment.aggregate,
  metadata = percent.mito,
  col.name= "percent.mito")
```
### Lets fix the sample names, reassign names with more meaningful factors

The original samples names (the names above in ids) can be found in the metadata slot, column orig.ident. Here we build a new metadata variable 'batchid' which can be used to specify treatment groups.

```r
samplename = experiment.aggregate@meta.data$orig.ident
table(samplename)
```

```
## samplename
## sample1 sample2 sample3 
##    6956    7855    6874
```

```r
batchid = rep("UCD_VitE_Def",length(samplename))
batchid[samplename %in% c("sample2")] = "UCD_Supp_VitE"
batchid[samplename %in% c("sample3")] = "UCD_Adj_VitE"
names(batchid) = rownames(experiment.aggregate@meta.data)

experiment.aggregate <- AddMetaData(
  object = experiment.aggregate,
  metadata = batchid,
  col.name = "batchid")

table(experiment.aggregate@meta.data$batchid)
```

```
## 
##  UCD_Adj_VitE UCD_Supp_VitE  UCD_VitE_Def 
##          6874          7855          6956
```
### Lets spend a little time getting to know the Seurat object.

The Seurat object is the center of each single cell analysis. It stores __all__ information associated with the dataset, including data, annotations, analyes, etc. The R function slotNames can be used to view the slot names within an object.


```r
slotNames(experiment.aggregate)
```

```
##  [1] "raw.data"     "data"         "scale.data"   "var.genes"   
##  [5] "is.expr"      "ident"        "meta.data"    "project.name"
##  [9] "dr"           "assay"        "hvg.info"     "imputed"     
## [13] "cell.names"   "cluster.tree" "snn"          "calc.params" 
## [17] "kmeans"       "spatial"      "misc"         "version"
```

We can then few the data within a slot with the `@` operator.

```r
head(experiment.aggregate@meta.data)
```

```
##                           nGene nUMI orig.ident percent.mito      batchid
## AAAACCTGAGATCACGG-sample1  1309 2435    sample1  0.075154004 UCD_VitE_Def
## AAAACCTGAGCATCATC-sample1  1331 2650    sample1  0.008301887 UCD_VitE_Def
## AAAACCTGAGCGCTCCA-sample1  1507 2817    sample1  0.038693646 UCD_VitE_Def
## AAAACCTGAGTGGGATC-sample1  1891 3625    sample1  0.005517241 UCD_VitE_Def
## AAAACCTGCAGACAGGT-sample1   786 1091    sample1  0.013748854 UCD_VitE_Def
## AAAACCTGGTATAGGTA-sample1  2505 7248    sample1  0.027179912 UCD_VitE_Def
```


```r
table(experiment.aggregate@meta.data$orig.ident)
```

```
## 
## sample1 sample2 sample3 
##    6956    7855    6874
```

#### Question(s)

1. What slots are empty, what slots have data?
2. What columns are available in meta.data?
3. Look up the help documentation for filter?

## Finally, save the original object, write out a tab-delimited table that could be read into excel, and view the object.

```r
# write.table(as.matrix(experiment.data),"raw.datatable.txt",sep="\t",col.names=T,row.names=T)
experiment.aggregate
```

```
## An object of class seurat in project scRNA workshop 
##  15606 genes across 21685 samples.
```

```r
## Original dataset in Seurat class, with no filtering
save(experiment.aggregate,file="original_seurat_object.RData")
```

## Session Information

```r
sessionInfo()
```

```
## R version 3.4.4 (2018-03-15)
## Platform: x86_64-apple-darwin15.6.0 (64-bit)
## Running under: macOS High Sierra 10.13.3
## 
## Matrix products: default
## BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] Seurat_2.2.1  Matrix_1.2-12 cowplot_0.9.2 ggplot2_2.2.1
## 
## loaded via a namespace (and not attached):
##   [1] diffusionMap_1.1-0   Rtsne_0.13           VGAM_1.0-5          
##   [4] colorspace_1.3-2     ggridges_0.4.1       class_7.3-14        
##   [7] modeltools_0.2-21    mclust_5.4           rprojroot_1.3-2     
##  [10] htmlTable_1.11.2     base64enc_0.1-3      proxy_0.4-21        
##  [13] rstudioapi_0.7       DRR_0.0.3            flexmix_2.3-14      
##  [16] prodlim_1.6.1        mvtnorm_1.0-7        lubridate_1.7.3     
##  [19] ranger_0.9.0         codetools_0.2-15     splines_3.4.4       
##  [22] R.methodsS3_1.7.1    mnormt_1.5-5         robustbase_0.92-8   
##  [25] knitr_1.20           tclust_1.3-1         RcppRoll_0.2.2      
##  [28] Formula_1.2-2        caret_6.0-78         ica_1.0-1           
##  [31] broom_0.4.3          ddalpha_1.3.1.1      cluster_2.0.6       
##  [34] kernlab_0.9-25       R.oo_1.21.0          sfsmisc_1.1-2       
##  [37] compiler_3.4.4       backports_1.1.2      assertthat_0.2.0    
##  [40] lazyeval_0.2.1       lars_1.2             acepack_1.4.1       
##  [43] htmltools_0.3.6      tools_3.4.4          bindrcpp_0.2        
##  [46] igraph_1.2.1         gtable_0.2.0         glue_1.2.0          
##  [49] reshape2_1.4.3       dplyr_0.7.4          Rcpp_0.12.16        
##  [52] trimcluster_0.1-2    gdata_2.18.0         ape_5.0             
##  [55] nlme_3.1-131.1       iterators_1.0.9      fpc_2.1-11          
##  [58] psych_1.7.8          timeDate_3043.102    gower_0.1.2         
##  [61] stringr_1.3.0        irlba_2.3.2          gtools_3.5.0        
##  [64] DEoptimR_1.0-8       MASS_7.3-49          scales_0.5.0        
##  [67] ipred_0.9-6          parallel_3.4.4       RColorBrewer_1.1-2  
##  [70] yaml_2.1.18          pbapply_1.3-4        gridExtra_2.3       
##  [73] segmented_0.5-3.0    rpart_4.1-13         latticeExtra_0.6-28 
##  [76] stringi_1.1.7        foreach_1.4.4        checkmate_1.8.5     
##  [79] caTools_1.17.1       lava_1.6             dtw_1.18-1          
##  [82] SDMTools_1.1-221     rlang_0.2.0          pkgconfig_2.0.1     
##  [85] prabclus_2.2-6       bitops_1.0-6         evaluate_0.10.1     
##  [88] lattice_0.20-35      ROCR_1.0-7           purrr_0.2.4         
##  [91] bindr_0.1.1          recipes_0.1.2        htmlwidgets_1.0     
##  [94] tidyselect_0.2.4     CVST_0.2-1           plyr_1.8.4          
##  [97] magrittr_1.5         R6_2.2.2             gplots_3.0.1        
## [100] Hmisc_4.1-1          dimRed_0.1.0         sn_1.5-1            
## [103] withr_2.1.2          pillar_1.2.1         foreign_0.8-69      
## [106] mixtools_1.1.0       survival_2.41-3      scatterplot3d_0.3-41
## [109] nnet_7.3-12          tsne_0.1-3           tibble_1.4.2        
## [112] KernSmooth_2.23-15   rmarkdown_1.9        grid_3.4.4          
## [115] data.table_1.10.4-3  FNN_1.1              ModelMetrics_1.1.0  
## [118] metap_0.8            digest_0.6.15        diptest_0.75-7      
## [121] numDeriv_2016.8-1    tidyr_0.8.0          R.utils_2.6.0       
## [124] stats4_3.4.4         munsell_0.4.3
```
