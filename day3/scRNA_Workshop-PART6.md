---
title: "single cell 10x single-cell analysis - part6"
author: "UC Davis Bioinformatics Core"
output:
  html_document:
    keep_md: true
---



#0. Setup
Load the final Seurat object, load libraries (also see additional required packages for each example)

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

#1. DE With Single Cell Data Using Limma
For differential expression using models more complex than those allowed by FindAllMarkers(), data from Seurat may be used in limma (https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf)

We illustrate by comparing sample 1 to sample 2 within cluster 0:

```r
library(limma)
cluster0 <- SubsetData(experiment.aggregate, ident.use = 0)
expr <- cluster0@data

# Filter out genes that are 0 for every cell in this cluster
bad <- which(rowSums(expr) == 0)
expr <- expr[-bad,]

mm <- model.matrix(~0 + orig.ident, data = cluster0@meta.data) 
fit <- lmFit(expr, mm)  
head(coef(fit)) # means in each sample for each gene
```

```
##         orig.identsample1 orig.identsample2 orig.identsample3
## Mrpl15        0.308056932       0.233220268       0.365672775
## Lypla1        0.296334196       0.264225566       0.325912665
## Tcea1         0.269425949       0.257152221       0.302355174
## Rgs20         0.006704523       0.007281172       0.011679030
## Atp6v1h       0.322142900       0.267904522       0.387039804
## Oprk1         0.001741248       0.001999901       0.001942202
```

```r
contr <- makeContrasts(orig.identsample2 - orig.identsample1, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contrasts = contr)
tmp <- eBayes(tmp)
topTable(tmp, sort.by = "P", n = 20) # top 20 DE genes
```

```
##                    logFC   AveExpr         t      P.Value    adj.P.Val        B
## Rps29          0.2323318 2.7382343  9.575604 2.176701e-21 2.482092e-17 37.68826
## Malat1         0.3985648 3.9071461  9.212669 6.141584e-20 3.501624e-16 34.40453
## 9130204L05Rik -0.2656994 0.4286414 -8.828500 1.851382e-18 7.037102e-15 31.05874
## Xist          -0.3786369 0.9744056 -8.533504 2.311746e-17 6.590210e-14 28.58098
## Snrpn         -0.2433459 2.2076308 -8.203731 3.539098e-16 8.071266e-13 25.90570
## Rpl41          0.1356469 3.6981405  7.849716 5.926983e-15 1.126423e-11 23.14557
## Ndn           -0.2352288 0.8125448 -7.606332 3.848212e-14 6.268737e-11 21.31550
## Fez1          -0.2333460 2.1579122 -7.483866 9.662257e-14 1.377234e-10 20.41555
## Calcb         -0.2422010 2.5339466 -7.342590 2.746832e-13 3.480236e-10 19.39478
## Oaz2          -0.2123888 1.0070645 -7.159490 1.035082e-12 1.180304e-09 18.09966
## Aldoa         -0.2216651 2.3634321 -6.768518 1.584710e-11 1.642768e-08 15.43981
## Zwint         -0.1978167 1.9348312 -6.745993 1.846446e-11 1.754585e-08 15.29097
## Pmm1          -0.2083454 1.3686732 -6.515075 8.611401e-11 7.171283e-08 13.79276
## Eif4a2        -0.1755853 2.3295583 -6.511694 8.804522e-11 7.171283e-08 13.77119
## Slc25a4       -0.1514549 3.0856518 -6.434632 1.455325e-10 1.106338e-07 13.28272
## Tuba1a        -0.1879961 3.1054677 -6.321886 3.005653e-10 2.142092e-07 12.57821
## Tuba1b        -0.2095709 1.7806285 -6.285291 3.793722e-10 2.394600e-07 12.35213
## Hspa8         -0.1943548 2.4676129 -6.280769 3.904116e-10 2.394600e-07 12.32429
## Ubb           -0.2208200 3.8414680 -6.277339 3.989949e-10 2.394600e-07 12.30318
## Rpl38          0.1756736 2.2450175  6.243086 4.954721e-10 2.824934e-07 12.09299
```
* logFC: log2 fold change (sample2/sample1)
* AveExpr: Average expression, in log2 counts per million, across all cells included in analysis (i.e. those in cluster 0)
* t: t-statistic, i.e. logFC divided by its standard error
* P.Value: Raw p-value from test that logFC differs from 0
* adj.P.Val: Benjamini-Hochberg false discovery rate adjusted p-value

The limma vignette linked above gives more detail on model specification.

# 2. Gene Ontology (GO) Enrichment of Genes Expressed in a Cluster

```r
library(topGO)
```

```
## Loading required package: BiocGenerics
```

```
## Loading required package: parallel
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ, clusterExport, clusterMap, parApply, parCapply, parLapply, parLapplyLB, parRapply, parSapply, parSapplyLB
```

```
## The following object is masked from 'package:limma':
## 
##     plotMA
```

```
## The following objects are masked from 'package:Matrix':
## 
##     as.vector, colMeans, colSums, rowMeans, rowSums, which
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, cbind, colMeans, colnames, colSums, do.call, duplicated, eval, evalq, Filter, Find, get, grep, grepl, intersect, is.unsorted, lapply, lengths, Map, mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank, rbind, Reduce, rowMeans, rownames, rowSums, sapply, setdiff, sort, table, tapply, union, unique, unsplit, which, which.max,
##     which.min
```

```
## Loading required package: graph
```

```
## Loading required package: Biobase
```

```
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with 'browseVignettes()'. To cite Bioconductor, see 'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```
## Loading required package: GO.db
```

```
## Loading required package: AnnotationDbi
```

```
## Loading required package: stats4
```

```
## Loading required package: IRanges
```

```
## Loading required package: S4Vectors
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following object is masked from 'package:Matrix':
## 
##     expand
```

```
## The following object is masked from 'package:base':
## 
##     expand.grid
```

```
## 
```

```
## Loading required package: SparseM
```

```
## 
## Attaching package: 'SparseM'
```

```
## The following object is masked from 'package:base':
## 
##     backsolve
```

```
## 
## groupGOTerms: 	GOBPTerm, GOMFTerm, GOCCTerm environments built.
```

```
## 
## Attaching package: 'topGO'
```

```
## The following object is masked from 'package:IRanges':
## 
##     members
```

```r
# install org.Mm.eg.db if not already installed (for mouse only)
# install.packages("org.Mm.eg.db") 
cluster0 <- SubsetData(experiment.aggregate, ident.use = 0)
expr <- cluster0@data
# Select genes that are expressed > 0 in at least 75% of cells (somewhat arbitrary definition)
n.gt.0 <- apply(expr, 1, function(x)length(which(x > 0)))
expressed.genes <- rownames(expr)[which(n.gt.0/ncol(expr) >= 0.75)]
all.genes <- rownames(expr)

# define geneList as 1 if gene is in expressed.genes, 0 otherwise
geneList <- ifelse(all.genes %in% expressed.genes, 1, 0)
names(geneList) <- all.genes

# Create topGOdata object
	GOdata <- new("topGOdata",
		ontology = "BP", # use biological process ontology
		allGenes = geneList,
		geneSelectionFun = function(x)(x == 1),
              annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "symbol")
```

```
## 
## Building most specific GOs .....
```

```
## Loading required package: org.Mm.eg.db
```

```
## 
```

```
## 	( 9613 GO terms found. )
```

```
## 
## Build GO DAG topology ..........
```

```
## 	( 13585 GO terms and 32452 relations. )
```

```
## 
## Annotating nodes ...............
```

```
## 	( 10752 genes annotated to the GO terms. )
```

```r
# Test for enrichment using Fisher's Exact Test
	resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
```

```
## 
## 			 -- Elim Algorithm -- 
## 
## 		 the algorithm is scoring 5401 nontrivial nodes
## 		 parameters: 
## 			 test statistic: fisher
## 			 cutOff: 0.01
```

```
## 
## 	 Level 19:	2 nodes to be scored	(0 eliminated genes)
```

```
## 
## 	 Level 18:	4 nodes to be scored	(0 eliminated genes)
```

```
## 
## 	 Level 17:	3 nodes to be scored	(8 eliminated genes)
```

```
## 
## 	 Level 16:	20 nodes to be scored	(8 eliminated genes)
```

```
## 
## 	 Level 15:	47 nodes to be scored	(10 eliminated genes)
```

```
## 
## 	 Level 14:	95 nodes to be scored	(55 eliminated genes)
```

```
## 
## 	 Level 13:	200 nodes to be scored	(91 eliminated genes)
```

```
## 
## 	 Level 12:	315 nodes to be scored	(580 eliminated genes)
```

```
## 
## 	 Level 11:	483 nodes to be scored	(645 eliminated genes)
```

```
## 
## 	 Level 10:	584 nodes to be scored	(1036 eliminated genes)
```

```
## 
## 	 Level 9:	712 nodes to be scored	(1417 eliminated genes)
```

```
## 
## 	 Level 8:	725 nodes to be scored	(1628 eliminated genes)
```

```
## 
## 	 Level 7:	794 nodes to be scored	(2010 eliminated genes)
```

```
## 
## 	 Level 6:	670 nodes to be scored	(2314 eliminated genes)
```

```
## 
## 	 Level 5:	437 nodes to be scored	(2862 eliminated genes)
```

```
## 
## 	 Level 4:	225 nodes to be scored	(3445 eliminated genes)
```

```
## 
## 	 Level 3:	64 nodes to be scored	(4692 eliminated genes)
```

```
## 
## 	 Level 2:	20 nodes to be scored	(4818 eliminated genes)
```

```
## 
## 	 Level 1:	1 nodes to be scored	(4818 eliminated genes)
```

```r
	GenTable(GOdata, Fisher = resultFisher, topNodes = 20, numChar = 60)
```

```
##         GO.ID                                                            Term Annotated Significant Expected  Fisher
## 1  GO:0006412                                                     translation       484          97    29.44 4.4e-18
## 2  GO:0015986                          ATP synthesis coupled proton transport        19          14     1.16 7.3e-14
## 3  GO:0000028                                ribosomal small subunit assembly        17          12     1.03 1.1e-11
## 4  GO:0006123        mitochondrial electron transport, cytochrome c to oxygen         8           7     0.49 2.3e-08
## 5  GO:0002181                                         cytoplasmic translation        47          16     2.86 7.4e-08
## 6  GO:0015991                         ATP hydrolysis coupled proton transport        18           8     1.09 4.6e-06
## 7  GO:0006120            mitochondrial electron transport, NADH to ubiquinone        13           6     0.79 5.9e-05
## 8  GO:0050821                                           protein stabilization       120          19     7.30 0.00011
## 9  GO:0034314                        Arp2/3 complex-mediated actin nucleation        27           8     1.64 0.00014
## 10 GO:0006122     mitochondrial electron transport, ubiquinol to cytochrome c        10           5     0.61 0.00016
## 11 GO:0046034                                           ATP metabolic process       156          47     9.49 0.00018
## 12 GO:1902255 positive regulation of intrinsic apoptotic signaling pathway...         6           4     0.36 0.00018
## 13 GO:0061077                              chaperone-mediated protein folding        51          11     3.10 0.00020
## 14 GO:0010592                   positive regulation of lamellipodium assembly        16           6     0.97 0.00023
## 15 GO:0006364                                                 rRNA processing       162          24     9.85 0.00026
## 16 GO:0045730                                               respiratory burst        11           5     0.67 0.00028
## 17 GO:0006900                                                membrane budding        47          10     2.86 0.00043
## 18 GO:0047497                       mitochondrion transport along microtubule        12           5     0.73 0.00045
## 19 GO:0006810                                                       transport      2768         254   168.37 0.00061
## 20 GO:0055114                                     oxidation-reduction process       583          69    35.46 0.00062
```
* Annotated: number of genes (out of all.genes) that are annotated with that GO term
* Significant: number of genes that are annotated with that GO term and meet our criteria for "expressed"
* Expected: Under random chance, number of genes that would be expected to be annotated with that GO term and meeting our criteria for "expressed"
* Fisher: (Raw) p-value from Fisher's Exact Test

#3. Weighted Gene Co-Expression Network Analysis (WGCNA)
WGCNA identifies groups of genes ("modules") with correlated expression.
WARNING: TAKES A LONG TIME TO RUN

```r
library(WGCNA)
```

```
## Loading required package: dynamicTreeCut
```

```
## Loading required package: fastcluster
```

```
## 
## Attaching package: 'fastcluster'
```

```
## The following object is masked from 'package:stats':
## 
##     hclust
```

```
## ==========================================================================
## *
## *  Package WGCNA 1.63 loaded.
## *
## *    Important note: It appears that your system supports multi-threading,
## *    but it is not enabled within WGCNA in R. 
## *    To allow multi-threading within WGCNA with all available cores, use 
## *
## *          allowWGCNAThreads()
## *
## *    within R. Use disableWGCNAThreads() to disable threading if necessary.
## *    Alternatively, set the following environment variable on your system:
## *
## *          ALLOW_WGCNA_THREADS=<number_of_processors>
## *
## *    for example 
## *
## *          ALLOW_WGCNA_THREADS=48
## *
## *    To set the environment variable in linux bash shell, type 
## *
## *           export ALLOW_WGCNA_THREADS=48
## *
## *     before running R. Other operating systems or shells will
## *     have a similar command to achieve the same aim.
## *
## ==========================================================================
```

```
## 
## Attaching package: 'WGCNA'
```

```
## The following object is masked from 'package:IRanges':
## 
##     cor
```

```
## The following object is masked from 'package:S4Vectors':
## 
##     cor
```

```
## The following object is masked from 'package:stats':
## 
##     cor
```

```r
options(stringsAsFactors = F)
datExpr <- t(experiment.aggregate@data)[,experiment.aggregate@var.genes]  # only use var.genes in analysis 

net <- blockwiseModules(datExpr, power = 10,
  corType = "bicor", # use robust correlation
	networkType = "signed", minModuleSize = 10,
	reassignThreshold = 0, mergeCutHeight = 0.15,
	numericLabels = F, pamRespectsDendro = FALSE,
	saveTOMs = TRUE,
	saveTOMFileBase = "TOM",
	verbose = 3)
```

```
##  Calculating module eigengenes block-wise from all genes
##    Flagging genes and samples with too many missing values...
##     ..step 1
##  ..Working on block 1 .
##     TOM calculation: adjacency..
##     ..will not use multithreading.
##      Fraction of slow calculations: 0.000000
##     ..connectivity..
##     ..matrix multiplication (system BLAS)..
##     ..normalization..
##     ..done.
##    ..saving TOM for block 1 into file TOM-block.1.RData
##  ....clustering..
##  ....detecting modules..
##  ....calculating module eigengenes..
##  ....checking kME in modules..
```

```
## Warning in bicor(structure(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.60022198957559, : bicor: zero MAD in variable 'x'. Pearson correlation was used for individual columns with zero (or missing) MAD.
```

```
## Warning in bicor(structure(c(1.63056832072081, 0, 0, 0, 0, 0, 0, 0, 1.74033425742062, : bicor: zero MAD in variable 'x'. Pearson correlation was used for individual columns with zero (or missing) MAD.
```

```
## Warning in bicor(structure(c(0, 1.5630975751754, 1.51509992585345, 0, 0, : bicor: zero MAD in variable 'x'. Pearson correlation was used for individual columns with zero (or missing) MAD.
```

```
## Warning in bicor(structure(c(0, 0, 0, 1.32405205224267, 0, 0, 2.00065018203098, : bicor: zero MAD in variable 'x'. Pearson correlation was used for individual columns with zero (or missing) MAD.
```

```
## Warning in bicor(structure(c(1.63056832072081, 2.77846763393956, 0, 0, 0, : bicor: zero MAD in variable 'x'. Pearson correlation was used for individual columns with zero (or missing) MAD.
```

```
##      ..removing 1 genes from module 5 because their KME is too low.
```

```
## Warning in bicor(structure(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.688381484486931, : bicor: zero MAD in variable 'x'. Pearson correlation was used for individual columns with zero (or missing) MAD.
```

```
## Warning in bicor(structure(c(1.63056832072081, 0, 0, 2.48777609221595, 2.31903926177167, : bicor: zero MAD in variable 'x'. Pearson correlation was used for individual columns with zero (or missing) MAD.
```

```
##      ..removing 1 genes from module 7 because their KME is too low.
```

```
## Warning in bicor(structure(c(3.24417925903256, 3.55426931273324, 3.10455324858405, : bicor: zero MAD in variable 'x'. Pearson correlation was used for individual columns with zero (or missing) MAD.
```

```
##      ..removing 4 genes from module 8 because their KME is too low.
```

```
## Warning in bicor(structure(c(1.63056832072081, 0, 0, 1.32405205224267, 0, : bicor: zero MAD in variable 'x'. Pearson correlation was used for individual columns with zero (or missing) MAD.
```

```
## Warning in bicor(structure(c(0, 1.5630975751754, 1.51509992585345, 1.32405205224267, : bicor: zero MAD in variable 'x'. Pearson correlation was used for individual columns with zero (or missing) MAD.
```

```
## Warning in (function (x, y = NULL, robustX = TRUE, robustY = TRUE, use = "all.obs", : bicor: zero MAD in variable 'x'. Pearson correlation was used for individual columns with zero (or missing) MAD.
```

```
##  ..merging modules that are too close..
##      mergeCloseModules: Merging modules whose distance is less than 0.15
##        Calculating new MEs...
```

```r
table(net$colors)
```

```
## 
##      blue     brown     green      grey turquoise    yellow 
##        39        26        17       974        46        21
```

```r
# Convert labels to colors for plotting
mergedColors = net$colors
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
"Module colors",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
```

![](scRNA_Workshop-PART6_files/figure-html/unnamed-chunk-4-1.png)<!-- -->
Genes in grey module are unclustered.

What genes are in the "blue" module?

```r
colnames(datExpr)[net$colors == "blue"]
```

```
##  [1] "Cd55"     "Smyd3"    "Cd82"     "Cd44"     "S100a6"   "Lpar3"    "Ugcg"     "Sh3bgrl3" "Pop5"     "Tmem233"  "Arpc1b"   "Cav1"     "Gm765"    "Cd9"      "Emp3"     "Mrgpra3"  "Nmb"      "Cd24a"    "Plxnc1"   "Dusp26"   "Adk"      "Prkcd"    "Lgals3"   "Cadm1"    "Scn11a"   "Tmem158"  "Kcnmb1"   "Skp1a"    "Atox1"    "Prkca"    "Ly86"     "Prkar2b"  "Crip2"    "Mal2"     "Carhsp1"  "Cystm1"   "Ctxn3"    "Rab27b"   "Gna14"
```

Each cluster is represented by a summary "eigengene".
Plot eigengenes for each non-grey module by clusters from Seurat:

```r
f <- function(module){
  eigengene <- unlist(net$MEs[paste0("ME", module)])
  means <- tapply(eigengene, experiment.aggregate@ident, mean, na.rm = T)
  return(means)
}
modules <- c("blue", "brown", "green", "turquoise", "yellow")
plotdat <- sapply(modules, f)
matplot(plotdat, col = modules, type = "l", lwd = 2, xaxt = "n", xlab = "Seurat Cluster",
        ylab = "WGCNA Module Eigengene")
axis(1, at = 1:16, labels = 0:15)
matpoints(plotdat, col = modules, pch = 21)
```

![](scRNA_Workshop-PART6_files/figure-html/unnamed-chunk-6-1.png)<!-- -->


```r
sessionInfo()
```

```
## R version 3.4.4 (2018-03-15)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 14.04.3 LTS
## 
## Matrix products: default
## BLAS: /usr/lib/openblas-base/libblas.so.3
## LAPACK: /usr/lib/lapack/liblapack.so.3.0
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] WGCNA_1.63            fastcluster_1.1.24    dynamicTreeCut_1.63-1 org.Mm.eg.db_3.4.1    topGO_2.28.0          SparseM_1.77          GO.db_3.4.1           AnnotationDbi_1.38.2  IRanges_2.10.5        S4Vectors_0.14.7      Biobase_2.36.2        graph_1.54.0          BiocGenerics_0.22.1   limma_3.32.10         Seurat_2.3.0          Matrix_1.2-3          cowplot_0.8.0         ggplot2_2.2.1        
## 
## loaded via a namespace (and not attached):
##   [1] snow_0.4-2            backports_1.1.1       Hmisc_4.0-3           VGAM_1.0-4            sn_1.5-0              plyr_1.8.4            igraph_1.1.2          lazyeval_0.2.1        splines_3.4.4         robust_0.4-18         digest_0.6.12         foreach_1.4.3         htmltools_0.3.6       lars_1.2              gdata_2.18.0          magrittr_1.5          checkmate_1.8.5       memoise_1.1.0         fit.models_0.5-14     doParallel_1.0.11    
##  [21] cluster_2.0.7         mixtools_1.1.0        ROCR_1.0-7            sfsmisc_1.1-1         recipes_0.1.0         gower_0.1.2           matrixStats_0.52.2    dimRed_0.1.0          R.utils_2.5.0         colorspace_1.3-2      rrcov_1.4-3           dplyr_0.7.4           impute_1.50.1         bindr_0.1             survival_2.41-3       zoo_1.8-1             iterators_1.0.8       ape_5.0               glue_1.2.0            DRR_0.0.2            
##  [41] gtable_0.2.0          ipred_0.9-6           kernlab_0.9-25        ddalpha_1.3.1         prabclus_2.2-6        DEoptimR_1.0-8        scales_0.5.0          mvtnorm_1.0-6         DBI_0.7               Rcpp_0.12.13          metap_0.9             dtw_1.18-1            htmlTable_1.9         tclust_1.3-1          foreign_0.8-66        proxy_0.4-19          mclust_5.3            preprocessCore_1.38.1 SDMTools_1.1-221      Formula_1.2-2        
##  [61] tsne_0.1-3            lava_1.5.1            prodlim_1.6.1         htmlwidgets_0.9       FNN_1.1               gplots_3.0.1          RColorBrewer_1.1-2    fpc_2.1-10            acepack_1.4.1         modeltools_0.2-21     ica_1.0-1             pkgconfig_2.0.1       R.methodsS3_1.7.1     flexmix_2.3-14        nnet_7.3-12           caret_6.0-77          rlang_0.1.4           reshape2_1.4.2        munsell_0.4.3         tools_3.4.4          
##  [81] RSQLite_1.1-2         ranger_0.8.0          ggridges_0.4.1        evaluate_0.10.1       stringr_1.2.0         yaml_2.1.14           ModelMetrics_1.1.0    knitr_1.20            fitdistrplus_1.0-9    robustbase_0.92-8     caTools_1.17.1        purrr_0.2.4           RANN_2.5.1            bindrcpp_0.2          packrat_0.4.8-1       pbapply_1.3-3         nlme_3.1-137          R.oo_1.21.0           RcppRoll_0.2.2        compiler_3.4.4       
## [101] png_0.1-7             tibble_1.3.4          pcaPP_1.9-73          stringi_1.1.5         lattice_0.20-33       trimcluster_0.1-2     diffusionMap_1.1-0    lmtest_0.9-36         data.table_1.10.4-3   bitops_1.0-6          irlba_2.3.1           R6_2.2.2              latticeExtra_0.6-28   KernSmooth_2.23-15    gridExtra_2.3         codetools_0.2-14      MASS_7.3-44           gtools_3.5.0          assertthat_0.2.0      CVST_0.2-1           
## [121] rprojroot_1.2         withr_2.1.0           mnormt_1.5-5          diptest_0.75-7        doSNOW_1.0.16         grid_3.4.4            rpart_4.1-10          timeDate_3012.100     tidyr_0.7.2           class_7.3-14          rmarkdown_1.8         segmented_0.5-2.2     Rtsne_0.13            numDeriv_2016.8-1     scatterplot3d_0.3-40  lubridate_1.7.1       base64enc_0.1-3
```




