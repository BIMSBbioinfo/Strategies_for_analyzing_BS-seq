The comparison of differentially methylated bases detection tools on simulated data
================
Katarzyna Wreczycka
2017-07-26

<!-- rnb-text-begin -->
<!-- rnb-text-end -->
<!-- rnb-text-begin -->
Goal
====

Here, we examined the performance of various differential methylation methods. We compared three classes of methods:

1.  t-test/linear regression,
2.  logistic regression,
3.  beta binomial regression.

Three different tools were used:

1.  DSS (beta binomial regression),
2.  limma (linear regression),
3.  methylKit (logistic regression with/without overdispersion correction).

Functions
=========

Here are functions to run limma and to run our simulation:

<!-- rnb-text-end -->
``` r
# Load libraries and functions
library("methylKit")
library(ggplot2)
library(gridExtra)
library(grid)
source("./functions/dataSim2.R")
source("./functions/limma.R")
```

<!-- rnb-text-begin -->
<!-- rnb-text-end -->
``` r

#' Calculate rates of models compared to simulation
#' 
#' @param simOutput this is the output of dataSim2 function
#' @param sub.methylDiff this is the q-value filtered methylDiff object
#'                       output of getMethylDiff()
#' @param methylDiff this is the non-filtered methylDiff object
#'                   output of calculateDiffMeth() and similar functions
#' @return returns a vector of accuracy metrics, TP, FP, Sensivity, etc
calc.rates<-function(simOutput,
                         sub.methylDiff,
                         methylDiff # this argument is not needed: TODO remove it to not be confused again
                         ){
  all=paste(simOutput[[1]][[1]],simOutput[[1]][[2]],
            simOutput[[1]][[3]])
  
  true.dm=all[simOutput[[2]]]
  true.ndm=all[-simOutput[[2]]]
  
  pred.dm=paste(sub.methylDiff[[1]],sub.methylDiff[[2]],
                sub.methylDiff[[3]])
  pred.ndm=all[! all %in% pred.dm]
  
  TP=sum(true.dm %in% pred.dm)
  
  FN=sum(pred.ndm %in% true.dm)
  
  FP=sum(pred.dm %in% true.ndm)
  
  TN=sum(pred.ndm %in% true.ndm)
  
  p = TP / (TP + FP)
  r = TP / (TP+FN)
  f_score = 2*((p*r)/(p+r))
  
  return(c(TP=TP,FN=FN,FP=FP,TN=TN,
           acc=(TP+TN)/length(all),
           spec=TN/(TN+FP) ,
           sens=TP/(TP+FN),
           f_score= f_score,
           precision=as.numeric(TP / (TP + FP)),
           recall=r,
           NPV= as.numeric(TN / (TN + FN))
           ) )
}


# ---------------------------------------------------------------------------- #
#' Run simulation
#' 
#' Call differentially methylated cytosines using methylKit, DSS and limma.
#' It calculate true positive positives (TP), false negatives (FN), false positives (FP),
#' accuracy (acc), specificity (spec), sensiticity (sens) and F-score (f_score).
#' 
#' @param sim.methylBase a methylBase object from the methylKit library
#' @param cores a number of cores
#' @param difference cutoff for absolute value of methylation percentage change
#'                   between test and control (default:5)
#' @param qvalue cutoff for qvalue of differential methylation statistic
#'               (default:0.01)
#' @return returns a matrix with TP, FN, FP, TN, acc, spec, sens, f_score (columns)
#'         using tools that calculate differentially methylated regions (rows)
run.models = function(sim.methylBase, cores=1,
                      difference=5, qvalue=0.01){
  
  require(methylKit)
  require(DSS)

  
  ## run methylkit
  combined = data.frame(test=c("F", "Chisq","F", "Chisq"),
                        adjust="qvalue",
                        overd=c("none","none", "MN", "MN"),
                        name=c("methylKit.F.qvalue.none",
                               "methylKit.Chisq.qvalue.none",
                               "methylKit.F.qvalue.MN",
                               "methylKit.Chisq.qvalue.MN"), 
                        stringsAsFactors = FALSE)
  diff.list = list()
  methylKit.list=list()
  for(i in 1:nrow(combined)){
    co = combined[i,]
    methylkit.obj <- calculateDiffMeth(sim.methylBase[[1]], 
                                       overdispersion=co$overd,
                                       adjust = co$adjust,
                                       test=co$test,
                                       mc.cores=cores)
    methylkit.obj.diff = getMethylDiff(methylkit.obj, 
                                       difference=difference,qvalue=qvalue)
    diff.list[[i]] <- methylkit.obj.diff
    methylKit.list[[i]]=calc.rates(sim.methylBase,
                                       methylkit.obj.diff,
                                       methylkit.obj)
    
  }
  names(methylKit.list) <- combined$name
  names(diff.list) <- combined$name
  
  
  ## run DSS
  dss.qvalue = calculateDiffMethDSS(sim.methylBase[[1]],
                                    adjust="qvalue",
                                    mc.cores=cores)
  dss.qvalue.diff = getMethylDiff(dss.qvalue, difference=difference,qvalue=qvalue)
  
  diff.list[["DSS"]]=dss.qvalue.diff
  methylKit.list[["DSS"]]=calc.rates(sim.methylBase,dss.qvalue.diff,
                                                dss.qvalue)

  limma.qvalue=limma.meth(sim.methylBase[[1]])
  limma.qvalue.diff = getMethylDiff(limma.qvalue, 
                                    difference=difference,qvalue=qvalue)
  
  
  ## run limma
  diff.list[["limma.qvalue"]] = limma.qvalue.diff
  methylKit.list[["limma.qvalue"]]=calc.rates(sim.methylBase,
                                                  limma.qvalue.diff,
                                                  limma.qvalue)

  do.call("rbind",methylKit.list)
}
```

<!-- rnb-text-begin -->
Simulation
==========

We simulated a dataset consisting of 6 samples (3 controls and 3 samples with treatment). The read coverage modeled by a binomial distribution. The methylation background followed a beta distribution with parameters \(alpha=0.4\), \(beta=0.5\) and \(theta=10\). We simulated 6 sets of 5000 CpG sites where methylation at 50% of the sites was affected by the treatment to varying degrees - specifically, methylation was elevated by 5%, 10%, 15%, 20% and 25% in the test sample respectively in each set.

To adjust p-values for multiple testing, we used q-value method and we defined differentially methylated CpG sites with q-values below 0.01 for all examined methods. We calculated sensitivity, specificity and F-score for each of the three methods above. Sensitivity measured the proportion of true differentially methylated CpGs that were correctly identified as such, specificity was calculated as the proportion of detected CpGs that were truly not differentially methylated and correctly identified as such and F-score refers to a way to measure sensitivity and specificity by calculating their harmonic mean.

Here, we calculate sensitivity, specificity and F-score of performance of tools for calling differentially methylated cytosines:

<!-- rnb-text-end -->
``` r

# variables
effects = c(5, 10, 15, 20, 25)
cores=20

models.res=list()
set.seed(111)
for(effect in effects){
  
  # Effect by the treatment
  print(effect)

  # Generate simulated data using methylKit library
  sim.methylBase = dataSim2(replicates=6,
                                 sites=5000,
                                 treatment=c(1,1,1,0,0,0),
                                 percentage=50,
                                 effect=effect,
                                 add.info=TRUE)
  
  # Run models 
  models.res[[as.character(effect)]] = run.models(sim.methylBase, cores=cores,
                                                   difference=5, qvalue=0.01)
  
}
## [1] 5
## Loading required package: DSS
## Warning in library(package, lib.loc = lib.loc, character.only = TRUE,
## logical.return = TRUE, : there is no package called 'DSS'
## two groups detected:
##  will calculate methylation difference as the difference of
## treatment (group: 1) - control (group: 0)
## two groups detected:
##  will calculate methylation difference as the difference of
## treatment (group: 1) - control (group: 0)
## two groups detected:
##  will calculate methylation difference as the difference of
## treatment (group: 1) - control (group: 0)
## two groups detected:
##  will calculate methylation difference as the difference of
## treatment (group: 1) - control (group: 0)
## Using internal DSS code... 
## [1] 10
## Loading required package: DSS
## Warning in library(package, lib.loc = lib.loc, character.only = TRUE,
## logical.return = TRUE, : there is no package called 'DSS'
## two groups detected:
##  will calculate methylation difference as the difference of
## treatment (group: 1) - control (group: 0)
## two groups detected:
##  will calculate methylation difference as the difference of
## treatment (group: 1) - control (group: 0)
## two groups detected:
##  will calculate methylation difference as the difference of
## treatment (group: 1) - control (group: 0)
## two groups detected:
##  will calculate methylation difference as the difference of
## treatment (group: 1) - control (group: 0)
## Using internal DSS code... 
## [1] 15
## Loading required package: DSS
## Warning in library(package, lib.loc = lib.loc, character.only = TRUE,
## logical.return = TRUE, : there is no package called 'DSS'
## two groups detected:
##  will calculate methylation difference as the difference of
## treatment (group: 1) - control (group: 0)
## two groups detected:
##  will calculate methylation difference as the difference of
## treatment (group: 1) - control (group: 0)
## two groups detected:
##  will calculate methylation difference as the difference of
## treatment (group: 1) - control (group: 0)
## two groups detected:
##  will calculate methylation difference as the difference of
## treatment (group: 1) - control (group: 0)
## Using internal DSS code... 
## [1] 20
## Loading required package: DSS
## Warning in library(package, lib.loc = lib.loc, character.only = TRUE,
## logical.return = TRUE, : there is no package called 'DSS'
## two groups detected:
##  will calculate methylation difference as the difference of
## treatment (group: 1) - control (group: 0)
## two groups detected:
##  will calculate methylation difference as the difference of
## treatment (group: 1) - control (group: 0)
## two groups detected:
##  will calculate methylation difference as the difference of
## treatment (group: 1) - control (group: 0)
## two groups detected:
##  will calculate methylation difference as the difference of
## treatment (group: 1) - control (group: 0)
## Using internal DSS code... 
## [1] 25
## Loading required package: DSS
## Warning in library(package, lib.loc = lib.loc, character.only = TRUE,
## logical.return = TRUE, : there is no package called 'DSS'
## two groups detected:
##  will calculate methylation difference as the difference of
## treatment (group: 1) - control (group: 0)
## two groups detected:
##  will calculate methylation difference as the difference of
## treatment (group: 1) - control (group: 0)
## two groups detected:
##  will calculate methylation difference as the difference of
## treatment (group: 1) - control (group: 0)
## two groups detected:
##  will calculate methylation difference as the difference of
## treatment (group: 1) - control (group: 0)
## Using internal DSS code...
```

<!-- rnb-text-begin -->
Visualisation
=============

The visualisation of sensitivity, specificity and F-score for different effect sizes.

<!-- rnb-text-end -->
``` r

# Convert list of matrices to a data.frame
models.res.ma = do.call("rbind", models.res)
models.res.df = data.frame(models.res.ma)
models.res.df = cbind(models.res.df, 
                      tool=rownames(models.res.ma),
                      effect=as.factor(as.numeric(sapply(effects, 
                                               function(x) 
                                                 rep(x, nrow(models.res[[1]])  )))))

# for the publication
models.res.df$tool = as.character(models.res.df$tool)
models.res.df = models.res.df[-which(models.res.df$tool=="methylKit.F.qvalue.none"),]
models.res.df$tool = as.factor(models.res.df$tool)

# Rename names of tools
levels(models.res.df$tool)[levels(models.res.df$tool)=="DSS.qvalue"] <- "DSS"
levels(models.res.df$tool)[levels(models.res.df$tool)=="limma.qvalue"] <- "limma"
levels(models.res.df$tool)[levels(models.res.df$tool)=="methylKit.Chisq.qvalue.MN"] <- "methylKit-Chisqtest-OC"
levels(models.res.df$tool)[levels(models.res.df$tool)=="methylKit.Chisq.qvalue.none"] <- "methylKit-Chisqtest"
levels(models.res.df$tool)[levels(models.res.df$tool)=="methylKit.F.qvalue.MN"] <- "methylKit-Ftest-OC"
# Sort names of tools
#levels(models.res.df$tool)[levels(models.res.df$tool)=="methylKit.F.qvalue.none"] <- "methylKit-Ftest"
#models.res.df$tool=relevel(models.res.df$tool, "methylKit-Ftest")
models.res.df$tool=relevel(models.res.df$tool, "methylKit-Ftest-OC")
models.res.df$tool=relevel(models.res.df$tool, "methylKit-Chisqtest-OC")
models.res.df$tool=relevel(models.res.df$tool, "methylKit-Chisqtest")
models.res.df$tool=relevel(models.res.df$tool, "limma")
models.res.df$tool=relevel(models.res.df$tool, "DSS")

  

# A palette
cbPalette <- c( "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

p_sens<-ggplot(models.res.df,aes(effect, sens, fill=tool))+
  geom_bar(stat="identity",position='dodge',colour="black",
           width=0.65)+
  coord_cartesian(ylim=c(0.0,1.00))+
  scale_fill_manual(values=cbPalette)+
  labs(y="Sensitivity", x="Effect size (methylation difference)",fill='Tool')

p_spec<-ggplot(models.res.df,aes(effect, spec, fill=tool))+
  geom_bar(stat="identity",
           position="dodge",
           colour="black",
           width=0.65)+
  coord_cartesian(ylim=c(0.0,1.00))+
  scale_fill_manual(values=cbPalette)+
  labs(y="Specificity", x="Effect size (methylation difference)",fill='Tool')

p_recall <- ggplot(models.res.df,aes(effect, recall, fill=tool))+
  geom_bar(stat="identity",position='dodge',colour="black",
           width=0.65)+
  coord_cartesian(ylim=c(0.0,1.00))+
  scale_fill_manual(values=cbPalette)+
  labs(y="Recall", x="Effect size (methylation difference)",fill='Tool')

p_fscore <- ggplot(models.res.df,aes(effect, f_score, fill=tool))+
  geom_bar(stat="identity",position='dodge',colour="black",
           width=0.65)+
  coord_cartesian(ylim=c(0.0,1.00))+
  scale_fill_manual(values=cbPalette)+
  labs(y="F-score", x="Effect size (methylation difference)",fill='Tool')


p_sens1 <- arrangeGrob(p_sens, top = textGrob("a", x=unit(0, "npc"),y=unit(1, "npc"),
                                              just=c("left","top"), gp=gpar(col="black", fontsize=14)))
p_spec1 <- arrangeGrob(p_spec, top = textGrob("b", x=unit(0, "npc"),y=unit(1, "npc"),just=c("left","top"), gp=gpar(col="black", fontsize=14)))
p_fscore1 <- arrangeGrob(p_fscore, top = textGrob("c", x=unit(0, "npc"),y=unit(1, "npc"),just=c("left","top"), gp=gpar(col="black", fontsize=14)))


grid.arrange(p_sens1, p_spec1, p_fscore1, ncol = 2, nrow=2)
```

<img src="diff_meth_figs/diff_meth-vis-1.png" width="960" />

<!-- rnb-text-begin -->
<!-- rnb-text-end -->
``` r
sessionInfo()
## R version 3.4.0 (2017-04-21)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 16.04.1 LTS
## 
## Matrix products: default
## BLAS: /usr/lib/libblas/libblas.so.3.6.0
## LAPACK: /usr/lib/lapack/liblapack.so.3.6.0
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=de_DE.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=de_DE.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=de_DE.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
##  [1] grid      parallel  stats4    stats     graphics  grDevices utils    
##  [8] datasets  methods   base     
## 
## other attached packages:
##  [1] gridExtra_2.2.1      ggplot2_2.2.1        qvalue_2.8.0        
##  [4] limma_3.32.3         emdbook_1.3.9        methylKit_1.3.3     
##  [7] GenomicRanges_1.28.4 GenomeInfoDb_1.12.2  IRanges_2.10.2      
## [10] S4Vectors_0.14.3     BiocGenerics_0.22.0  rmarkdown_1.6       
## 
## loaded via a namespace (and not attached):
##  [1] SummarizedExperiment_1.6.3 gtools_3.5.0              
##  [3] reshape2_1.4.2             splines_3.4.0             
##  [5] lattice_0.20-34            tcltk_3.4.0               
##  [7] colorspace_1.3-2           htmltools_0.3.6           
##  [9] rtracklayer_1.36.4         yaml_2.1.14               
## [11] base64enc_0.1-3            XML_3.98-1.9              
## [13] rlang_0.1.1                R.oo_1.21.0               
## [15] R.utils_2.5.0              BiocParallel_1.10.1       
## [17] fastseg_1.22.0             matrixStats_0.52.2        
## [19] GenomeInfoDbData_0.99.0    plyr_1.8.4                
## [21] stringr_1.2.0              zlibbioc_1.22.0           
## [23] Biostrings_2.44.1          munsell_0.4.3             
## [25] gtable_0.2.0               R.methodsS3_1.7.1         
## [27] coda_0.19-1                evaluate_0.10.1           
## [29] labeling_0.3               Biobase_2.36.2            
## [31] knitr_1.16                 Rcpp_0.12.12              
## [33] scales_0.4.1               backports_1.1.0           
## [35] DelayedArray_0.2.7         jsonlite_1.5              
## [37] XVector_0.16.0             Rsamtools_1.28.0          
## [39] digest_0.6.12              stringi_1.1.5             
## [41] numDeriv_2016.8-1          rprojroot_1.2             
## [43] tools_3.4.0                bitops_1.0-6              
## [45] bbmle_1.0.19               magrittr_1.5              
## [47] lazyeval_0.2.0             RCurl_1.95-4.8            
## [49] tibble_1.3.3               MASS_7.3-45               
## [51] Matrix_1.2-7.1             data.table_1.10.4         
## [53] mclust_5.3                 GenomicAlignments_1.12.1  
## [55] compiler_3.4.0
```

<!-- rnb-text-begin -->
<!-- rnb-text-end -->
