# ---
# title: "BSmooth"
# author: "Katarzyna Wreczycka"
# date: "1/17/2017"
# output: html_document
# ---
  
.libPaths(c("/home/kwreczy/Rlibs/3.4/"))
library("methylKit")
library("bsseq") #BSmooth


## Tutorial


a = readRDS("/data/akalin/Projects/AAkalin_Catalog_RI/Results/ENSG00000263874_Benchmarki1.overlap.rds")

## Simulation
# sim.methylBase = get.simdata(25)
# simobj=sim.methylBase$all
# diff=sim.methylBase$diff
# no.diff=sim.methylBase$nodiff
# print(simobj)
# 
# oobj = getData(sim.methylBase[[1]])
# Covs <- as.matrix(sim.d[, sim.methylBase[[1]]@coverage.index ])
# Meths <- as.matrix(sim.d[, sim.methylBase[[1]]@numCs.index ])
# nbr.st = length(simobj@treatment)/2
# colData = data.frame(Type=simobj@treatment,
#            Pair=c(paste0("sample", 1:nbr.st),paste0("sample", 1:nbr.st)))
# rownames(colData) <- c(paste0(rep("T",nbr.st), 1:nbr.st), paste0(rep("NONT", nbr.st), 1:nbr.st))
# print(colData)
# BS1 <- BSseq(chr = oobj$chr, pos = oobj$start,
#                M = Meths, Cov = Covs, sampleNames = rownames(colData),
#              pData=colData)
# BS1.fit <- BSmooth(BS1, mc.cores = 30, verbose = TRUE)
# 
# 
# # remove CpGs with little or no coverage
# # If this is not done, you
# # may find many DMRs in areas of the genome with very little coverage, which are most likely false
# # positives
# BS.cov <- getCoverage(BS1.fit)
# # keep CpGs where at least 2 cancer samples and at least 2 normal samples have at least 2x in
# # coverage
# keepLoci.ex <- which(rowSums(BS.cov[, BS1$Type == 1] >= 2) >= 2 &
#                     rowSums(BS.cov[, BS1$Type == 0] >= 2) >= 2)
# length(keepLoci.ex)
# BS1.fit <- BS1.fit[keepLoci.ex,]
# 
# # compute t-statistic
# # T-statistics are formed as the difference in means between group 1 and group 2 divided by an estimate of the standard deviation, assuming that the variance in the two groups are the same (same), that we have paired samples (paired) or only estimate the variance based on group 2 (group2)
# BS.tstat <- BSmooth.tstat(BS1.fit,
#                                     group1 = paste0(rep("T",nbr.st), 1:nbr.st),
#                                     group2 = paste0(rep("NONT",nbr.st), 1:nbr.st),
#                                     estimate.var = "same",
#                                     local.correct = FALSE,
#                                     verbose = TRUE)
# 
# plot(BS.tstat)
# 
# 
# # compute differentially methylated regions (DMRs) by
# #  thresholding the t-statistics
# 
# dmrs0 <- dmrFinder(BS.tstat, cutoff = NULL)

#################################


convert.methylBase2BSseq.obj = function(methylBase.obj){
  require(methylKit)
  require(bsseq)
  
  # assumtion: there is the same number of treated and non treated samples
  nbr.samples = length(methylBase.obj@treatment)
  n.t = nbr.samples/2
  
  methylBase.d = getData(methylBase.obj)
  Covs <- as.matrix(methylBase.d[, methylBase.obj@coverage.index ])
  Meths <- as.matrix(methylBase.d[, methylBase.obj@numCs.index ])
  colData = data.frame(Type=methylBase.obj@treatment,
                 Pair=c(paste0("sample", 1:n.t),paste0("sample", 1:n.t)))
  rownames(colData) <- c(paste0(rep("C",n.t), 1:n.t), paste0(rep("T", n.t), 1:n.t))
  BS.obj <- BSseq(chr = methylBase.d$chr, pos = methylBase.d$start,
               M = Meths, Cov = Covs, sampleNames = rownames(colData),
               pData=colData)
  return(BS.obj)
}

runBSmooth.smooth = function(methylBase.obj, 
                             estimate.var="group2",
                             cores=10){
  
  nmbr.treat = sum(methylBase.obj@treatment==1)
  nmbr.notreat = sum(methylBase.obj@treatment==0)
  BSseq.obj = convert.methylBase2BSseq.obj(methylBase.obj)
  
  BS1.fit <- BSmooth(BSseq.obj, 
                     mc.cores = cores, 
                     verbose = TRUE)
  BS1.tstat <- try(BSmooth.tstat(BS1.fit,
                                 group1 = paste0(rep("T", nmbr.treat), 1:nmbr.treat),
                                 group2 = paste0(rep("C",nmbr.notreat), 1:nmbr.notreat),
                                 estimate.var = estimate.var,
                                 local.correct = TRUE,
                                 verbose = TRUE,
                                 mc.cores=cores))
  if(class( BS1.tstat)=="try-error"){
    print(BS1.tstat)
    return(NA)
  }
  return(BS1.tstat)
}



runBSmooth = function(methylBase.obj, 
                      estimate.var="group2",
                      cutoff=c(-4.5, 4.5),
                      qcutoff=NULL,
                      cores=10){
  
    nmbr.treat = sum(methylBase.obj@treatment==1)
    nmbr.notreat = sum(methylBase.obj@treatment==0)
    BSseq.obj = convert.methylBase2BSseq.obj(methylBase.obj)
    
    BS1.fit <- BSmooth(BSseq.obj, 
                   mc.cores = cores, 
                   verbose = TRUE)
    
    # row.has.na <- apply(BS1.fit@assays[[3]], 1, function(x){any(is.na(x))})
    # wh.row.na = which(row.has.na)
    # BS1.fit@assays = BS1.fit@assays[-wh.row.na,]
    # BS1.fit@assays = BS1.fit@assays[ - which(BS1.fit@assays[[1]][,1]==0 &     BS1.fit@assays[[1]][,2]==0), ]
    # 
    #It is open to personal preferences exactly which CpGs to remove, but for this analysis we
    #will only keep CpGs where at least 2 treatment samples and at least 2 control samples have at least 2x in
    #coverage.
    #samples. (2) Within each CpG island we considered regions to re-
    #ceive a DMR only if each of its CpG sites was covered in at least
    #half of all samples. Those regions are referred to as ‘covered is-
    # BS.cov <- getCoverage(BS1.fit)
    # keepLoci.ex <- which(rowSums(BS.cov[, BS1.fit$Type == "1"] >= 2) >= 2 &
    #                        rowSums(BS.cov[, BS1.fit$Type == "0"] >= 2) >= 2)
    # #length(keepLoci.ex)
    # BS1.fit <- BS1.fit[keepLoci.ex,]
    #group1: A vector of sample names or indexes for the ‘treatment’ group.
    #group2: A vector of sample names or indexes for the ‘control’ group.
    BS1.tstat <- try(BSmooth.tstat(BS1.fit,
                           group1 = paste0(rep("T", nmbr.treat), 1:nmbr.treat),
                           group2 = paste0(rep("C",nmbr.notreat), 1:nmbr.notreat),
                            estimate.var = estimate.var,
                           local.correct = TRUE,
                           verbose = TRUE))
    if(class( BS1.tstat)=="try-error"){
      print(BS1.tstat)
      return(NA)
    }
    #We rank DMRs by the column areaStat which is the sum of the t-statistics in each CpG. This is
    #kind of the area of the DMR, except that it is weighted by the number of CpGs and not by genomic
    #length. This is currently the best statistic we know, although it is far from perfect (we would like to do
    #something better).
    # ! Essentially identifies regions of the genome where all methylation loci have 
    # an associated t-statistic that is beyond a (low, high) cutoff.
    dmrs0 <- try(dmrFinder(BS1.tstat, 
                   cutoff = cutoff,
                   qcutoff = qcutoff)) #cutoff using these quantiles of the t-statistic.
    # Error in quantile.default(dmrStat, qcutoff) : 
    #   missing values and NaN's not allowed if 'na.rm' is FALSE
    if(class( dmrs0)=="try-error"){
      print(dmrs0)
      return(NA)
    }
    if(is.null(dmrs0)){
      return(GRanges())
    }
    
    # Here, we filter out DMRs that do not have at least 3 CpGs in them and at least a mean difference
    # (across the DMR) in methylation between normal and cancers of at least 0.1
    # dmrs <- subset(dmrs0, n >= 3 & meanDiff >= 0.1)
    # #dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
    # 
    # #mean methylation difference below 0.1 and contain-
    #  # ing fewer than 3 CpGs
    # if(nrow(dmrs)==0){
    #  return(GRanges())
    # }
    dmrs.gr = makeGRangesFromDataFrame(dmrs0,
                                           keep.extra.columns=TRUE)
    return(dmrs.gr)
}

# 
# 
# tmp=get.simdata(25)
# sim.methylBase = tmp[[1]]
# sim.diff= tmp$diff
# sim.effect= tmp$treatment.effect.size
# 
# sim.d = getData(sim.methylBase[[1]])
# oobj = getData(sim.methylBase[[1]])
# Covs <- as.matrix(sim.d[, sim.methylBase[[1]]@coverage.index ])
# Meths <- as.matrix(sim.d[, sim.methylBase[[1]]@numCs.index ])
# a = data.frame(Type=c(rep("control",3), rep("treatment",3)),
#            Pair=c(paste0("sample", 1:3),paste0("sample", 1:3)))
# rownames(a) <- c(paste0(rep("C",3), 1:3), paste0(rep("T", 3), 1:3))
# BS1 <- BSseq(chr = oobj$chr, pos = oobj$start,
#                M = Meths, Cov = Covs, sampleNames = rownames(a),
#              pData=a)
# BS1.fit <- BSmooth(BS1, mc.cores = 30, verbose = TRUE)
# #It is open to personal preferences exactly which CpGs to remove, but for this analysis we
# #will only keep CpGs where at least 2 treatment samples and at least 2 control samples have at least 2x in
# #coverage.
# #samples. (2) Within each CpG island we considered regions to re-
# #ceive a DMR only if each of its CpG sites was covered in at least
# #half of all samples. Those regions are referred to as ‘covered is-
# BS.cov <- getCoverage(BS1.fit)
# keepLoci.ex <- which(rowSums(BS.cov[, BS1$Type == "treatment"] >= 2) >= 2 &
#                      rowSums(BS.cov[, BS1$Type == "control"] >= 2) >= 2)
# length(keepLoci.ex)
# BS1.fit <- BS1.fit[keepLoci.ex,]
# # compute t-statistic
# BS1.tstat <- BSmooth.tstat(BS1.fit,
#                     group1 = paste0(rep("T", 3), 1:3),
#                     group2 = paste0(rep("C",3), 1:3),
#                     estimate.var = "paired",
#                     local.correct = TRUE,
#                     verbose = TRUE)
# plot(BS1.tstat)
# #We rank DMRs by the column areaStat which is the sum of the t-statistics in each CpG. This is
# #kind of the area of the DMR, except that it is weighted by the number of CpGs and not by genomic
# #length. This is currently the best statistic we know, although it is far from perfect (we would like to do
# #something better).
# dmrs0 <- dmrFinder(BS1.tstat, 
#                    #cutoff = c(-4.5, 4.5))
#                    qcutoff = c(0.1, 0.90))
#                    #qcutoff = c(0.01, 0.99))
# 
# # Here, we filter out DMRs that do not have at least 3 CpGs in them and at least a mean difference
# # (across the DMR) in methylation between normal and cancers of at least 0.1
# dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
# 
# plotRegion(BS1.tstat, dmrs[1,], extend = 5000, addRegions = dmrs)
# 
# 
# 
