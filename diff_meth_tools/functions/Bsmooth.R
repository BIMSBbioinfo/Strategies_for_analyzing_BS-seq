
require("methylKit")
require("bsseq") #BSmooth


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

    dmrs.gr = makeGRangesFromDataFrame(dmrs0,
                                           keep.extra.columns=TRUE)
    return(dmrs.gr)
}

