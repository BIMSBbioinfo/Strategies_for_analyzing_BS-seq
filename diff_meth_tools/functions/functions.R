

# ---------------------------------------------------------------------------- #
#' Run limma for differential methylation calling
#' 
#' @param sim.methylBase a methylBase object from the methylKit library
#' @param transform a boolean indicating whether percent of methylation values
#'                  will be converted using a logistic transformation (Default: TRUE)
#' @return returns a methylDiff object from the methylKit library
limma.meth<-function(methylBase.obj,transform=TRUE){
  
  require(limma)
  require(qvalue)
  
  group<-methylBase.obj@treatment
  design<-model.matrix(~group)
  
  # do the test in limma
  p.meth=percMethylation(methylBase.obj) # get percent methylation values
  if(transform){
    p.meth2=log2((p.meth+1)/(100-p.meth+1))
    fit <- lmFit(p.meth2, design = design)
    fit2 <- eBayes(fit)
  }else{
    fit <- lmFit(p.meth, design = design)
    fit2 <- eBayes(fit)
  }
  # make the data for methylKit object
  df=cbind(getData(methylBase.obj)[,1:4],pvalue=fit2$p.value[,2],
           qvalue=qvalue(fit2$p.value[,2])$qvalues,
           meth.diff=rowMeans(p.meth[,methylBase.obj@treatment==1])-rowMeans(p.meth[,methylBase.obj@treatment==0])  )
  
  # create a new methylDiff object
  obj=new("methylDiff",df,sample.ids=methylBase.obj@sample.ids,
          assembly=methylBase.obj@assembly,
          context=methylBase.obj@context,treatment=methylBase.obj@treatment,
          destranded=methylBase.obj@destranded,resolution=methylBase.obj@resolution)
}



# ---------------------------------------------------------------------------- #
# DSS functions

run.DSS = function(methylBase.obj, difference=5, cores=1){
  
  print(cores)
  
  dss.qvalue = calculateDiffMethDSS(methylBase.obj,
                                    adjust="qvalue",
                                    mc.cores=cores)
  dss.fdr = calculateDiffMethDSS(methylBase.obj,
                                 adjust="fdr",
                                 mc.cores=cores)
  dss.slim = calculateDiffMethDSS(methylBase.obj,
                                  adjust="SLIM",
                                  mc.cores=cores)
  mylist.dss = list(dss.qvalue, dss.fdr, dss.slim)
  names( mylist.dss) = c("qvalue", "fdr", "slim")
  
  mylist.dss.diff.gr = mclapply(mylist.dss, 
                                function(x){
                                  as(getMethylDiff(x, difference=difference),"GRanges")},
                                mc.cores=3)
  names(mylist.dss.diff.gr) = c("qvalue", "fdr", "slim")
  
  return(list(mylist.dss, mylist.dss.diff.gr))
  
}

# ---------------------------------------------------------------------------- #
# methylKit functions

run.methylkit = function(sim.methylBase, cores=1,
                         difference=5){
  
  require("methylKit")
  adjusts = c("SLIM", "fdr", "qvalue")
  overdispersion = c("none",
                     "MN",
                     "shrinkMN"
  )
  tests = c("F", "Chisq")
  
  combined <- expand.grid(test=tests,
                          adjust=adjusts,
                          overd=overdispersion,
                          stringsAsFactors = FALSE)
  combined = cbind(combined,
                   name=with(combined, paste("methylKit",
                                             test, 
                                             adjust, 
                                             overd, sep=".")))
  combined$name = as.character(combined$name)
  
  methylKit.list=list()
  for(i in 1:nrow(combined)){
    co = combined[i,]
    methylkit.obj <- calculateDiffMeth(sim.methylBase, 
                                       overdispersion=co$overd,
                                       adjust = co$adjust,
                                       test=co$test,
                                       mc.cores=cores)
    methylkit.obj.diff = getMethylDiff(methylkit.obj, difference=difference)
    methylKit.list[[i]] <- as(methylkit.obj.diff, "GRanges")
  }
  names(methylKit.list) <- combined$name
  
  return(methylKit.list)
}

# ---------------------------------------------------------------------------- #
# BSmooth functions

smoothed.dmrFinder = function(smoothed){

    x=smoothed

    if( class( x@stats[,"tstat.corrected"] ) == "character"){
      x@stats = apply(x@stats,2,as.numeric)}
    if( any(is.na(x@stats)) ){
      #b = apply(x@stats, 1, function(y) any(is.na(y)))
      b = is.na(x@stats[,"tstat.corrected"])
      x@stats = x@stats[-which(b),]}

    bsmooth.cutoff4.5 <- function.error(try(my_dmrFinder(x,
                                             cutoff = c(-4.5, 4.5),
                                             qcutoff = NULL)) )
    bsmooth.cutoff4.5 = subset(bsmooth.cutoff4.5,
                                    n >= 3 & abs(meanDiff) >= 0.1)

    bsmooth.qcutoff1 <- function.error(try(my_dmrFinder(x,
                                             cutoff = NULL,
                                qcutoff = c(0.01, 0.99)) ))
    bsmooth.qcutoff1 = subset(bsmooth.qcutoff1,
                                    n >= 3 & abs(meanDiff) >= 0.1)

    bsmooth.qcutoff10 <- function.error(try(my_dmrFinder(x,
                                             cutoff = NULL,
                                 qcutoff = c(0.1, 0.9))) )
    bsmooth.qcutoff10 = subset(bsmooth.qcutoff10,
                                    n >= 3 & abs(meanDiff) >= 0.1)

    mylist = list(bsmooth.cutoff4.5, bsmooth.qcutoff1, bsmooth.qcutoff10)
    names(mylist) = c("bsmooth.cutoff4.5","bsmooth.qcutoff1","bsmooth.qcutoff10")
    return(mylist)
}

function.error = function(dmrs0){
  if(class( dmrs0)=="try-error"){
    print(dmrs0)
    return(NA)
  }
  if(is.null(dmrs0)){
    return(GRanges())
  }
  return(dmrs0)
}

smoothed.dmrFinder.paper = function(smoothed){
  
  x=smoothed
  
  if( class( x@stats[,"tstat.corrected"] ) == "character"){
    x@stats = apply(x@stats,2,as.numeric)}
  if( any(is.na(x@stats)) ){
    #b = apply(x@stats, 1, function(y) any(is.na(y)))
    b = is.na(x@stats[,"tstat.corrected"])
    x@stats = x@stats[-which(b),]}
  
  bsmooth.qcutoff <- function.error(try(my_dmrFinder(x,
                                                     cutoff = c(0.025, 0.975),
                                                     qcutoff = NULL)) )
  bsmooth.qcutoff.default_sub = subset(bsmooth.qcutoff,
                                       n >= 3 & abs(meanDiff) >= 0.1)
  
  list(bsmooth.same.default = bsmooth.qcutoff.default_sub)
}


run.BSmooth = function(methylBase.obj.list, cores=1,
                       difference=5){
  
  ctcf.bsmooth.smooth.same = lapply(methylBase.obj.list, function(x) {
    methylBase.obj = rm.small.chr.methylBase(x)
    runBSmooth.smooth(methylBase.obj,
                      estimate.var="same",
                      cores=60)
  })
  ctcf.bsmooth.smooth.group2 = lapply(methylBase.obj.list, function(x) {
    methylBase.obj = rm.small.chr.methylBase(x)
    runBSmooth.smooth(methylBase.obj,
                      estimate.var="group2",
                      cores=60)
  })
  
}


