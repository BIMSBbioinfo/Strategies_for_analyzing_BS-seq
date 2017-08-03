# ---------------------------------------------------------------------------- #
#' Run limma for differential methylation calling. 
#' 
#' It uses limma R package to fit linear model for each cytosine within given samples 
#' (limma::lmFit) and computes moderated t-statistics, moderated F-statistic, and 
#' log-odds of differential methylation by empirical Bayes moderation of the 
#' standard errors towards a common value (limma::eBayes).
#' 
#' @param sim.methylBase a methylBase object from the methylKit library
#' @param transform a boolean indicating whether percent of methylation values
#'                  will be converted using a logistic transformation 
#'                  (Default: TRUE)
#' @param ... argument(s) that are passed to the limma::eBayes function, 
#'            such as trend, robust, etc.
#' 
#' @author Altuna Akalin, Katarzyna Wreczycka
#' 
#' @return a methylDiff object from the methylKit library
limma.meth<-function(methylBase.obj,transform=TRUE, ...){
  
  require(limma)
  require(qvalue)
  
  group<-methylBase.obj@treatment
  design<-model.matrix(~group)
  
  # do the test in limma
  p.meth=percMethylation(methylBase.obj) # get percent methylation values
  if(transform){
    p.meth2=log2((p.meth+1)/(100-p.meth+1))
    fit <- lmFit(p.meth2, design = design)
    fit2 <- eBayes(fit, ...)
  }else{
    fit <- lmFit(p.meth, design = design)
    fit2 <- eBayes(fit, ...)
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

