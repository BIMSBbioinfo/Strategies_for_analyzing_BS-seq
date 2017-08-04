

#' Returns percentage of peaks from peaks.gr overlapping with at least
#' one region from diff.gr
overlp.diff.with.supp.peaks = function(diff.gr, peaks.gr){
  if(class(diff.gr)!="GRanges" | class(peaks.gr)!="GRanges"){return(NA)}
  overlp = findOverlaps(
    peaks.gr,
    diff.gr )
  length(unique(queryHits(overlp))) / length(peaks.gr) * 100
}
create.empty.matrix = function(axisname){
  matrix(0, nrow = length(axisname), 
         ncol = length(axisname),
         dimnames = list(axisname,axisname))
}

#' Returns numeric matrix with percent of differentially methylated region
#' calculated by a tool with at least one DMC
matrix.perc.with.dmc = function(models, comb, 
                                names.all.celllines, 
                                supp.tbl2.gr,
                                rrbs.methylBase.obj.list){
  my.matrix.both1 = create.empty.matrix(names.all.celllines)
  my.matrix.one1one0 = create.empty.matrix(names.all.celllines)
  for(i in 1:nrow(comb)){
    c1=comb[i,]$cellline1
    c2=comb[i,]$cellline2
    #print(paste(c1, c2))
    
    # find peaks 
    peaksc1c2 = aggregate.peaks(c1, c2, supp.tbl2.gr)
    # take peaks that have RRBS coverage
    mi = as(rrbs.methylBase.obj.list[[i]],"GRanges")
    peaksc1c2$both.one = subsetByOverlaps(peaksc1c2$both.one, mi)
    peaksc1c2$onezero.oneone = subsetByOverlaps(peaksc1c2$onezero.oneone, mi)
    
    a=overlp.diff.with.supp.peaks(models[[i]], peaksc1c2$both.one)
    b=overlp.diff.with.supp.peaks(models[[i]], peaksc1c2$onezero.oneone)
    my.matrix.both1[which(rownames(my.matrix.both1)==c1),
                    which(colnames(my.matrix.both1)==c2)] = a
    my.matrix.one1one0[which(rownames(my.matrix.one1one0)==c1),
                       which(colnames(my.matrix.one1one0)==c2)] = b
    
  }
  list(my.matrix.both1, my.matrix.one1one0 )
}

#' check accuracy of models compared to tools for DM detection
#' 
#' @param all.peaks a GRanges object indicating all peaks the the suppl. table
#'                  no matter if they are 'truly' sensitive or non-sensitive to methylation
#' @param peak.change a GRanges object indicating which genomic regions 
#'                    are 'truly' sensitive to methylation
#' @param pred.DM a GRanges object indicating which genomic regions 
#'                    are estimated to be sensitive to methylation
#' @return returns a vector of accuracy metrics, TP, FP, Sensivity, etc
check.accuracy.CTCF<-function(all.peaks, peaks.change,
                              pred.DM){
  
  all = paste0( seqnames(all.peaks),
                start(all.peaks),
                end(all.peaks))
  
  true.dm=paste0( seqnames(peaks.change),
                  start(peaks.change),
                  end(peaks.change))
  true.ndm = all[! all %in% true.dm]
  
  pred.dm=paste0( seqnames(pred.DM),
                  start(pred.DM),
                  end(pred.DM))
  pred.ndm =all[! all %in% pred.dm]
  
  TP=sum(true.dm %in% pred.dm)
  
  FN=sum(pred.ndm %in% true.dm)
  
  FP=sum(pred.dm %in% true.ndm)
  
  TN=sum(pred.ndm %in% true.ndm)
  
  return(c(TP=TP,FN=FN,FP=FP,TN=TN,
           acc=(TP+TN)/length(all),
           spec=TN/(TN+FP) ,
           sens=TP/(TP+TN) ) )
}


#' Aggregate peaks from the suppl. table S2
#' to peaks that change (1:0/0:1) and peaks that doesn't change (1:1)
#'
#' @c1, @c2 characters indicating cell line
#' @gr a GRanges object
#' @return returns a list of 
aggregate.peaks = function(c1, c2, gr){
  
  gr.c1c2 = gr[ ,which( colnames(mcols(gr)) %in% c(c1, c2) ) ]
  which.both.one = which(mcols(gr.c1c2[,1])[,1]==1 & mcols(gr.c1c2[,2])[,1]==1)
  both.one = gr.c1c2[which.both.one,]
  which.both.zero = which(mcols(gr.c1c2[,1])[,1]==0 & mcols(gr.c1c2[,2])[,1]==0)
  both.zero = gr.c1c2[which.both.zero,]
  onezero.oneone = gr.c1c2[-c(which.both.one, which.both.zero)]
  
  list(both.one=both.one, 
       both.zero=both.zero,
       onezero.oneone=onezero.oneone)
}

#' @param methylBase.obj.list a list of methylBase objects for each combination of cell lines
#' @param comb a data frame that have combinations of names of cell lines
#' @param supp.tbl2.gr a granges object witch columns containing 1 or 0 for a cell line 
#' @return returns a list of vectors of accuracy metrics, TP, FP, Sensivity, etc
rates.cell.lines = function(methylBase.obj.list, comb,
                            supp.tbl2.gr,
                            pred.dm.per.cell.line,
                            cores=20){
  
  my.list = mclapply(1:nrow(comb), function(i){
    
    co1 = comb[i,]
    c1=co1[1];c2=co1[2]
    
    # extract peaks (1:0/0:1, 1:1 and 0:0)
    peaksc1c2 = aggregate.peaks(c1, c2, supp.tbl2.gr)
    # take peaks that are covered by RRBS from the ucsc/encode
    mi = as(methylBase.obj.list[[i]],"GRanges")
    peaksc1c2$both.one = subsetByOverlaps(mi, peaksc1c2$both.one)
    peaksc1c2$both.zero = subsetByOverlaps(mi, peaksc1c2$both.zero)
    peaksc1c2$onezero.oneone = subsetByOverlaps(mi, peaksc1c2$onezero.oneone)
    all.peaks=unlist(GRangesList(peaksc1c2))
    
    ca = check.accuracy.CTCF(all.peaks=all.peaks,
                             peaks.change=peaksc1c2$onezero.oneone,
                             pred.DM=pred.dm.per.cell.line[[i]])
    ca = c(ca,
           pred.DM=length(pred.dm.per.cell.line[[i]]),
           true.DM=length(peaksc1c2$onezero.oneone),
           all=length(all.peaks))
    my.list[[i]] = ca
  }, mc.cores=cores)
  do.call("rbind", my.list)
}





calc.metrics.atleast1DM.GR = function(methylBase.obj, 
                                      co,
                                      supp.tbl2.gr,
                                      pred.dm.per.cell.line,
                                      sample_no_change_peaks=FALSE){
  c1=co[1];c2=co[2]
  
  # find peaks
  peaksc1c2 = aggregate.peaks(c1, c2, supp.tbl2.gr)
  # take peaks that have RRBS coverage
  mi = as(methylBase.obj,"GRanges")
  peaksc1c2$both.one = subsetByOverlaps(peaksc1c2$both.one, mi, ignore.strand=TRUE)
  peaksc1c2$both.zero = subsetByOverlaps(peaksc1c2$both.zero, mi, ignore.strand=TRUE)
  peaksc1c2$onezero.oneone = subsetByOverlaps(peaksc1c2$onezero.oneone, mi, ignore.strand=TRUE)
  
  all.peaks=unlist(GRangesList(peaksc1c2))
  peaks.dont.change = c(peaksc1c2$both.one,peaksc1c2$both.zero)
  peaks.change = peaksc1c2$onezero.oneone
  
  if(sample_no_change_peaks){
    # subset peaks that dont change to have the same number of elements
    # like set of peaks that change
    # because if we have way larger set peaks that dont change than change
    # then we will detect, e.g. more FN than TP
    indx.sampled = sample(1:length(peaks.dont.change), length(peaks.change))
    peaks.dont.change.o = peaks.dont.change
    peaks.dont.change = peaks.dont.change[indx.sampled]
  }
  
  pred.dm = pred.dm.per.cell.line
  # propably I should fixed somewhere earlier in this analysis
  # but for now I do it here
  # if pred.dm is not a methylDiff object, is a BSmooth object, add a column
  # that would imitate the object
  if( !("meth.diff" %in% colnames(mcols(pred.dm))) )
  {
    if(!(is.null(pred.dm)))
      if(length(pred.dm)!=0)
          pred.dm$meth.diff = mcols(pred.dm)$meanDiff*100
  }
  
  fi.peakschange.overlp.pred.dm = findOverlaps(peaks.change,
                                               pred.dm,ignore.strand=TRUE)
  indx.peakschange.overlp.pred.dm = unique(queryHits(fi.peakschange.overlp.pred.dm))
  fi.peaksnochange.overlp.pred.dm = findOverlaps(peaks.dont.change,
                                                 pred.dm,ignore.strand=TRUE)
  indx.peaksnochange.overlp.pred.dm = unique(queryHits(fi.peaksnochange.overlp.pred.dm))
  
  if(length(indx.peakschange.overlp.pred.dm)!=0){
    # peak change + DM
    TP.gr = peaks.change[  indx.peakschange.overlp.pred.dm  ]
    # peak change + no-DM
    FN.gr = peaks.change[-indx.peakschange.overlp.pred.dm  ]
  }else{
    # peak change + DM
    TP.gr = integer(0)
    # peak change + no-DM
    FN.gr = peaks.change
  }
  if(length(indx.peaksnochange.overlp.pred.dm)!=0){
    # no peak change  + DM
    FP.gr = peaks.dont.change[ indx.peaksnochange.overlp.pred.dm  ]
    # no peak change + no-DM
    TN.gr = peaks.dont.change[- indx.peaksnochange.overlp.pred.dm  ]
  }else{
    # no peak change  + DM
    FP.gr = integer(0)
    # no peak change + no-DM
    TN.gr = peaks.dont.change
  }
  return(
    list(TP.peaks = TP.gr,
         FP.peaks = FP.gr,
         FN.peaks = FN.gr,
         TN.peaks = TN.gr,
         TP.meth = pred.dm[unique(subjectHits(fi.peakschange.overlp.pred.dm))],
         FP.meth = pred.dm[unique(subjectHits(fi.peaksnochange.overlp.pred.dm))]
    )
  )
  
}



calc.metrics.atleast1DM = function(methylBase.obj, co,
                                   supp.tbl2.gr,
                                   pred.dm.per.cell.line,
                                   sample_no_change_peaks=FALSE){
  c1=co[1];c2=co[2]
  
  # find peaks
  peaksc1c2 = aggregate.peaks(c1, c2, supp.tbl2.gr)
  # take peaks that have RRBS coverage
  mi = as(methylBase.obj,"GRanges")
  
  peaksc1c2$both.one = subsetByOverlaps(peaksc1c2$both.one, mi, ignore.strand=TRUE)
  peaksc1c2$both.zero = subsetByOverlaps(peaksc1c2$both.zero, mi, ignore.strand=TRUE)
  peaksc1c2$onezero.oneone = subsetByOverlaps(peaksc1c2$onezero.oneone, mi, ignore.strand=TRUE)
  
  all.peaks=unlist(GRangesList(peaksc1c2))
  peaks.dont.change = c(peaksc1c2$both.one,peaksc1c2$both.zero)
  peaks.change = peaksc1c2$onezero.oneone
  
  if(sample_no_change_peaks){
    # subset peaks that dont change to have the same number of elements
    # like set of peaks that change
    # because if we have way larger set peaks that dont change than change
    # then we will detect, e.g. more FN than TP
    indx.sampled = sample(1:length(peaks.dont.change), length(peaks.change))
    peaks.dont.change.o = peaks.dont.change
    peaks.dont.change = peaks.dont.change[indx.sampled]
  }
  
  pred.dm = pred.dm.per.cell.line
  
  indx.peakschange.overlp.pred.dm = unique(queryHits(findOverlaps(peaks.change,
                                                                  pred.dm,ignore.strand=TRUE)))
  indx.peaksnochange.overlp.pred.dm = unique(queryHits(findOverlaps(peaks.dont.change,
                                                                    pred.dm,ignore.strand=TRUE)))
  
  if(length(indx.peakschange.overlp.pred.dm)!=0){
    # peak change + DM
    TP.gr = peaks.change[  indx.peakschange.overlp.pred.dm  ]
    # peak change + no-DM
    FN.gr = peaks.change[-indx.peakschange.overlp.pred.dm  ]
  }else{
    # peak change + DM
    TP.gr = integer(0)
    # peak change + no-DM
    FN.gr = peaks.change
  }
  if(length(indx.peaksnochange.overlp.pred.dm)!=0){
    # no peak change  + DM
    FP.gr = peaks.dont.change[ indx.peaksnochange.overlp.pred.dm  ]
    # no peak change + no-DM
    TN.gr = peaks.dont.change[- indx.peaksnochange.overlp.pred.dm  ]
  }else{
    # no peak change  + DM
    FP.gr = integer(0)
    # no peak change + no-DM
    TN.gr = peaks.dont.change
  }
  TP = length(TP.gr); FN=length(FN.gr); FP=length(FP.gr); TN = length(TN.gr)
  
  # https://en.m.wikipedia.org/wiki/F1_score
  # https://en.wikipedia.org/wiki/Precision_and_recall
  # precision is the number of correct positive results divided by
  # the number of all positive results
  # = how many selected items are relevant
  p = TP / (TP + FP)
  # recall is the number of correct positive results divided by
  # the number of positive results that should have been returned.
  # = how many relevant items are selected
  r = TP / (TP+FN)
  f_score = 2*((p*r)/(p+r))
  
  my.vec = c(TP=TP,FN=FN,FP=FP,TN=TN,
             acc=(TP+TN)/length(all.peaks),
             spec=TN/(TN+FP) ,
             sens=TP/(TP+FN), ### here was a buuug  TP/(TP+TN)
             f_score=f_score,
             all.peaks=length(all.peaks),
             pred.dm=length(pred.dm),
             precision=as.numeric(TP / (TP + FP)),
             NPV=as.numeric(TN / (TN + FN)),
             recall=as.numeric(r))
  my.vec
}
