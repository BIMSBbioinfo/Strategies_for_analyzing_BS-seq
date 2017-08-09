

#' Read methylation data from bigWig files
#'
#' @param methbw path/url to methylation bigwig file
#' @param covbw path/url to coverage bigwig file
#' @param chr a list of chromosome ids to extract data from
#' @param as.GRanges if TRUE, a GRanges object with meth and coverage columns 
#' are returned. If False a data.frame with chr,start,end,strand,coverage,numCs
#' and numTs columns are returned. This is essentially data.frame part
#' of a methylRaw object
#' @param skip.even some methylation bigwig files have coverage over both the 
#' C and the G of a CpG dinucleotide. Roadmap epigenomics files seems to be like that
#' If skip.even is TRUE only the C coordinate is retained, as keeping the G
#' coordinate and information will be redundant. 
#'
readbw<-function(methbw,covbw,chr="chr21",as.GRanges=FALSE,
                 skip.even=FALSE){
  
  require(rtracklayer)
  
  meth= BigWigFile(methbw)
  cov = BigWigFile(covbw)
  if(is.null(chr)){
    
    mbw = import(meth )
    cbw = import(cov )
  }else{
    
    chrlens=seqinfo(meth) # get all chrs
    len=seqlengths(chrlens[chr]) # length of chr
    mbw = import(meth, which=GRanges(seqnames=names(len),ranges=IRanges(1,len)))
    cbw = import(cov, which=GRanges(seqnames=names(len),ranges=IRanges(1,len)))
    
    
  }
  
  if(length(mbw) != length(cbw)){
    stop("lengths of coverage and methylation bigWig files do not match")
  }
  
  if(skip.even){
    mbw=mbw[seq(1,length(mbw),by=2),]
    cbw=cbw[seq(1,length(cbw),by=2),]
  }
  
  if(as.GRanges){
    mcols(mbw)$coverage=cbw$score
    colnames(mcols(mbw))[1]="meth"
    strand(mbw)="+"
    return(mbw)
  }else{
    numCs=round(mbw$score*cbw$score)
    numTs=cbw$score-numCs
    
    data.frame(chr=seqnames(mbw),start=start(mbw)
               ,end=end(mbw),strand="+",
               coverage=cbw$score,numCs=numCs,numTs=numTs)
  }
}

#' Join directly neighbouring segments of same class
#' 
#'
#' @param res object returned from a methSeg call
#'
#' @return res object
# @export
#'
# @examples
joinSegmentNeighbours <- function(res) {
  
  require(data.table)
  
  if (length(unique(seqnames(res))) > 1 ) {
    gr <- lapply(split(res,seqnames(res)),joinSegmentNeighbours)
    gr <- do.call(c, unlist(gr,use.names = FALSE) )
    return( gr )
  } 
  else if (length(res)<=1) { 
    return(res)
  }
  else{
    
    group.neighbours <- rle(res$seg.group)
    N = length(group.neighbours$lengths)
    # res_joined <- vector(mode="list", length=N)
    k <- numeric()
    k[1] <- 1
    l <- group.neighbours$lengths - 1
    for ( i in 2:N ) { k[i] = k[i-1] + l[ i -1] + 1 }
    #toJoin <- l==0
    
    # print(paste("k=",k,"l=",l))
    
    res_dt <- copy(as.data.table(res))
    
    for (i in which(l!=0)) {
      res_dt[k[i]:(k[i]+l[i]),`:=`(seqnames=unique(seqnames),
                                   start=min(start),
                                   end=max(end),
                                   strand = "*",
                                   width = sum(as.numeric(width)),
                                   ID = unique(ID),
                                   num.mark = sum(as.numeric(num.mark)),
                                   seg.mean = mean(seg.mean),
                                   startRow = min(startRow),
                                   endRow = max(endRow),
                                   seg.group = unique(seg.group))]
    }
    
    res_dt <- unique(res_dt)
    
    return(makeGRangesFromDataFrame(res_dt,keep.extra.columns = TRUE))
  }
}



