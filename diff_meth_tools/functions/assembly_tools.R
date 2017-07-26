#' vec a vector of characters
find.combinations = function(vec){
  
  all.comb = expand.grid(vec,
                         vec,
                         KEEP.OUT.ATTRS = FALSE,
                         stringsAsFactors = FALSE)
  comb = all.comb[!duplicated(apply(all.comb,1,function(x) paste(sort(x),collapse=''))),]
  # Remove rows where cell line are the same
  wh.re = which(comb[,1]==comb[,2])
  comb = comb[-wh.re,] 
  
  return(comb)
  
}


#' Find differentially methylated CpG sites in at least 2 tools
#' 
#' ... a list of GRanges or GRangesList and GRanges indicates genomic sites
assembly_tools = function(grl, cores=1){
  
  # Find all DMC from all methods
  # to them check which of them appear in output
  # of at least 2 tools
  un.grl = unlist(grl)
  uniq.un.grl = unique(unlist(grl))
  
  if( any(width(uniq.un.grl) > 2) ){
    stop("GRanges need to have width not bigger than 2")
  }
  strand(uniq.un.grl) = "*"
  comb = find.combinations( names(grl) )
  
  # Convert Granges to characters like chr.start.end
  # It's easier to work on such characters than on GRanges, i.e.
  # overlap granges and if DMV are next to each other they will be 
  # merged together.
  all.dms.char=paste(seqnames(uniq.un.grl),start(uniq.un.grl), sep=".")
  grl.char = lapply(grl, function(gr){
    paste0(seqnames(gr),start(gr), end(gr))    
  })
  names(grl.char) <- names(grl)
  
  # Take output of 2 tools and find common DMCs
  dmcs.comb.char = mclapply(1:nrow(comb), function(i){
    v1 = as.character(comb[i,][1])
    v2 = as.character(comb[i,][2])
    dmcs.2.tools = unique(intersect(grl[[v1]], grl[[v2]]))
    dmcs.2.tools.char = paste(seqnames(dmcs.2.tools),start(dmcs.2.tools), sep=".")
    return(dmcs.2.tools.char)
  }, mc.cores=cores)
  
  # Check if DMC is in output of at least 2 tools (TRUE) or not (FALSE)
  dmcs.in.atleast2tools.binary = unlist(mclapply(1:length(all.dms.char), function(i){
    dmc = all.dms.char[i]
    for(j in 1:length(dmcs.comb.char)){
      if(dmc %in% dmcs.comb.char[[j]]){
        return(TRUE)
        break
      }}
    return(FALSE)
  }, mc.cores=cores))
  dmcs.in.atleast2tools = all.dms.char[which(dmcs.in.atleast2tools.binary)]
  
  # Convert characters like chr1.1.3 to GRanges object
  cha2df = lapply(1:length(dmcs.in.atleast2tools), function(i){
    unlist(strsplit(dmcs.in.atleast2tools[i], "[.]"))
  })
  df = as.data.frame(do.call("rbind", cha2df),
                     stringsAsFactors = FALSE)
  df[,3] = df[,2] # end=start
  colnames(df) <- c("chr", "start", "end")
  df$start = as.numeric(df$start)
  df$end = as.numeric(df$start) + 1 # end=start+1
  
  makeGRangesFromDataFrame(df)
  
}  

# An example:
# f <- function(..., cores=1) {
#   list(...)
# }
# f(a = 1, b = 2)
# toy dataset
# gr1 <-
#   GRanges(seqnames = "chr1", ranges = IRanges(start=c(7,15,20), width = 1),
#           strand = "*")
# gr2 <-
#   GRanges(seqnames = "chr1",
#           ranges = IRanges(start=c(7,13,20), width = 1),
#           strand = "*")
# gr3 <-
#   GRanges(seqnames = "chr1",
#           ranges = IRanges(start=c(15,21,30), width = 1),
#           strand = "*")
# grl <- GRangesList("gr1" = gr1, "gr2" = gr2, "gr3" = gr3)
# grl
# assembly_tools(grl)