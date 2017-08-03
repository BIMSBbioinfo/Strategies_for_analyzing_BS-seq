
#' Visualise intersections between DMCs from given tools for calling differentiall methylated cytosines.
#'
#' @param list_toolsDMCs list of methylDiff diff or GRanges objects fom the methylKit library. The list should be named by the tools names.
#' @param vis boolean indicating whether intersection of DMCs between given tools will be visualizated using UpSetR::upset function.
#'
#' @return a matrix with first column that indicate positions of DMCs and other columns presence or absence of the DMCs from given tools.
plot.intersetion.tools = function(list_toolsDMCs, vis=TRUE, ...){
  
  require(dplyr)
  require(tidyr)
  
  DMCs2char = lapply(list_toolsDMCs, function(x){
    if(class(x)=="methylDiff"){
      paste0( x$chr,x$start,x$end )
    }else if(class(x)=="GRanges"){
      paste0( seqnames(x), start(x), end(x) )
    }
    
  } )
  
  charDMCs2df = lapply(1:length(DMCs2char), 
                       function(j){
                         # if there are DCs
                         if(length(DMCs2char[[j]])!=0){
                           data.frame(tool=names(DMCs2char)[j], pos=DMCs2char[[j]])
                         }
                       })
  d = do.call('rbind',charDMCs2df)
  binary_matrix <- d %>% mutate(value =1) %>% spread(tool, value, fill = 0 )
  
  if(vis){
    require(UpSetR)
    return( upset(binary_matrix, 
          order.by = "freq",
          mainbar.y.label = "DMCs intersections",
          sets.x.label = "Set size",
          ...)
    )
  }
  
  return( invisible(binary_matrix) )
}
