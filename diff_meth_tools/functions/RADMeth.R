

.libPaths(c("/home/kwreczy/Rlibs/3.4/"))

# RADMeth needs two input files
# 1. The proportion table, such as:
#control a control b control c case a case b case c
#chr1:108:109 9 6 10 8 1 1 2 2 2 1 14 1
#chr1:114:115 17 7 10 0 14 3 5 1 9 1 7 1
#chr1:160:161 12 8 10 5 17 4 15 14 13 6 4 4
#chr1:309:310 1 1 1 0 17 12 12 8 2 1 19 8
#chr1:499:500 8 4 6 5 15 6 14 10 14 11 15 1
#chr1:510:511 0 0 0 0 14 8 4 0 5 3 5 1
# 2. The design matrix for this dataset describes the structure of the experiment
# base case
# control a 1 0
# control b 1 0
# control c 1 0
# case a 1 1
# case b 1 1
# case c 1 1

create_proportion_table = function(simobj, outputpath){
  
  simd = getData(simobj)
  col0 = paste0(simd$chr,":", simd$start,":", simd$end)
  noT = simd[,-c(1,2,3,4,simobj @ numTs.index)]
  noT = cbind(col0, noT)
  write.table(noT, file = outputpath,
              append = FALSE, quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE,
              col.names = FALSE, qmethod = c("escape", "double"),
              fileEncoding = "")
  
  add1stline2input=paste0("sed -i ", "'1i ", 
                          paste0(simobj@sample.ids, collapse=" "), "' ",
                          outputpath)
  system(add1stline2input)
}

create_design_matrix = function(simobj, outputpath){
  
  df = data.frame(base=rep(1,length(simobj@ treatment)),
                  treatment=simobj@ treatment, 
                  row.names = simobj@sample.ids)
  write.table(df, file = outputpath,
              append = FALSE, quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = TRUE,
              col.names = TRUE, qmethod = c("escape", "double"),
              fileEncoding = "")
}


runRADMeth = function(methylBase.obj){
  # Input is a methylBase object from methylkit
  # It return granges with regions that are differentially methylated
  
  
  prop.table.path = "~/projects/Strategies_for_analyzing_bisulfite-seq_data/proportion_table.txt"
  design.matrix.path = "~/projects/Strategies_for_analyzing_bisulfite-seq_data/design_matrix.txt"
  regr.path = "~/projects/Strategies_for_analyzing_bisulfite-seq_data/regression.bed"
  regr.adjusted.path = "~/projects/Strategies_for_analyzing_bisulfite-seq_data/regression_adjusted.bed"
  
  create_proportion_table(methylBase.obj,
                          prop.table.path)
  create_design_matrix(methylBase.obj,
                       design.matrix.path)
  
  
  RADMeth.bin = "/home/kwreczy/programs/radmeth/bin/"
  
  # The -factor parameter specifies the factor with respect to
  # which we want to test for the differential methylation
  # The test factor is case, meaning that we are testing for
  # differential methylation between cases and controls.
  # /home/kwreczy/programs/radmeth/bin/wand -factor treatment ./design_matrix.txt ./proportion_table.txt > ./tmp.bed
  system(paste0(RADMeth.bin,"wand ",
                "-factor treatment ",
                design.matrix.path, " ",
                prop.table.path, " ",
                " > ",
                regr.path))
  # header: 
  # chrom start end c:log-odds-ratio:mean-meth-diff pval
  
  # We do not use these p-values directly, but instead we adjust the p-value of each CpG site based on the p-values
  # of the neighboring CpGs.
  # Here, the only required parameter, besides the input file, is -bins whose value is set to 1:100:1. This means
  # that for each n = 1, 2, . . . 99, RADMeth will compute correlation between p-values of cpgs located at distance
  # n from each other.
  # /home/kwreczy/programs/radmeth/bin/adjust -bins 1:100:1 ./tmp.bed > ./tmp.adjusted.bed
  system(paste0(RADMeth.bin,"adjust ",
                "-bins 1:100:1 ",
                regr.path, " ",
                " > ",
                regr.adjusted.path))
  
  # header
  # chrom start end c:log-odds-ratio:mean-meth-diff:pval:combined-pval fdr-pval
  # where all of the parameters, except combined-pval and fdr-pval, are as before; combined-pval is the p-value
  # given by the Z test which combines pvals from proximal CpG sites and fdr-pval is the FDR corrected
  # combined-pval.
  
  # Invidual differentially methylated sites 
  # To get all CpGs with FDR-corrected p-value below 0.01
  # awk '$5 < 0.01 "{ print $0; $}"' tmp.adjusted.bed > dm_tmp.bed
  tbl = read.table(regr.adjusted.path,
                   col.names = c("chrom", "start", "end", 
                                 "c:log-odds-ratio:mean-meth-diff:pval:combined-pval",
                                 "fdr-pval"))
  tbl.corr = tbl[which(tbl$fdr.pval < 0.01),]

  # Differentially methylated regions
  # RADMeth can also join individually differentially methylated CpGs into
  # differentially methylated regions. it joins neighboring differentially methylated sites with p-value
  # below 0.01 (set by the -p paramter).
  # ∼/radmeth/bin/dmrs -p 0.01 cpgs.adjusted.bed > dmrs.bed
  
  tbl.corr$strand = "*"
  tbl.corr.gr = makeGRangesFromDataFrame(tbl.corr,keep.extra.columns=TRUE)
  
  # clean up
  # file.remove(prop.table.path)
  # file.remove(design.matrix.path)
  # file.remove(regr.path)
  # file.remove(regr.adjusted.path)

  return(tbl.corr.gr)

}


#### Test

# library("methylKit")
# sim.methylBase = get.simdata(25, sites=10000)
# simobj=sim.methylBase$methylbase
# diff=sim.methylBase$diff
# no.diff=sim.methylBase$nodiff
# methylBase.obj=simobj
# 
# start.time <- Sys.time()
# runRADMeth(methylBase.obj)
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
# #Time difference of 2.246054 mins
# 






