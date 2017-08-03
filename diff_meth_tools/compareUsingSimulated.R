#'---
#'title: "Strategies for analyzing bisulfite-seq data - on simulated data"
#'author: "Katarzyna Wreczycka"
#'output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#'    fig_width: 10
#'    fig_height: 10
#'---
#'
.libPaths(c("/home/kwreczy/Rlibs/3.5/"))
setwd("/home/kwreczy/projects/Strategies_for_analyzing_bisulfite-seq_data/")
library("methylKit")


############## 
############## FUNCTIONS
############## 


source("./functions/dataSim2.R")
source("./functions/assembly_tools.R")
source("./functions/compareUsingSim_func.R")



############## 
############ MAIN
############## 

#' # Calculate sensitivity and specificity of tools for differential methylation

# variables
effects = c(5, 10, 15, 20, 25)
models.res=list()
set.seed(111)

for(effect in effects){
  
  print(effect)
  
  
  # Generate simulated data
  sim.methylBase = dataSim2(replicates=6,
                                 sites=5000,
                                 treatment=c(1,1,1,0,0,0),
                                 percentage=50,
                                 effect=effect,
                                 add.info=TRUE)
  
  # Run models 
   models.res[[as.character(effect)]] = run.models(sim.methylBase, cores=20,
                                                   difference=5, qvalue=0.01)
  
}
models.res.save = models.res
saveRDS(models.res,"~/models.res.sim.21.07.2017.rds")
models.res = readRDS("~/models.res.sim.21.07.2017.rds")
#models.res = readRDS("~/models.res.sim.assembly.rds")



#' # Visualize sensitivity and specificity for different effect sizes

models.res.ma = do.call("rbind", models.res)
models.res.df = data.frame(models.res.ma)
models.res.df = cbind(models.res.df, 
                      tool=rownames(models.res.ma),
                      effect=as.factor(as.numeric(sapply(effects, 
                                               function(x) 
                                                 rep(x, nrow(models.res[[1]])  )))))

models.res.df = models.res.df[-which(models.res.df$tool=="methylKit.F.qvalue.none"),]
models.res.df$tool <- factor(models.res.df$tool)

levels(models.res.df$tool)[levels(models.res.df$tool)=="DSS.qvalue"] <- "DSS"
levels(models.res.df$tool)[levels(models.res.df$tool)=="limma.qvalue"] <- "limma"
levels(models.res.df$tool)[levels(models.res.df$tool)=="methylKit.Chisq.qvalue.MN"] <- "methylKit-Chisqtest-OC"
levels(models.res.df$tool)[levels(models.res.df$tool)=="methylKit.Chisq.qvalue.none"] <- "methylKit-Chisqtest"
levels(models.res.df$tool)[levels(models.res.df$tool)=="methylKit.F.qvalue.MN"] <- "methylKit-Ftest-OC"
#levels(models.res.df$tool)[levels(models.res.df$tool)=="assembly"] <- "DMCs in >= 2 tools"

#levels(models.res.df$tool)[levels(models.res.df$tool)=="methylKit.F.qvalue.none"] <- "methylKit-Ftest"
#models.res.df$tool=relevel(models.res.df$tool, "DMCs in >= 2 tools")
models.res.df$tool=relevel(models.res.df$tool, "methylKit-Ftest-OC")
models.res.df$tool=relevel(models.res.df$tool, "methylKit-Chisqtest-OC")
models.res.df$tool=relevel(models.res.df$tool, "methylKit-Chisqtest")
models.res.df$tool=relevel(models.res.df$tool, "limma")
models.res.df$tool=relevel(models.res.df$tool, "DSS")


p = models.res.df$TP / (models.res.df$TP + models.res.df$FP)
r = models.res.df$TP / (models.res.df$TP+models.res.df$FN)
models.res.df$f1_score =  2*((p*r)/(p+r))  
models.res.df$precision=as.numeric(p)
models.res.df$recall = as.numeric(r)
models.res.df$NPV = as.numeric(models.res.df$TN / (models.res.df$TN + models.res.df$FN))

# The palette with grey:
cbPalette <- c( "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


library(ggplot2)
p_sens<-ggplot(models.res.df,aes(effect, sens, fill=tool))+
  geom_bar(stat="identity",position='dodge',colour="black",
           width=0.65)+
  coord_cartesian(ylim=c(0.0,1.00))+
  scale_fill_manual(values=cbPalette)+
  labs(y="Sensitivity", x="Effect size (methylation difference)",fill='Tool')

p_spec<-ggplot(models.res.df,aes(effect, spec, fill=tool))+
  geom_bar(stat="identity",
           position="dodge",
           #position = "position_dodge(width = 0.9)",
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

p_fscore <- ggplot(models.res.df,aes(effect, f1_score, fill=tool))+
  geom_bar(stat="identity",position='dodge',colour="black",
           width=0.65)+
  coord_cartesian(ylim=c(0.0,1.00))+
  scale_fill_manual(values=cbPalette)+
  labs(y="F-score", x="Effect size (methylation difference)",fill='Tool')

# p_acc <- ggplot(models.res.df,aes(effect, acc, fill=tool))+
#   geom_bar(stat="identity",position='dodge',colour="black")+
#   coord_cartesian(ylim=c(0.0,1.00))+
#   labs(y="Accuracy", x="Effect",fill='Tool')
# p_precision <- ggplot(models.res.df,aes(effect, precision, fill=tool))+
#   geom_bar(stat="identity",position='dodge',colour="black")+
#   coord_cartesian(ylim=c(0.0,1.00))+
#   labs(y="Precision", x="Effect",fill='Tool')
# p_npv <- ggplot(models.res.df,aes(effect, NPV, fill=tool))+
#   geom_bar(stat="identity",position='dodge',colour="black")+
#   coord_cartesian(ylim=c(0.0,1.00))+
#   labs(y="Negative predictive value", x="Effect",fill='Tool')


grid.arrange(p_sens, p_spec, p_fscore, ncol = 2, nrow=2)


library(gridExtra)
library(grid)
p_sens1 <- arrangeGrob(p_sens, top = textGrob("a", x=unit(0, "npc"),y=unit(1, "npc"),
                                              just=c("left","top"), gp=gpar(col="black", fontsize=14)))
p_spec1 <- arrangeGrob(p_spec, top = textGrob("b", x=unit(0, "npc"),y=unit(1, "npc"),just=c("left","top"), gp=gpar(col="black", fontsize=14)))
#p_acc1 <- arrangeGrob(p_acc, top = textGrob("c", x=unit(0, "npc"),y=unit(1, "npc"),just=c("left","top"), gp=gpar(col="black", fontsize=14)))
p_fscore1 <- arrangeGrob(p_fscore, top = textGrob("c", x=unit(0, "npc"),y=unit(1, "npc"),just=c("left","top"), gp=gpar(col="black", fontsize=14)))

# https://www.elsevier.com/journals/electronic-journal-of-biotechnology/0717-3458/guide-for-authors
# Arial, Courier, Times New Roman, Symbol

pdf("/data/akalin/kwreczy/projects/Strategies_for_analyzing_bisulfite-seq_data/simulated_rates.pdf",
    width = 10, height = 7)
library(gridExtra)
grid.arrange(p_sens1, p_spec1, p_fscore1, ncol = 2, nrow=2)
dev.off()

# ensemble, DMC in at least 2 of tools, see if ensemble method is better than 
# f-score general measure of acc



#rmarkdown::render("/home/kwreczy/projects/Strategies_for_analyzing_bisulfite-seq_data/compareUsingSimulated.R")



