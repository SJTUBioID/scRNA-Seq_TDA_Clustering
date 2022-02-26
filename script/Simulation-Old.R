library(splatter)
library(scater)
library(ggplot2)
library(SingleCellExperiment)

#browseVignettes("splatter")
setwd("/Users/mac/Desktop/R_Splatter")
real = read.table("splatter-paper/data/Klein.csv",sep=",", header = T,row.names = 1,check.names = F)
real <- real[rowSums(real) > 0, ]
real <- as.matrix(real)
params <- splatEstimate(real)





## Download sampling
Simulation1 = function(seed)
{
  colors1=c("DodgerBlue4","MediumVioletRed","DarkTurquoise","MediumSeaGreen","DarkMagenta","Snow3","Maroon")
  sim <- splatSimulate(params, 
                       batchCells      = 20000,
                       group.prob      = c(0.5,0.5),
                       de.prob         = 0.01,
                       de.facLoc       = 0.1,
                       de.facScale     = 0.3,
                       method = "groups")
  sim <- calculateQCMetrics(sim)
  sim@metadata$Params
  sim <- normalize(sim)
  plotPCA(sim, colour_by = "Group")
  #plotTSNE(sim, colour_by = "Group")
  #plotMDS(sim,colour_by="Group")
  Count = data.frame(counts(sim))

  write.table(Count,     "Splatter1_20000.txt",quote=F,sep="\t")
  write.table(sim$Group, "Splatter1_20000.label",quote=F,row.names = F,sep="\t")
  write.table(Count[,1:18000],     "Splatter1_18000.txt",quote=F,sep="\t")
  write.table(sim$Group[1:18000],  "Splatter1_18000.label",quote=F,row.names = F,sep="\t")
  write.table(Count[,1:16000],     "Splatter1_16000.txt",quote=F,sep="\t")
  write.table(sim$Group[1:16000],  "Splatter1_16000.label",quote=F,row.names = F,sep="\t")
  write.table(Count[,1:14000],     "Splatter1_14000.txt",quote=F,sep="\t")
  write.table(sim$Group[1:14000],  "Splatter1_14000.label",quote=F,row.names = F,sep="\t")
  write.table(Count[,1:12000],     "Splatter1_12000.txt",quote=F,sep="\t")
  write.table(sim$Group[1:12000],  "Splatter1_12000.label",quote=F,row.names = F,sep="\t")
  write.table(Count[,1:10000],     "Splatter1_10000.txt",quote=F,sep="\t")
  write.table(sim$Group[1:10000],  "Splatter1_10000.label",quote=F,row.names = F,sep="\t")
  write.table(Count[,1:8000],      "Splatter1_8000.txt",quote=F,sep="\t")
  write.table(sim$Group[1:8000],   "Splatter1_8000.label",quote=F,row.names = F,sep="\t")
  write.table(Count[,1:6000],      "Splatter1_8000.txt",quote=F,sep="\t")
  write.table(sim$Group[1:6000],   "Splatter1_8000.label",quote=F,row.names = F,sep="\t")
  write.table(Count[,1:4000],      "Splatter1_4000.txt",quote=F,sep="\t")
  write.table(sim$Group[1:4000],   "Splatter1_4000.label",quote=F,row.names = F,sep="\t")
  write.table(Count[,1:2000],      "Splatter1_2000.txt",quote=F,sep="\t")
  write.table(sim$Group[1:2000],   "Splatter1_2000.label",quote=F,row.names = F,sep="\t")
  
} 

Simulation1(seed=100)


## Download sampling
Simulation1 = function(seed)
{
  colors1=c("DodgerBlue4","MediumVioletRed","DarkTurquoise","MediumSeaGreen","DarkMagenta","Snow3","Maroon")
  sim <- splatSimulate(params, 
                       batchCells      = 10000,
                       group.prob      = c(0.5,0.5),
                       de.prob         = 0.01,
                       de.facLoc       = 0.1,
                       de.facScale     = 0.3,
                       method = "groups")
  sim <- logNormCounts(sim)
  sim <- runPCA(sim)
  plotPCA(sim, colour_by = "Group")
  sim <- runTSNE(sim)
  plotTSNE(sim, colour_by = "Group")
  #plotMDS(sim,colour_by="Group")

  Count <- as.data.frame(as.matrix(counts(sim)))
  write.table(Count,     "Splatter1_10000.txt",quote=F,sep="\t")
  write.table(sim$Group, "Splatter1_10000.label",quote=F,row.names = F,sep="\t")
  write.table(Count[,1:9000],     "Splatter1_9000.txt",quote=F,sep="\t")
  write.table(sim$Group[1:9000],  "Splatter1_9000.label",quote=F,row.names = F,sep="\t")
  write.table(Count[,1:8000],     "Splatter1_8000.txt",quote=F,sep="\t")
  write.table(sim$Group[1:8000],  "Splatter1_8000.label",quote=F,row.names = F,sep="\t")
  write.table(Count[,1:7000],     "Splatter1_7000.txt",quote=F,sep="\t")
  write.table(sim$Group[1:7000],  "Splatter1_7000.label",quote=F,row.names = F,sep="\t")
  write.table(Count[,1:6000],     "Splatter1_6000.txt",quote=F,sep="\t")
  write.table(sim$Group[1:6000],  "Splatter1_6000.label",quote=F,row.names = F,sep="\t")
  write.table(Count[,1:5000],     "Splatter1_5000.txt",quote=F,sep="\t")
  write.table(sim$Group[1:5000],  "Splatter1_5000.label",quote=F,row.names = F,sep="\t")
  write.table(Count[,1:4000],      "Splatter1_4000.txt",quote=F,sep="\t")
  write.table(sim$Group[1:4000],   "Splatter1_4000.label",quote=F,row.names = F,sep="\t")
  write.table(Count[,1:3000],      "Splatter1_3000.txt",quote=F,sep="\t")
  write.table(sim$Group[1:3000],   "Splatter1_3000.label",quote=F,row.names = F,sep="\t")
  write.table(Count[,1:2000],      "Splatter1_2000.txt",quote=F,sep="\t")
  write.table(sim$Group[1:2000],   "Splatter1_2000.label",quote=F,row.names = F,sep="\t")
  
} 
Simulation1(seed=100)