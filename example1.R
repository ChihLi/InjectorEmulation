### Loading libraries
library(doParallel)
library(snowfall)
library(foreach)
library(bigmemory)
library(biganalytics)
library(rARPACK)
library(nloptr)
library(huge)
library(SLHD)
library(colorRamps)
library(shape)
library(doMPI)
library(InjectorEmulationHelper)

### Settings
Rawdata.folder.path  <- "data1/rawdata"
Testdata.folder.path <- "data1/testdata"
design.file.path     <- "data1/DoE.csv"
new.design.file.path <- "data1/newX.csv"
output.folder.path   <- "output1"
time.step            <- 1:3
Response.ColumnIndex <- 4:9          # indicate which columns present responses in the raw data     
Response.ColumnNames <- c("u", "v", "w", "P", "T", "rho")       
POD.energy           <- 0.99         # the proportion of total flow energy for POD truncation
confidence.level     <- 0.80         # indicate the level of confidence interval
dompi                <- FALSE        # do mpi parallel or not  

### set up for parallel
if(dompi){
  cl <- startMPIcluster(count = 2, verbose = TRUE, maxcores = detectCores())
  registerDoMPI(cl)
}else{
  sfInit(parallel = TRUE, cpus = detectCores())
  cl <- sfGetCluster()
  registerDoParallel(cl)
}

### run
source("CommonGrid.R")
source("PODexpansion.R")
source("CoKriging.R")
source("Emulation.R")
source("Validation.R")

### turn off parallel
if(dompi){
  closecluster(cl)
}else{
  sfStop()
}
