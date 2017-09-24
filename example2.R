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
Testdata.folder.path <- "data2/testdata"
design.file.path     <- "data2/DoE.csv"
new.design.file.path <- "data2/newX.csv"
output.folder.path   <- "output2"
time.step            <- 1:3          # represent the real data at time 725, 775, and 825
Response.ColumnIndex <- 4:9          # indicate which columns present responses in the raw data     
Response.ColumnNames <- c("u", "v", "w", "P", "T", "rho")  
POD.energy           <- 0.99         # the proportion of total flow energy for POD truncation
confidence.level     <- 0.80         # indicate the level of confidence interval
dompi                <- FALSE        # do mpi parallel or not  

### the global objects for skipping Common Grid.R and PODexpansion.R ###
des <- read.csv(design.file.path)        # read in design points
p <- ncol(des)                           # number of parameters
n <- nrow(des)                           # number of cases
TT <- length(time.step)                  # number of time steps 
R <- length(Response.ColumnIndex)        # number of responses
N  <-  n * TT
Lval  <- des[,1]/1000                    # length of injector
Rval  <- des[,2]/1000                    # radius of injector
DLval <- des[,5]/1000                    # dL of injector
densest_case <- 22
K <- c(47, 141, 45, 7, 655, 142)
J <- 84544
Number.NN <- 10
load(paste0(output.folder.path, "/CommonGrid.RData"))
load(paste0(output.folder.path, "/MeanFlow.RData"))
##########################################################################

### set up for parallel
#sfInit(parallel = TRUE, cpus = 20, type = "MPI")
if(dompi){
  cl <- startMPIcluster(count = 2, verbose = TRUE, maxcores = 20)
  registerDoMPI(cl)
}else{
  sfInit(parallel = TRUE, cpus = 64)
  cl <- sfGetCluster()
  registerDoParallel(cl)
}

### run
source("CoKriging.R")
source("Emulation.R")
source("Validation.R")
source("plotPOD.R")

### turn off parallel
if(dompi){
  closecluster(cl)
}else{
  sfStop()
}

