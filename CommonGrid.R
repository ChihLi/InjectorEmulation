##########################################################
#         *** PRODUCING COMMON GRID ***                  #
#    MADE BY CHIH-LI and SIMON ON 10/5/2015              #
#    NOTE: PLEASE READ THE APPENDIX A.1 Common Grid      #
##########################################################

cat("Producing Common Grid ...\n")
###     Step 0: setting
des <- read.csv(design.file.path)        # read in design points
p <- ncol(des)                           # number of parameters
n <- nrow(des)                           # number of cases
TT <- length(time.step)                  # number of time steps 
R <- length(Response.ColumnIndex)        # number of responses
N  <-  n * TT
Lval  <- des[,1]/1000                    # length of injector
Rval  <- des[,2]/1000                    # radius of injector
DLval <- des[,5]/1000                    # dL of injector
if(!dir.exists(paste0(output.folder.path))) # for storing the new flow with the common geometry 
  dir.create(paste0(output.folder.path))
if(!dir.exists(paste0(output.folder.path, "/flow_commongrid"))) # for storing the new flow with the common geometry 
  dir.create(paste0(output.folder.path, "/flow_commongrid"))
for(r in 1:R){
  if(!dir.exists(paste0(output.folder.path, "/flow_commongrid/", Response.ColumnIndex[r])))
    dir.create(paste0(output.folder.path, "/flow_commongrid/", Response.ColumnIndex[r]))
}

###     Step 1: Identify the densest grid (i.e., with the most grid points) among the n simulation runs, and set this as the common reference grid.
###                                     See Appendix A.1.1.
Coordinate <- vector("list", n)
Coordinate.index <- vector("list", n)
for (i in 1:n){
  Coordinate[[i]] <- read.table(paste0(Rawdata.folder.path, "/", i, "/1.dat"), header = TRUE)[,1:2]
  # Cut X-axis by Lval + 4 * Rval and Y-axis by 2 * Rval: See Figure A.1
  select.fg <- (Coordinate[[i]][,1] <= (Lval[i] + 4 * Rval[i])) & (Coordinate[[i]][,2] <= (2 * Rval[i]))
  Coordinate[[i]] <- Coordinate[[i]][select.fg, ]
  Coordinate.index[[i]] <- which(select.fg)
}
grid.vt <- sapply(Coordinate, FUN = nrow)
densest_case <- which.max(grid.vt)
J <- max(grid.vt)
cat("   Densest grid is in case", densest_case, "\n")

###     Step 2: For each simulation, partition the grid into the following four parts
###                    See Appendix A.1.2 and Figure A.1.

# Divide into Four parts
CutCoordInd <- vector("list",n)
for (i in 1:n){
  tmpList <- vector("list",4)
  tmpList[[1]] <- which((Coordinate[[i]][,1]<=DLval[i]))
  tmpList[[2]] <- which((Coordinate[[i]][,1]>DLval[i]) & (Coordinate[[i]][,1]<=Lval[i]))
  tmpList[[3]] <- which((Coordinate[[i]][,1]>Lval[i]) & (Coordinate[[i]][,2]<=Rval[i]))
  tmpList[[4]] <- which((Coordinate[[i]][,1]>Lval[i]) & (Coordinate[[i]][,2]>Rval[i]))

  CutCoordInd[[i]] <- tmpList
}

# Four parts of injector with densest grid
a.coord <- Coordinate[[densest_case]][CutCoordInd[[densest_case]][[1]],]
b.coord <- Coordinate[[densest_case]][CutCoordInd[[densest_case]][[2]],]
d.coord <- Coordinate[[densest_case]][CutCoordInd[[densest_case]][[3]],]
c.coord <- Coordinate[[densest_case]][CutCoordInd[[densest_case]][[4]],]
save(a.coord, b.coord, c.coord, d.coord, des, file = paste0(output.folder.path, "/CommonCoord.RData"))

###     Step 3: Linearly rescale each part of the partition to the common grid by corresponding geometry parameters
###                        See Appendix A.1.3.

CommonCoord <- vector("list", n)
for (i in 1:n){
  tmpList <- vector("list", 4)
  #part (a)
  tmpList[[1]] <- cbind(DLval[i]*a.coord[,1]/DLval[densest_case], Rval[i]*a.coord[,2]/Rval[densest_case])
  #part (b)
  tmpList[[2]] <- cbind(DLval[i] + (Lval[i] - DLval[i])* (b.coord[,1] - DLval[densest_case])/(Lval[densest_case] - DLval[densest_case]), Rval[i]*b.coord[,2]/Rval[densest_case])
  #part (c)
  tmpList[[3]] = cbind(Lval[i] + (d.coord[,1]-Lval[densest_case]) * Rval[i]/Rval[densest_case], Rval[i]*d.coord[,2]/Rval[densest_case]) 
  #part (d)
  tmpList[[4]] = cbind(Lval[i] + (c.coord[,1]-Lval[densest_case]) * Rval[i]/Rval[densest_case], Rval[i] + (c.coord[,2]-Rval[densest_case]) * Rval[i]/Rval[densest_case]) 
  CommonCoord[[i]] <- rbind(tmpList[[1]], tmpList[[2]], tmpList[[3]], tmpList[[4]])
}

###     Step 4:For each simulation, interpolate the original flow data onto the spatial grid of the common geometry
###                        See Appendix A.1.4.

Number.NN <- 10 # number of nearest neighbours for inverse distance weighting interpolation method
Y.ls <- vector("list", n)
Y.mean <- matrix(0, nrow = nrow(CommonCoord[[densest_case]]), ncol = R) # for POD expansion use
cat("   Interpolate the original flow data onto the spatial grid: Case ")
for (i in 1:n){
  cat(i, ", ")
  RealCoord.i <- Coordinate[[i]]
  CommonCoord.i <- CommonCoord[[i]]
  coords <- t(RealCoord.i[,1:2])
  out.Matrix <- foreach(k = 1:nrow(CommonCoord.i), .combine = rbind) %dopar% {
    sort.out <- sort.int(colSums((coords - CommonCoord.i[k,])^2),
                         decreasing = FALSE, index.return = TRUE)
    if(sort.out$x[1] < 10^(-16)){ #if exactly the same point, weight 1 o.w 0
      min.dist.ind <- sort.out$ix[1:Number.NN]
      output <- c(min.dist.ind, c(1, rep(0, Number.NN-1)))
    }else{
      min.dist <- sort.out$x[1:Number.NN]
      weight <- 1/min.dist
      weight <- weight/sum(weight)
      min.dist.ind <- sort.out$ix[1:Number.NN]
      output <- c(min.dist.ind, weight)
    }
    return(output)
  }
  min.dist.mx <- out.Matrix[,1:Number.NN]
  weight.mx <- out.Matrix[,(Number.NN + 1) : (2 * Number.NN)]

  rownames(min.dist.mx) <- NULL
  rownames(weight.mx) <- NULL
  
  # use bigmemory to store the flow with common grid
  Y.ls[[i]] <- vector("list", R)
  for(r in 1:R){
    response.index <- Response.ColumnIndex[r]
    Y.ls[[i]][[r]] <- big.matrix(nrow(CommonCoord.i), length(time.step),
                             backingpath = paste0(output.folder.path, "/flow_commongrid/", response.index),
                             backingfile = paste0("Y_", i,".bck"),
                             descriptorfile = paste0("Y_", i,".dsc"))
  }
  Ydes <- lapply(Y.ls[[i]],  describe)
  
  foreach.out <- foreach(tm = time.step, .packages = "bigmemory") %dopar% {
    #Read snapshot for case i at timestep tm
    flow <- read.table(paste0(Rawdata.folder.path, "/", i, "/", tm, ".dat"), header = TRUE)
    flow <- flow[Coordinate.index[[i]],Response.ColumnIndex]
    #interpolation: inverse distance weighting with nearest k = 10 neighbours
    for(r in 1:R){
      Y <- attach.big.matrix(Ydes[[r]], backingpath = paste0(output.folder.path, "/flow_commongrid/", Response.ColumnIndex[r]))
      interpolate.val <- rep(0, nrow(min.dist.mx))
      for(num in 1:ncol(min.dist.mx)){
        interpolate.val <- interpolate.val + flow[min.dist.mx[,num], r] * weight.mx[,num]
      }
      Y[, tm] <- interpolate.val
    }
  }
  
  Y.mean <- Y.mean + 
    foreach(r = 1:R, .packages = c("bigmemory", "biganalytics"), .combine = cbind) %dopar% {
      Y <- attach.big.matrix(Ydes[[r]], backingpath = paste0(output.folder.path, "/flow_commongrid/", Response.ColumnIndex[r]))
      apply(Y, 1, mean)/n
    }
}

# save mean values of Y
save(Y.mean, file = paste0(output.folder.path, "/MeanFlow.RData"))


# Do centering for POD expansion use
for(r in 1:R){
  foreach(i = 1:n, .packages = "bigmemory") %dopar% {
    ##### Compute mean at each grid for each case
    Y <- attach.big.matrix(paste0(output.folder.path, "/flow_commongrid/", Response.ColumnIndex[r], "/Y_", i, ".dsc"))
    Y[] <- Y[] - Y.mean[,r]
  }
}

cat("\n")

