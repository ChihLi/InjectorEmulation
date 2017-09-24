##########################################################
#              *** VALIDATION ***                        #
#    MADE BY CHIH-LI and SIMON ON 10/5/2015              #
#    NOTE: COMPARISON BY FIGURES AND MRE                 #
##########################################################

cat("Validating new design ...\n")

###     Step 0: setting
if(!dir.exists(paste0(output.folder.path, "/validation"))) # for storing the validation results
  dir.create(paste0(output.folder.path, "/validation"))
if(!dir.exists(paste0(output.folder.path, "/validation/comparison"))) # for storing the validation results
  dir.create(paste0(output.folder.path, "/validation/comparison"))
for(r in 1:R){
  if(!dir.exists(paste0(output.folder.path, "/validation/", Response.ColumnIndex[r])))
    dir.create(paste0(output.folder.path, "/validation/", Response.ColumnIndex[r]))
}

###     Step 1: interpolate the original flow data onto the spatial grid of the common geometry
###             so that we can do point-by-point comparison
###             See CommonGrid.R
Coordinate.new <- vector("list", n.newx)
Coordinate.new.index <- vector("list", n.newx)
for (i in 1:n.newx){
  Coordinate.new[[i]] <- read.table(paste0(Testdata.folder.path, "/", i, "/1.dat"), header = TRUE)[,1:2]
  # Cut X-axis by Lval + 4 * Rval and Y-axis by 2 * Rval: See Figure A.1
  select.fg <- (Coordinate.new[[i]][,1] <= (Lval.new[i] + 4 * Rval.new[i])) & (Coordinate.new[[i]][,2] <= (2 * Rval.new[i]))
  Coordinate.new[[i]] <- Coordinate.new[[i]][select.fg, ]
  Coordinate.new.index[[i]] <- which(select.fg)
}
Ytest.ls <- vector("list", n.newx)
cat("   interpolate validation run ")
for (i in 1:n.newx){
  cat(i, ", ")
  RealCoord.i <- Coordinate.new[[i]]
  CommonCoord.i <- CommonCoord.new[[i]]
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
  Ytest.ls[[i]] <- vector("list", R)
  for(r in 1:R){
    response.index <- Response.ColumnIndex[r]
    Ytest.ls[[i]][[r]] <- big.matrix(nrow(CommonCoord.i), length(time.step),
                                 backingpath = paste0(output.folder.path, "/validation/", response.index),
                                 backingfile = paste0("Ytest_", i,".bck"),
                                 descriptorfile = paste0("Ytest_", i,".dsc"))
  }
  Ydes <- lapply(Ytest.ls[[i]],  describe)
  
  foreach.out <- foreach(tm = time.step, .packages = "bigmemory") %dopar% {
    #Read snapshot for case i at timestep tm
    flow <- read.table(paste0(Testdata.folder.path, "/", i, "/", tm, ".dat"), header = TRUE)
    flow <- flow[Coordinate.new.index[[i]], Response.ColumnIndex]
    #interpolation: inverse distance weighting with nearest k = 10 neighbours
    for(r in 1:R){
      Y <- attach.big.matrix(Ydes[[r]], backingpath = paste0(output.folder.path, "/validation/", Response.ColumnIndex[r]))
      interpolate.val <- rep(0, nrow(min.dist.mx))
      for(num in 1:ncol(min.dist.mx)){
        interpolate.val <- interpolate.val + flow[min.dist.mx[,num], r] * weight.mx[,num]
      }
      Y[, tm] <- interpolate.val
    }
  }
}

cat("\n   anaylysis for temperature...\n")
r <- 5  # temperature

for(i in n.newx){
  cat("i = ", i, "\n")
  inlet.fg <- CommonCoord.new[[i]][,1] < (DLval.new + Delta.new/2) & CommonCoord.new[[i]][,1] > (DLval.new - Delta.new/2) #inlet
  transition.fg <- CommonCoord.new[[i]][,1] > (DLval.new + Rval.new) & CommonCoord.new[[i]][,1] < Lval.new  & CommonCoord.new[[i]][,2] < Rval.new  & CommonCoord.new[[i]][,2] > Rval.new*0.75 # fluid transition region
  exit.fg <- CommonCoord.new[[i]][,1] < (Lval.new + Rval.new * 0.25) & CommonCoord.new[[i]][,1] > Lval.new  & CommonCoord.new[[i]][,2] < Rval.new*1.25  #injector exit
  
  for(tt in time.step){
    cat("   t = ", tt, "\n")
    Y.hat <- read.table(paste0(output.folder.path, "/prediction/newX", i, "/pred_", tt,'.dat'), header = TRUE)
    Y.true <- attach.big.matrix(paste0(output.folder.path, "/validation/", Response.ColumnIndex[r], "/Ytest_", i, ".dsc"))
    ### Compute MRE at injector inlet
    yhat.S <- Y.hat[inlet.fg,r+2]
    y.S <- Y.true[inlet.fg, tt]
    cat("   MRE at injector inlet =", round(sum(abs(yhat.S - y.S))/sum(abs(y.S)) * 100, 2), "%\n")
    
    ### Compute MRE at fluid transition region
    yhat.S <- Y.hat[transition.fg,r+2]
    y.S <- Y.true[transition.fg, tt]
    cat("   MRE at fluid transition region =", round(sum(abs(yhat.S - y.S))/sum(abs(y.S)) * 100, 2), "%\n")
    
    ### Compute MRE at injector exit
    yhat.S <- Y.hat[exit.fg,r+2]
    y.S <- Y.true[exit.fg, tt]
    cat("   MRE at injector exit =", round(sum(abs(yhat.S - y.S))/sum(abs(y.S)) * 100, 2), "%\n")
  }
}

cat("\n   plot comparison figures for temperature...\n")

numCol <- 1000
colors <- matlab.like(numCol+1)

for(i in n.newx){
  for(tt in time.step){
    Range <- c(20, 280)
    png(paste0(output.folder.path, "/validation/comparison/i", i, "-t", tt, "-r", r, ".png"), width = 600, height = 400)
    Y.hat <- read.table(paste0(output.folder.path, "/prediction/newX", i, "/pred_", tt,'.dat'), header = TRUE)
    if(r == 4) Y.hat[,r+2] <- Y.hat[,r+2]*1e3   ## unit is different when r = pressure
    zcolor <- rep(0, nrow(Y.hat))
    zcolor[Y.hat[,r+2] < Range[1]] <-  colors[1]
    zcolor[Y.hat[,r+2] > Range[2]] <-  colors[numCol+1]
    select.fg <- Y.hat[,r+2] > Range[1] & Y.hat[,r+2] < Range[2]
    zcolor[select.fg] <- colors[(Y.hat[select.fg,r+2] - Range[1])/diff(Range)*numCol + 1] 
    Coordinate <- Y.hat[,1:2]
    plot(Coordinate[,1], Coordinate[,2], col = zcolor, pch=15, cex=1.2,
         xlab = 'x (m)', ylab = 'y (m)', main = paste0("Emulated (upper) and Simulated (lower) flow i", i, "-t", tt, "-r", r), 
         ylim = c(-max(Coordinate[,2]), max(Coordinate[,2]))) 
    
    Y.true <- attach.big.matrix(paste0(output.folder.path, "/validation/", Response.ColumnIndex[r], "/Ytest_", i, ".dsc"))
    zcolor <- rep(0, nrow(Y.true))
    zcolor[Y.true[,tt] < Range[1]] <-  colors[1]
    zcolor[Y.true[,tt] > Range[2]] <-  colors[numCol+1]
    select.fg <- Y.true[,tt] > Range[1] & Y.true[,tt] < Range[2]
    zcolor[select.fg] <- colors[(Y.true[select.fg,tt] - Range[1])/diff(Range)*numCol + 1] 
    points(Coordinate[,1], -Coordinate[,2], col = zcolor, pch=15, cex=1.2,
           xlab = 'x (m)', ylab = 'y (m)') 
    colorlegend(col = colors, zlim=range(Range))
    dev.off()
  }
}

