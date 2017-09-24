##########################################################
#              *** EMULATION ***                         #
#    MADE BY CHIH-LI and SIMON ON 10/5/2015              #
#    NOTE: PLEASE SEE EQUATION (9) AND (10)              #
##########################################################

cat("Emulating new design ...\n")

###     Step 0: setting
if(!dir.exists(paste0(output.folder.path, "/prediction"))) # for storing the prediction results
  dir.create(paste0(output.folder.path, "/prediction"))
des.new <- read.csv(new.design.file.path)        # read in new design points
Lval.new  <- des.new[,1]/1000                    # length of injector
Rval.new  <- des.new[,2]/1000                    # radius of injector
Delta.new <- des.new[,4]/1000                    # delta of injector
DLval.new <- des.new[,5]/1000                    # dL of injector
des.new.norm <- as.matrix(des.new)
for(i in 1:p) des.new.norm[,i] <- (des.new[,i] - min(des[,i])) / diff(range(des[,i]))  # normalize parameters
n.newx <- nrow(des.new.norm)                 # number of new cases
for(i in 1:n.newx){
  if(!dir.exists(paste0(output.folder.path, "/prediction/newX", i)))
    dir.create(paste0(output.folder.path, "/prediction/newX", i))
}

###     Step 1: create the coordinate for new X
###             See Appendix A.1.3.
CommonCoord.new <- vector("list", n.newx)
for (i in 1:n.newx){
  tmpList <- vector("list", 4)
  #part (a)
  tmpList[[1]] <- cbind(DLval.new[i]*a.coord[,1]/DLval[densest_case], Rval.new[i]*a.coord[,2]/Rval[densest_case])
  #part (b)
  tmpList[[2]] <- cbind(DLval.new[i] + (Lval.new[i] - DLval.new[i])* (b.coord[,1] - DLval[densest_case])/(Lval[densest_case] - DLval[densest_case]), Rval.new[i]*b.coord[,2]/Rval[densest_case])
  #part (c)
  tmpList[[3]] = cbind(Lval.new[i] + (d.coord[,1]-Lval[densest_case]) * Rval.new[i]/Rval[densest_case], Rval.new[i]*d.coord[,2]/Rval[densest_case]) 
  #part (d)
  tmpList[[4]] = cbind(Lval.new[i] + (c.coord[,1]-Lval[densest_case]) * Rval.new[i]/Rval[densest_case], Rval.new[i] + (c.coord[,2]-Rval[densest_case]) * Rval.new[i]/Rval[densest_case]) 
  CommonCoord.new[[i]] <- rbind(tmpList[[1]], tmpList[[2]], tmpList[[3]], tmpList[[4]])
}

###     Step 2: Emulation
###             See equation (9) and (10)

foreach(tt = time.step, .packages = c("Rcpp", "InjectorEmulationHelper")) %dopar% {
  
  # load estimation results
  load(paste0(output.folder.path, "/estimation/allparams_t", tt, ".RData"))
  
  # predict POD mode coefficients 
  predcoef <- GPpred(des.new.norm, des.norm, params$yMat, params$corr, params$mu, as.matrix(params$Tcov), 2)  
  beta.hat <- predcoef$pred  # mean:     equation (8)
  beta.var <- predcoef$var   # variance: equation (8)
  
  for (i in 1:n.newx){
    # reconstruct prediction results to flow field prediction equation (9) and (10)
    Y.hat <- CI.width <- matrix(0, nrow = J, ncol = R)
    colnames(Y.hat) <- colnames(CI.width) <- Response.ColumnNames
    for (r in 1:R){
      load(paste0(output.folder.path, "/POD_expansion/", Response.ColumnIndex[r], "/mode.phi.RData"))
      Y.hat[, r]    <- Y.mean[, r] + c(mode.phi %*% beta.hat[i,(K.cum[r]+1):K.cum[r+1]])
      CI.width[, r] <- qnorm(1-(1-confidence.level)/2) * sqrt( c(mode.phi^2 %*% (beta.var[i, (K.cum[r]+1):K.cum[r+1]])))
    }
    Y.hat    <- cbind(CommonCoord.new[[i]], Y.hat)      # combine with coordinates
    CI.width <- cbind(CommonCoord.new[[i]], CI.width)   # combine with coordinates
    colnames(Y.hat)[1:2] <- colnames(CI.width)[1:2] <- c("x", "y")
    # write predictions and UQ
    write.table(Y.hat,    file = paste0(output.folder.path, "/prediction/newX", i, "/pred_", tt, ".dat"), row.names = FALSE)
    write.table(CI.width, file = paste0(output.folder.path, "/prediction/newX", i, "/uq_",   tt, ".dat"), row.names = FALSE)
  }
}
