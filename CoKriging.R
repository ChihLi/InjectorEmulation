##########################################################
#                *** COKRIGING ***                       #
#    MADE BY CHIH-LI and SIMON ON 10/5/2015              #
#    NOTE: PLEASE READ THE ALGORITHM 1 IN PAPER          #
##########################################################

cat("Fitting Kriging Model ...\n")
###     Step 0: setting
if(!dir.exists(paste0(output.folder.path, "/estimation"))) # for storing the estimation results
  dir.create(paste0(output.folder.path, "/estimation"))
des.norm <- apply(des, 2, function(x) (x - min(x))/diff(range(x)) )    # normalize each column of the design parameters into [0,1]
B.mx <- matrix(0, nrow = N, ncol = sum(K)) # matrix of POD coefficients
K.cum <- c(0, cumsum(K))
for (r in 1:R){
  load(paste0(output.folder.path, "/POD_expansion/", Response.ColumnIndex[r], "/mode.B.RData"))
  B.mx[, (K.cum[r]+1):K.cum[r+1]] <- mode.B
}

foreach(tt = time.step, .packages = c("Rcpp", "huge", "SLHD", "nloptr", "InjectorEmulationHelper")) %dopar% {
  
  select.index <- seq(from = tt, by = TT, length.out = n)  # find indices for each case
  B <- B.mx[select.index,]

  ###     Step 1: BCD algorithm for maximum likelihood estimation: Optimizing mu and correlation parameters
  ###                  See Algorithm 1 in paper
  # lower bound and upper bound for correlation parameters
  lower <- 0.01
  upper <- 0.99
  # the estimates for the correlation parameters
  corr <- lbfgs(rep(0.5, p), fn = GPdev, control = list(xtol_rel = 1e-2),
                          lower = rep(lower, p), upper = rep(upper, p), xMat = des.norm, yMat = B, pw = 2)$par
  pm <- computeS(corr, des.norm, B, 2);   # pm$mu = mu, pm$Rinv = Rinv, pm$Rsd = data matrix for T
  vars <- apply(pm$Rsd, 2, function(xx){mean(xx^2)})
  
  ###     Step 2: BCD algorithm for maximum likelihood estimation: Optimizing T
  ###             See Algorithm 1 in paper
  glhuge <- huge(x = pm$Rsd, method = 'glasso', lambda.min.ratio = 0.1)
  glsel  <- huge.select(glhuge)
  optind <- which.min(glsel$ebic.score[-1]) + 1 #Find the best model which gives correlations
  invMat <- cov2cor(as.matrix(glsel$icov[[optind]]))
  Tinv   <- cor2cov( (invMat + t(invMat))/2, 1/sqrt(vars*diag(glsel$icov[[optind]])) )
  Tcov   <- solve(Tinv)
  
  # Identify which correlations are significant
  testmat <- (Tinv)
  idx <- which(testmat!=0, arr.ind=T)
  offdiag <- apply(idx,1,function(ind){
    if (ind[1]==ind[2]){
      return (0)
    }
    else{
      return (1)
    }
  })
  idxnz <- idx[which(offdiag==1),]
  
  params <- list(mu = pm$mu, Ticov = Tinv, Tcov = Tcov, corr = corr, yMat = B, idxnz = idxnz)
  save(params, file = paste0(output.folder.path, "/estimation/allparams_t", tt, ".RData"))
}
