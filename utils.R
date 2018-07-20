# install.packages(c("spatstat", "foreach", "RandomFields", "stringr", "doMC"))
library(spatstat)
library(spatstat.utils)
library(stats)
library(stringr)
# install.packages("RandomFields")
library(RandomFields)
set.seed(1)

do.rMultiGRF <- function( R, BETA, win) {
  # generate a multi variate GRF on region ``win''
  # Input:
  #     R:     The covariance matrix of Z(s) with size M x M
  #     BETA:  The spatial covariance function coefficient 
  #     win:   The region object, created by spatstat.owin()
  #
  # Return:
  #    Z(s), the desired M dimension GRF
  
  R.eig = eigen(R)
  rho = RMexp(var = 1, scale = 1/BETA)
  model = RMmatrix(M = R.eig$vectors[,1]*sqrt(R.eig$values[1]), rho)
  if (nrow(R) > 1){
    for (i in 2:length(R.eig$values)){
      model = model + RMmatrix(M = R.eig$vectors[,i]*sqrt(R.eig$values[i]), rho)
    }
  }
  w <- as.mask(w=win, eps=NULL, dimyx=20, xy=NULL)
  xcol <- w$xcol;yrow <- w$yrow;dim <- w$dim
  RFoptions(spConform=FALSE)
  z <- RFsimulate(model, xcol, yrow, grid = TRUE)
  RFoptions(spConform=TRUE)
  return(z)
}



do.Intensity <- function(mu,GRF,win){
  # Generate multi variate LGCP with mu and GRF on region win
  #
  # Input:
  #       mu:    A M-dimensional vector
  #       GRF:   A realization of M-dimensional GRF
  #       win:   The region object, created by spatstat.owin()
  #
  # Return:
  #       Lambda(s)  The conditional intensity function
  
  M = length(mu)
  Intensities = list()
  if (M > 1){
  for (i in 1:M){
    Intensities[[i]] =  as.im(exp(mu[i] + GRF[, , i]), W=win)
  }}
  else{
    Intensities[[1]] =  as.im(exp(mu[i] + GRF), W=win)
  }
  return(Intensities)
}


do.Generate.Process <- function(mu, R, BETA, region, J, omega){
  # Generate Point Process
  #
  # Input:
  #       mu:   A M-dimensional vector
  #       R:     The covariance matrix of Z(s) with size M x M
  #       BETA:  The spatial covariance function coefficient 
  #       region:   The region object, created by spatstat.owin()
  #       J:    The mapping function
  #       omega:    The link coefficients
  
  GRF <- do.rMultiGRF(R = R, BETA = BETA, win = region)
  Intensities = do.Intensity(mu=mu, GRF=GRF,win=region)
  ProcessData = data.frame(x=NULL,y=NULL,cate = NULL)
  
  for (i in 1:length(J)){
    tmp = Intensities[[J[i]]]
    X = rpoispp(tmp*omega[i])[region]
    ProcessData = rbind(ProcessData,
                        data.frame(x = X$x,y = X$y,cate = i))}
  X = ppp(ProcessData$x, ProcessData$y, window = region, marks = as.factor(ProcessData$cate))
  return(X)
  }



cache.emp.noparallel <- function(X, rvals){
  cached_KIJ = list()
  total_processs = length(levels(X$marks))
  for (i in 1:total_processs){
    cached_KIJ[[i]] = list()
    for (j in 1:total_processs){
      cached_KIJ[[i]][[j]] = get.empirical(X, i, j, rvals)
    }
  }
  return(cached_KIJ)
}

cache.emp <- function(X, rvals){
  total_process = length(levels(X$marks))
  grp_ind1 = c()
  grp_ind2 = c()
  for (mi in 1:total_process){
    for (mj in 1:total_process){
      grp_ind1 = c(grp_ind1, mi)
      grp_ind2 = c(grp_ind2, mj)
    }
  }
  
  collection = foreach(i=1:length(grp_ind1), .combine=rbind) %dopar%{
    mi = grp_ind1[i]; mj = grp_ind2[i];
    result = list(get.empirical(X, mi, mj, rvals))
    result
  }
  
  cached_KIJ = list()
  for (process_i in 1:length(levels(X$marks))){
    cached_KIJ[[process_i]] = list()
  }
  
  for (i in 1:length(grp_ind1)){
    ii = grp_ind1[i]; jj = grp_ind2[i];
    cached_KIJ[[ii]][[jj]] = unlist(collection[i])
  }
  return(cached_KIJ)
}

get.theoretical <- function (sigma12, rvals, beta){
  # get theoretical kross k   sigma12 * rho(h)
  # Input:
  #     rvals: points to evaluate the K function, e.g., seq(0, 0.1*1, 0.1*1/20)
  #     sigma12:  cross covariance element
  #     beta:   scalar
  #
  # Output:
  #     K function evaluated at rvals
  integrand <- function(r, sigma12, ...) 2 * pi * r * exp(sigma12 * exp(-r * beta))
  
  nr <- length(rvals)
  th <- numeric(nr)
  th[1L] <- if (rvals[1L] == 0) 
    0
  else integrate(integrand, lower = 0, upper = rvals[1L], 
                 sigma12 = sigma12)$value
  for (i in 2:length(rvals)) {
    delta <- integrate(integrand, lower = rvals[i - 1L], 
                       upper = rvals[i], sigma12 = sigma12)
    th[i] = th[i - 1L] + delta$value
  }
  return(th)
}


get.empirical <- function(X, i, j, rvals){
  hatKij = Kcross(X,
                  i = levels(X$marks)[i],
                  j = levels(X$marks)[j],
                  r = rvals,
                  correction = "isotropic")
  # rmax <- quantile(rvals,0.999); rmin <- quantile(rvals,0.001)
  # em <- hatKij$iso
  # sub <- (rvals >= rmin) & (rvals <= rmax)
  # rvals <- rvals[sub]; em <- em[sub] 
  return(hatKij$iso)
}


cal.beta_fitted <- function(X){
  pfm2 = lgcp.estK(X[X$marks==1], rmin=0, rmax=0.15)
  return( 1 / pfm2$par[2])
}



contrast.sum <- function(theta , allargs, ...){
  diff = 0
  objarg = allargs[1]
  objarg = unlist(objarg,recursive = FALSE)
  par = c(var = theta, scale = objarg$betafitted)
  theo <- get.theoretical(sigma12=par[1], rvals = objarg$rvals, beta=par[2])
  for (process_i in 1:length(objarg$Gi)){
    for (process_j in 1:length(objarg$Gj)){  
      objarg = allargs[(process_i-1) * length(objarg$Gj) + process_j]
      objarg = unlist(objarg,recursive = FALSE)
      discrep <- (abs(theo^0.25 - objarg$hatKij^0.25))^2
      diff = diff + mean(discrep)
    }
  }
  return(diff)
}



cal.sigmaij <-  function(X, Gi, Gj, rvals, betafitted){
    # calc the loss of merged Group i & j
    # Input:
    #        X:   marked point process
    #        Gi:  a group, contains list of process ids, length(Gi) = M_i
    #        Gj:  a group, contains list of process ids, length(Gj) = M_j
    #        rvals: evaluate K function at rvals
    #        betafitted: pre-fitted beta
    # 
    # Return:
    #        estimated sigmaij (scalar)
    allargs = list()
    for (process_i in Gi){
      for (process_j in Gj){
        hatKij = get.empirical(X, process_i, process_j, rvals)
        objargs <- list(rvals = rvals,
                        betafitted = betafitted,
                        Gi = Gi, Gj = Gj, 
                        hatKij = hatKij)
        allargs = c(allargs,list(objargs))
        length(allargs) 
      }
    }
    startpar = 1
    minimum <- optim(startpar, 
                     fn = contrast.sum, 
                     allargs = allargs,
                     method = "Brent", 
                     lower = 0.1,
                     upper = 4,
                     control = list(reltol=1e-6, maxit=100));
   return(minimum$par)
}

cal.group.loss <- function(X, Gi, Gj, rvals, betafitted, cached_KIJ){
  # calc the loss of merged Group i & j
  # Input:
  #        X:   marked point process
  #        Gi:  a group, contains list of process ids, length(Gi) = M_i
  #        Gj:  a group, contains list of process ids, length(Gj) = M_j
  #        rvals: evaluate K function at rvals
  #        betafitted: pre-fitted beta
  # 
  # Return:
  #        A losses matrix of size len(Gi), len(Gj), each elements represent the 
  #        loss of using ``linked'' model
  total_processes = length(levels(X$marks))
  Loss = matrix(0, total_processes, total_processes)
  allargs = list()
  for (process_i in Gi){
    for (process_j in Gj){
      hatKij = cached_KIJ[[process_i]][[process_j]]
      objargs <- list(rvals = rvals,
                      betafitted = betafitted,
                      Gi = Gi, Gj = Gj, 
                      hatKij = hatKij)
      allargs = c(allargs,list(objargs))
      length(allargs) 
    }
  }
  startpar = 1
  minimum <- optim(startpar, 
                   fn = contrast.sum, 
                   allargs = allargs,
                   method = "Brent", 
                   lower = 0.1,
                   upper = 4,
                   control = list(reltol=1e-6, maxit=100));
  final_theo = get.theoretical(minimum$par, rvals, betafitted) 
  
  for (process_i in Gi){
    for (process_j in Gj){

      hatKij = cached_KIJ[[process_i]][[process_j]]
      
      discrep <- (abs(final_theo^0.25 - hatKij^0.25))^2
      Loss[process_i,process_j] = mean(discrep)
    }
  }
  return(Loss)
}


mergeGrp <- function(working_G, i, j){
  tmp = working_G[j]
  working_G = working_G[-j]
  working_G[i] = list(sort(c(unlist(working_G[i]),unlist(tmp))))
  return(working_G)
}



measure <- function(X, J, estimators, beta, mu, Sigma, omega){
  sum_squared_mu = 0
  for (i in 1:length(mu)){
    sum_squared_mu = sum_squared_mu + (estimators$mu[i] - mu[i])^2
  }
  
  sum_squared_omega = 0
  for (i in 1:length(omega)){
    sum_squared_omega = sum_squared_omega + (estimators$omega[i] - omega[i])^2
  }
  
  sum_squared_Sigma = 0
  for (i in 1:nrow(Sigma)){
    for (j in 1:ncol(Sigma)){
      sum_squared_Sigma = sum_squared_Sigma + (Sigma[i, j] - estimators$Sigma[i, j])^2
    }
  }
  
  
  result = list(M=max(J),
                size=(X$window$xrange[2] * X$window$yrange[2]),
                N=length(J),
                n=X$n, 
                rmse_beta=(estimators$beta - beta)^2,
                rmse_mu=sum_squared_mu / length(mu),
                rmse_Sigma=sum_squared_Sigma / nrow(Sigma) / ncol(Sigma),
                rmse_omega=sum_squared_omega / length(omega)
                )
}



report_table1 <- function(data){
  data = data.frame(data)
  result = c(
    M = unlist(data$M[1]),
    size = unlist(data$size[1]),
    rho = unlist(data$rho[1]),
    N = unlist(data$N[1]),
    n = mean(unlist(data$n)),
    ufp = mean(unlist(data$ufp)),
    cfp = mean(unlist(data$cfp)),
    ofp = mean(unlist(data$ofp)),
    tpgp = mean(unlist(data$tpgp)),
    tngp = mean(unlist(data$tngp)),
    time = median(unlist(data$time))
  )
  return(result)
}

report_table2 <- function(data){
  data = data.frame(data)
  result = c(
    M = unlist(data$M[1]),
    size = unlist(data$size[1]),
    rho = unlist(data$COV_DECAY[1]),
    N = unlist(data$N[1]),
    n = mean(unlist(data$n)),
    beta = sqrt(mean(unlist(data$rmse_beta))),
    mu = sqrt(mean(unlist(data$rmse_mu))),
    sigma = sqrt(mean(unlist(data$rmse_Sigma))),
    omega = sqrt(mean(unlist(data$rmse_omega))),
    time = median(unlist(data$time))
  )
  return(result)
}

