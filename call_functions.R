source('utils.R')


calc.parameters.est <- function(X, rvals, J, beta=-1){
  # calculate parameters of process X, assume mapping J() is known
  #  Input:
  #       X: The spatial Point Process
  #       rvals: points for K-function
  #       J:  mapping scheme
  #
  ind_name = levels(X$marks)
  G = list()
  for (i in 1:max(J)){
    G[i] = list((1:length(J))[J==i])
  } 
  if (beta == -1){
    betas = c()
    for (group in 1:max(J)){
      process_id = min(unlist(G[group]))
      tmp = lgcp.estK(X[X$marks==process_id], rmin=0, rmax=0.15)
      betas = c(betas,  1 / tmp$par[2])
    }
  }else{
    betas = c(beta)
  }
  sigma_hat = matrix(0, length(G), length(G))
  for (i in 1:length(G)){
    for (j in i:length(G)){
      sigma_hat[i,j] = cal.sigmaij(X, unlist(G[i]), unlist(G[j]), rvals, mean(betas)) 
    }
  }
  if (length(G) > 1){
    for (i in 2:length(G)){
      for (j in 1:(i-1)){
        sigma_hat[i,j] = sigma_hat[j, i]
      }
    }
  }
  lambda_hat = rep(0, length(ind_name))
  for (process in 1:length(J)){
    lambda_hat[process] =  X[X$marks==ind_name[process]]$n / (X$window$xrange[2] * X$window$yrange[2]) 
  }
  mu_hat = c()
  for (group in 1:length(G)){
    mu_hat = c(mu_hat, log(max(lambda_hat[unlist(G[group])])) - sigma_hat[group, group] / 2 )
  }
  
  oemga_h = rep(0, length(ind_name))
  for (group in G){
    for (process_id in unlist(group)){
      oemga_h[process_id] =  lambda_hat[process_id] / max(lambda_hat[unlist(group)])
    }
  }
  
  return(list(beta=mean(betas), mu=mu_hat, Sigma=sigma_hat, omega=oemga_h))
}
  


wrapper.params.est <- function(MU_BASE, R, BETA, region, J, COV_DECAY){
  mu = MU_BASE + runif(max(J))
  omega = runif(n=length(J),min = 0.5,max = 1)
  G = list()
  for (i in 1:max(J)){
    G[i] = list((1:length(J))[J==i])
  } 
  for (group in G){
    omega[unlist(group)[1]] = 1
  }
  X = do.Generate.Process(mu=mu, R=R, BETA=BETA, region=region, J=J, omega=omega)
  tstart = Sys.time()
  estimators = calc.parameters.est(X, rvals, J)
  performance = measure(X, J,estimators, BETA, mu, R, omega)
  performance$time = as.numeric(Sys.time() - tstart, units="secs")
  performance$COV_DECAY = COV_DECAY
  return(performance)
}



calc.grouping.result <- function(X, rvals){
  #  Generate merging paths
  #  Input:
  #       X:  The spatial point Process
  #       rvals:  points for K-function
  #
  #  Return:
  #       A structure, c(merge_path=merge_path, loss_curve=loss_curve)
  #
  cached_KIJ = cache.emp.noparallel(X, rvals)
  betafitted = cal.beta_fitted(X)
  num_process = length(levels(X$marks))
  working_G = list()
  for (i in 1:num_process){
    working_G = c(working_G,list(i))
  }
  
  Loss = matrix(0,num_process,num_process)
  for (gi in 1:(length(working_G)-1)){
    for (gj in (gi + 1):length(working_G)){
      Loss = Loss + cal.group.loss(X,unlist(working_G[gi]),unlist(working_G[gj]), rvals, betafitted, cached_KIJ=cached_KIJ)
    }
  }
  
  Loss0 = Loss
  
  merge_path = list()
  rep_str = ""
  for (i in working_G){
    rep_str = str_c(rep_str,i,sep=',')
  }
  merge_path[[1]] = rep_str
  loss_curve = c()
  loss_curve = c(loss_curve, -log(sum(Loss)))
  
  for (loop in 1:(num_process -1)){
    min_scalar_loss = 9999
    min_matrix_loss = matrix(0,num_process,num_process)
    for (mi in 1:(length(working_G)-1)){
      for (mj in (mi+1):length(working_G)){
        new_G = mergeGrp(working_G,mi,mj)
        Loss = Loss0
        
        is_change = matrix(0,num_process,num_process)
        for (ii in unlist(working_G[mi])){
          is_change[ii,] = 1
          is_change[,ii] = 1
        }
        for (ii in unlist(working_G[mj])){
          is_change[ii,] = 1
          is_change[,ii] = 1
        }
        
        for (gi in 1:length(new_G)){
          for (gj in gi:length(new_G)){
            flag = FALSE
            for (ii in unlist(new_G[gi])){
              for (jj in unlist(new_G[gj])){
                if (is_change[ii,jj] == 1){
                  flag = TRUE
                }
              }
            }
            if (flag){
              for (ii in unlist(new_G[gi])){
                for (jj in unlist(new_G[gj])){
                  Loss[ii,jj] = 0
                }
              }
              Loss = Loss + cal.group.loss(X, unlist(new_G[gi]),unlist(new_G[gj]), rvals, betafitted, cached_KIJ=cached_KIJ)
            }
          }
        }
        scalar_loss = sum(Loss - Loss0)
        if (scalar_loss < min_scalar_loss){
          min_scalar_loss = scalar_loss
          MI = mi
          MJ = mj
          min_matrix_loss = Loss
        }
      }
    }
    working_G = mergeGrp(working_G ,MI,MJ)
    rep_str = c()
    for (i in working_G){
      rep_str = c(rep_str,paste(i,collapse = '-'))
    }
    merge_path[[loop]] = paste(rep_str, collapse = ',')
    Loss0 = min_matrix_loss
    loss_curve = c(loss_curve, - log(sum(min_matrix_loss)))
  }
  return(list(merge_path=merge_path, loss_curve=loss_curve))
}



wrapper.grouping <- function(MU_BASE, R, BETA, region, J, COV_DECAY){
  mu = MU_BASE + runif(max(J))
  omega = runif(n=length(J),min = 0.5,max = 1)
  G = list()
  for (i in 1:max(J)){
    G[i] = list((1:length(J))[J==i])
  } 
  for (group in G){
    omega[unlist(group)[1]] = 1
  }
  X = do.Generate.Process(mu=mu, R=R, BETA=BETA, region=region, J=J, omega=omega)
  tstart = Sys.time()
  grouping_result = calc.grouping.result(X, rvals)
  
  tau = c()
  for (i in 2:(length(grouping_result$loss_curve)-1)){
    tau = c(tau, (grouping_result$loss_curve[i+1] - grouping_result$loss_curve[i]) / (grouping_result$loss_curve[i] - grouping_result$loss_curve[i-1]))
  }
  M_hat = length(J) - which.max(tau)
  fitting_probability = c(0,0,0)
  if (M_hat < M) fitting_probability[1] =1
  if (M_hat == M) fitting_probability[2] = 1
  if (M_hat > M) fitting_probability[3] = 1
  
  
  grouping_schema = unlist(grouping_result$merge_path[length(J) - M_hat])
  temp = strsplit(grouping_schema,',')
  grouping_schema_str = unlist(temp[order(temp)])
  J_hat = c()
  
  for (i in 1:length(grouping_schema_str)){
    est_group = grouping_schema_str[i]
    temp = unlist(strsplit(est_group,'-'))
    for (j in temp){
      J_hat[as.numeric(j)] = i
    }
  }
  
  
  # calculate true positive
  right=0; total=0;
  for (process_i in 1:length(J)){
    for (process_j in 1:length(J)){
      if (J[process_i] == J[process_j]){
        total = total + 1
        if (J_hat[process_i] == J_hat[process_j]){
          right = right + 1
        }
      }
    }
  }
  tpgp = right / total
  
  # calculate true negative
  
  right=0; total=0;
  for (process_i in 1:length(J)){
    for (process_j in 1:length(J)){
      if (J[process_i] != J[process_j]){
        total = total + 1
        if (J_hat[process_i] != J_hat[process_j]){
          right = right + 1
        }
      }
    }
  }
  tngp = right / total
  
  performance = c(M=max(J), 
                 size=(X$window$xrange[2] * X$window$yrange[2]),
                 rho = COV_DECAY,
                 N=length(J),
                 n=X$n, 
                 ufp=fitting_probability[1],
                 cfp=fitting_probability[2],
                 ofp=fitting_probability[3],
                 tpgp=tpgp,
                 tngp=tngp,
                 time=as.numeric(Sys.time() - tstart, units="secs"))
  return(performance)
}



