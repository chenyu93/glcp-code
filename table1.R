rm(list = ls())
set.seed(0)
library(foreach)
library(doMC)
registerDoMC(4)
source('utils.R')
source('call_functions.R')

options(warn=1)

all_simulate_data = data.frame()

BETA = 12
MU_BASE = 5

for (WIDTH in c(1, 2)){
  for (COV_DECAY in c(0.2, 0.5, 0.8)){
    for (M in c(2, 3, 4)){ 
      for (N in c(2, 3)){
        J = rep(1:M, N)
        cat(WIDTH, COV_DECAY,M,N,'\n')
        region = owin(c(0, WIDTH),c(0,WIDTH))

        
        R = matrix(0,nrow = max(J), ncol = max(J))
        for (i in 1:max(J)){
          for (j in 1:max(J)){
            R[i,j] = COV_DECAY ^ abs(i - j)
          }
        }
        rvals = seq(0, 0.15, 0.15/20)
        
        collection = foreach(i=1:100, .combine=rbind) %dopar%{
          result = wrapper.grouping(MU_BASE, R, BETA, region, J, COV_DECAY)
          result
        }
        temp = report_table1(collection)
        cat(temp, '\n')
        all_simulate_data = rbind(all_simulate_data, temp)
      }
    }
  }
}


colnames(all_simulate_data) = c('M', '|D|', 'rho', 'N', 'E(X)', 'ufp', 'cfp', 'ofp', 'tpgp', 'tngp')


print(xtable(all_simulate_data[order(all_simulate_data$M, all_simulate_data$rho),c(1, 3, 2, 4, 5,6,7,8,9,10, 11)],digits=c(0, 0 ,1,0,0,2,2,2,2,2,2, 1)), include.rownames = F)


# tstart = Sys.time()
# grouping_result = calc.grouping.result(X, rvals)
# grouping_result$time = as.numeric(Sys.time() - tstart)


