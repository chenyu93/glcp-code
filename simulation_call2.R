
rm(list = ls())
set.seed(1)
library(foreach)
library(doMC)
registerDoMC(4)
source('utils.R')
source('call_functions.R')


all_simulate_data = data.frame()

BETA = 12
MU_BASE = 5
COV_DECAY = 0.5

for (M in c(1,4)){
  for (COV_DECAY in c(0.2, 0.5, 0.8)){
    for (WIDTH in c(1, 2)){
      for (N in 1:3){
        J = rep(1:M, N)
        region = owin(c(0, WIDTH),c(0,WIDTH))
        R = matrix(0,nrow = max(J), ncol = max(J))
        for (i in 1:max(J)){
          for (j in 1:max(J)){
            R[i,j] = COV_DECAY ^ abs(i - j)
          }
        }
        rvals = seq(0, 0.15, 0.15/20)
        
        collection = foreach(i=1:100, .combine=rbind) %dopar%{
          result = wrapper.params.est(MU_BASE, R, BETA, region, J, COV_DECAY)
          result
        }
        temp = report_table2(collection)
        all_simulate_data = rbind(all_simulate_data, temp)
        cat(temp, '\n')
      }
    }
  }
}

colnames(all_simulate_data) = c('M', '|D|', 'rho', 'N', 'E(X)', 'beta', 'mu', 'sigma', 'omega', 'time')

all_simulate_data = all_simulate_data[c(1,3,2,4,5,6,7,8,9)]
library(xtable)

print(xtable(all_simulate_data,digits=c(0, 0,1,0,0,2,2,2,2,2)), include.rownames = F)

# tstart = Sys.time()
# grouping_result = calc.grouping.result(X, rvals)
# grouping_result$time = as.numeric(Sys.time() - tstart)






