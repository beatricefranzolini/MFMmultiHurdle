rm(list = ls())
################################################################################
source("MFMmultiHurdle.R")
library(readxl)

library(lpSolve) #to load mcclust
library(mcclust.ext) #to estimate clustering


#which_analysis = "real"  #to reproduce results on real data
which_analysis = "synt" #to reproduce results on synthetic data

if (which_analysis == "real"){
  dat <- read_excel("data_clean.xlsx")
  n = 1154; d=7; Tmax=7
  y = array(, dim = c(n, d, Tmax))
  data = dat[1:49] 
  n = dim(data)[1]
  tab = summary(data)
  
  count = 1
  for(jj in (1:7)){
    y[, jj, ] = as.matrix(data[, count:(count+6)])
    count = count + 7
  }
  
  #compute n, d, Tmax
  n = dim(y)[1]
  d = dim(y)[2]
  Tmax = dim(y)[3]
  
  data = rowSums(y, dim=2)
  
  #compute data summaries 
  #bernoulli part
  ybin = y
  ybin[ybin > 1] = 1
  ybin_count = rowSums(ybin, dim=2)
  
  #shifted binomial part
  yminus1 = y - 1
  yminus1[yminus1 < 0] = 0
  yminus1 = rowSums(yminus1, dim=2)
  ypos_count = rowSums((y > 0),dim=2)
  
}else if(which_analysis == "synt"){ 
  #simulated data 
  set.seed(0)
  n = 1000; Tmax = 7; d = 7; dact = c(2,4,5,7); nact = setdiff(c(1:7), dact)
  y = array(, dim = c(n, d, Tmax))
  
  temp_cor_bin = sample(2, n, replace = TRUE, prob = c(3/7, 4/7))
  temp_cor_bin[temp_cor_bin==1] = sample(c(0,1), sum(temp_cor_bin==1), replace = TRUE)
  temp_cor_bin[temp_cor_bin==2] = sample(seq(2,9), sum(temp_cor_bin==2), replace = TRUE)
  
  y[temp_cor_bin==1,dact,1] = c(1,1,1,1)    
  y[temp_cor_bin==2,dact,1] = c(1,1,1,0)    
  y[temp_cor_bin==3,dact,1] = c(1,1,0,1)
  y[temp_cor_bin==4,dact,1] = c(1,0,1,1)
  y[temp_cor_bin==5,dact,1] = c(0,1,1,1)
  
  
  y[temp_cor_bin==0,dact,1] = c(0,0,0,0)
  y[temp_cor_bin==6,dact,1] = c(0,0,0,1)
  y[temp_cor_bin==7,dact,1] = c(0,0,1,0)
  y[temp_cor_bin==8,dact,1] = c(0,1,0,0)
  y[temp_cor_bin==9,dact,1] = c(1,0,0,0)
  
  cluster = rep(0,n)
  cluster[temp_cor_bin==1] = 1
  cluster[temp_cor_bin==2] = 1
  cluster[temp_cor_bin==3] = 1
  cluster[temp_cor_bin==4] = 1
  cluster[temp_cor_bin==5] = 1
  
  for (tt in 2:Tmax){
    temp1 = sample(5, sum(cluster), prob= c(3/7, rep(1/7,4)), replace = TRUE)
    temp0 = sample(c(0,6,7,8,9), n - sum(cluster), prob= c(3/7, rep(1/7,4)), replace = TRUE)
    
    y[cluster==1,dact,tt][temp1==1,] = c(1,1,1,1)    
    y[cluster==1,dact,tt][temp1==2,] = c(1,1,1,0)    
    y[cluster==1,dact,tt][temp1==3,] = c(1,1,0,1)
    y[cluster==1,dact,tt][temp1==4,] = c(1,0,1,1)
    y[cluster==1,dact,tt][temp1==5,] = c(0,1,1,1)
    
    y[cluster==0,dact,tt][temp0==0,] = c(0,0,0,0)
    y[cluster==0,dact,tt][temp0==6,] = c(0,0,0,1)
    y[cluster==0,dact,tt][temp0==7,] = c(0,0,1,0)
    y[cluster==0,dact,tt][temp0==8,] = c(0,1,0,0)
    y[cluster==0,dact,tt][temp0==9,] = c(1,0,0,0)
  }
  
  for (tt in 1:Tmax){
    for (jj in nact){
      y[,jj,tt] = rbinom(n, 1, 0.5)
    }
  }
  
  y[y>0] = rpois(sum(y), 3) + 1
  
  
  #compute data summaries 
  #bernoulli part
  ybin = y
  ybin[ybin > 1] = 1
  ybin_count = rowSums(ybin, dim=2)
  
  #shifted binomial part
  yminus1 = y - 1
  yminus1[yminus1 < 0] = 0
  yminus1 = rowSums(yminus1, dim=2)
  ypos_count = rowSums((y > 0),dim=2)
}
################################################################################
#c.hyperparam
a = 1 #1 #not used
b = 1 #2 #not used
#prior over Bernoulli parameter of the zero/non-zero part
alpha = 0.5 #hyper3 #beta shape1
beta = 0.5 #hyper4 #beta shape 2  
#prior over the shifted negative binom
zeta = 0.5 #hyper5 #geometric over 1,2,... success probability
eta = 1 #hyper6 #beta shape1 
lambda = 1 #7 #beta sahape2
mu = 1 #hyper8 Poisson prior over M
nu = 1 #hyper9 Poisson prior over S_m
sigma = 1 #hyper10 #gamma prior over un_norm jumps both out and in
hyper = c(a, b, alpha, beta, zeta, eta, lambda, mu, nu, sigma)

################################################################################
#d.initialization
M = 4 #number of outer components (obs + unbos)
S = rep(2, M) #number of inner components (obs + unbos)
w_un = rep(1, M) #unnormalized upper weights
q_un = rep(1, sum(S)) #unnormalized lower weights
p = matrix(rbeta(d * M, hyper[1], hyper[2]), nrow = d) #outer atoms
r = matrix(rgeom(d * sum(S), hyper[5]), nrow = d) + 1 #inner atoms 1
theta = matrix(rbeta(d * sum(S), hyper[6], hyper[7]), nrow = d) #inner atoms 2
#cluster allocations randomly inizialised
s = rep(0, n)
m = sample(M, n, replace=TRUE, prob = w_un / sum(w_un))
for(mm in unique(m)){
  nm = length(s[m==mm])
  if(mm == 1){a = 0; b = S[1]}else{ a = sum(S[1:(mm-1)]) ; b = a + S[mm]}
  q_temp = q_un[(a + 1):b]
  s[m==mm] = sample(seq((a + 1),b), nm, replace = TRUE, 
                    prob = q_temp / sum(q_temp))
}
out = rearrange_ms(m, s, S) #reorder components
m = out$m; s=out$s; stom = out$stom; S=out$S; k=out$k; h_m=out$h_m
################################################################################
#e.MCMC
burnin = 5000
totiter = 10000
x_comp = seq(1:10000) - 1
geodist = dgeom(seq(1,1000)-1, prob=hyper[5], log = TRUE)

for (iter in (1:(burnin + totiter))){
  
  #1. messages
  if(iter==1){
    print(paste("MCMC has started, a message is displayed every 100 iterations showing the iteration number and the outer clusters frequencies"))
    start_time <- Sys.time()}
  if(iter==2){
    end_time <- Sys.time()
    print(paste("the first iteration was completed in", round(end_time - start_time, 3), "sec" ))
    }
  if((iter/100) == floor(iter/100)){print(paste("iter:",c(iter))); print(c("freq:",table(m)/n))}
  
  #2.sampling outer and inner clustering allocation
  out = sample_ms(y, m, s, w_un, q_un, p, 
                  r, theta, M, S, Tmax, d)
  m = out$m; s = out$s
  
  #3.sampling non parametric mixing distributions
  out = rearrange_ms(m, s, S) #reorder components
  m = out$m; s=out$s; stom = out$stom; S=out$S; k=out$k; h_m=out$h_m
  
  u_up = rgamma(1, n, sum(w_un)) #auxiliary variable outer
  M = sampleM(k, u_up, x_comp, hyper) #sample M
  
  w_un = sample_w_un(M, m, u_up, hyper)
  p = sample_p(ybin_count, hyper, M, d, n, m)
  
  u_low = sample_u_low(M, m, q_un) #auxiliary variables inner
  S = sampleS(c(h_m, rep(0,(M-k))), u_low, x_comp, hyper) #sample S_1, ..., S_M
  out = rearrange_ms(m, s, S) #reorder components (m is already ordered, code can be optimized here)
  m = out$m; s=out$s; stom = out$stom; S=out$S; k=out$k; h_m=out$h_m
  
  q_un = sample_q_un(M, S, m, s, u_low, hyper)
  r = sampler(y, s, S, hyper)
  theta = sampletheta(yminus1, ypos_count, s, S, hyper, r, d)
  
  write.table(t(m), "upper_clustering.csv", append = TRUE, sep = ", ", row.names = FALSE, col.names = FALSE)
  write.table(t(s), "lower_clustering.csv", append = TRUE, sep = ", ", row.names = FALSE, col.names = FALSE)
  write.table(p, "bernoulli_param.csv", append = TRUE, sep = ", ", row.names = FALSE, col.names = FALSE)
  write.table(r, "negbin_param1.csv", append = TRUE, sep = ", ", row.names = FALSE, col.names = FALSE)
  write.table(theta, "negbin_param2.csv", append = TRUE, sep = ", ", row.names = FALSE, col.names = FALSE)
  write.table(M, "M.csv", append = TRUE, sep = ", ", row.names = FALSE, col.names = FALSE)
  write.table(t(S), "S.csv", append = TRUE, sep = ", ", row.names = FALSE, col.names = FALSE)
  write.table(t(w_un), "upper_weights.csv", append = TRUE, sep = ", ", row.names = FALSE, col.names = FALSE)
  write.table(t(q_un), "lower_weights.csv", append = TRUE, sep = ", ", row.names = FALSE, col.names = FALSE)
}
###############################################################################################################
#f. conditional MCMC
set.seed(0)
a = burnin + 1
b = burnin + totiter

m_chain = read.csv("upper_clustering.csv", header = FALSE)
m_fin = m_chain[ a : b , ]
psm_up = comp.psm( as.matrix( m_fin ) )
clusters_up = minbinder(psm_up, as.matrix( m_fin ), method = "all") 
m = clusters_up$cl[1,] 
M = 3 #To be set to MAP of previous algorithm

S = rep(2, M) #number of lower clusters (obs + unbos)
q_un = rep(1, sum(S)) #unnormalized lower weights
r = matrix(rgeom(d * sum(S), hyper[5]), nrow = d) + 1 #inner atoms 1
theta = matrix(rbeta(d * sum(S), hyper[6], hyper[7]), nrow = d) #inner atoms 2
#cluster allocations
s = rep(0, n)
for(mm in unique(m)){
  nm = length(s[m==mm])
  if(mm == 1){a = 0; b = S[1]}else{ a = sum(S[1:(mm-1)]) ; b = a + S[mm]}
  q_temp = q_un[(a + 1):b]
  s[m==mm] = sample(seq((a + 1),b), nm, replace = TRUE, 
                    prob = q_temp / sum(q_temp))
}
out = rearrange_ms(m, s, S) #reorder components
m = out$m; s=out$s; stom = out$stom; S=out$S; k=out$k; h_m=out$h_m

for (iter in (1:(burnin + totiter))){
  
  #1. messages
  if(iter==1){
    print(paste("conditional MCMC has started, a message is displayed every 100 iterations showing the iteration number and the inner clusters frequencies"))
    start_time <- Sys.time()}
  if(iter==2){
    end_time <- Sys.time()
    print(paste("the first iteration was completed in", round(end_time - start_time, 3), "sec" ))
  }
  if((iter/100) == floor(iter/100)){print(paste("iter:",c(iter))); print(c("freq:",table(s)/n))}
  #2.sampling outer and lower clustering allocation
  out = cond_sample_s(y, m, s, w_un, q_un, p, 
                  r, theta, M, S, Tmax, d)
  s = out$s
  
  #3.sampling non parametric mixing distributions
  out = rearrange_ms(m, s, S) #reorder components
  m = out$m; s=out$s; stom = out$stom; S=out$S; k=out$k; h_m=out$h_m

  p = sample_p(ybin_count, hyper, M, d, n, m)
  
  u_low = sample_u_low(M, m, q_un) #auxiliary variables inner
  S = sampleS(c(h_m, rep(0,(M-k))), u_low, x_comp, hyper) #sample S_1, ..., S_M
  out = rearrange_ms(m, s, S) #reorder components (m is already ordered, code can be optimized here)
  m = out$m; s=out$s; stom = out$stom; S=out$S; k=out$k; h_m=out$h_m
  
  q_un = sample_q_un(M, S, m, s, u_low, hyper)
  r = sampler(y, s, S, hyper)
  theta = sampletheta(yminus1, ypos_count, s, S, hyper, r, d)
  
  write.table(t(s), "cond_lower_clustering.csv", append = TRUE, sep = ", ", row.names = FALSE, col.names = FALSE)
  write.table(r, "cond_negbin_param1.csv", append = TRUE, sep = ", ", row.names = FALSE, col.names = FALSE)
  write.table(theta, "cond_negbin_param2.csv", append = TRUE, sep = ", ", row.names = FALSE, col.names = FALSE)
  write.table(t(S), "cond_S.csv", append = TRUE, sep = ", ", row.names = FALSE, col.names = FALSE)
  write.table(t(q_un), "cond_lower_weights.csv", append = TRUE, sep = ", ", row.names = FALSE, col.names = FALSE)
}
