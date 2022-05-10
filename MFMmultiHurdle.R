#Functions used for the paper 
#"Bayesian clustering of multiple zero-inflated outcomes"
#by Beatrice Franzolini, Andrea Cremaschi, Willem van den Boom, and Maria De Iorio
#Enriched MFM multiHurdle

MBern <-function(y, hyper, log=TRUE){
  n0 = sum(y==0)
  n1 = prod(dim(y)) - n0
  out = lbeta(hyper[3] + n1, hyper[4] + n0) - lbeta(hyper[3], hyper[4])
  if(!log){ out = exp(out) }
  return (out)
}

Mtheta <- function(y, hyper, geodist, log = TRUE){
  x = y - 1; n = sum(y>0); x = x[x>=0]
  out = lbeta(hyper[6] + sum(x) , hyper[7] + n * c(1:length(geodist))) - lbeta(hyper[6], hyper[7])
  for (xx in x){
    temp = lgamma(xx + c(1:length(geodist))) - lgamma(c(1:length(geodist))) - lgamma(xx + 1)
    out = out + temp
  }
  out = sum(exp(out + geodist))
  if(!log){out = log(out)}
  return(out)
}


f_0arrow1 <-function(jj, y, m, s, hyper, geodist){
  out = 0
  for (mm in unique(m)){
    y_temp = y[(m==mm),jj, ,drop = FALSE]
    out = out + MBern(y_temp, hyper, log=TRUE)
    for (ss in unique(s[m==mm])){
      y_temp = y[(m==mm & s==ss),jj,] 
      out = out + Mtheta(y_temp, hyper, geodist, log=TRUE)
    }}
  y_temp = y[,jj, ,drop = FALSE]
  out = out - MBern(y_temp, hyper, log=TRUE) - 
    Mtheta(y_temp, hyper, geodist, log=TRUE)
  return(out)
}


f_1arrow0 <-function(jj, y, m, s, hyper, geodist){
  out = 0
  for (mm in unique(m)){
    y_temp = y[(m==mm),jj, ,drop = FALSE]
    out = out - MBern(y_temp, hyper, log=TRUE)
    for (ss in unique(s[m==mm])){
      y_temp = y[(m==mm & s==ss),jj,] 
      out = out - Mtheta(y_temp, hyper, geodist, log=TRUE)
    }}
  y_temp = y[,jj, ,drop = FALSE]
  out = out + MBern(y_temp, hyper, log=TRUE) +
    Mtheta(y_temp, hyper, geodist, log=TRUE)
  return(out)
}

#k for 1 obs
kernelk <-function(y, p, r, theta, Tmax, d, log = TRUE){
  y_temp = y
  y_temp[y_temp>0] = 1
  d = dim(y)[2]
  n1 = rowSums(y_temp, dim = 2); n0 = Tmax - n1
  out = sum(n0 * log(1 - p) + n1 * log(p))
  x = y - 1
  for (jj in 1:d){
    temp = x[,jj,]
    if (length(temp[temp>=0])>0){
      out = out + sum(dnbinom(temp[temp>=0], r[jj], theta[jj], log=TRUE))
      
    }
  }
  if(!log){ out = exp(out) }
  return (out)
} 

sample_ms <-function(y, m, s, w_un,q_un, p, r, theta, M, S, Tmax, d){
  for (i in 1:dim(y)[1]){
    prob = NULL
    for (ss in 1:sum(S)){
      mm = stom[ss] 
      prob = c(prob, log(w_un[mm]) + log(q_un[ss]) + kernelk(y[i,,,drop = FALSE], p[, mm], r[, ss],
                                                             theta[, ss], Tmax, d, log=TRUE) )
    }
    if(sum(exp(prob)) == 0){prob = prob - min(prob)}
    if(sum(exp(prob)) == Inf){prob = prob - max(prob)}
    s[i] = sample(sum(S), 1, prob = exp(prob) / sum(exp(prob)))
    m[i] = stom[s[i]]
  }
  return(list(m=m, s=s))
}

cond_sample_s <-function(y, m, s, w_un,q_un, p, r, theta, M, S, Tmax, d){
  for (i in 1:dim(y)[1]){
    prob = NULL
    for (ss in 1:sum(S)){
      mm = stom[ss]
      if(m[i] == mm){
      prob = c(prob, log(w_un[mm]) + log(q_un[ss]) + kernelk(y[i,,,drop = FALSE], p[, mm], r[, ss],
                                                             theta[, ss], Tmax, d, log=TRUE) )
      }else{
        prob = c(prob, NA)
      }
    }
    if(sum(exp(prob), na.rm = TRUE) == 0){prob = prob - min(prob, na.rm = TRUE)}
    if(sum(exp(prob), na.rm = TRUE) == Inf){prob = prob - max(prob, na.rm = TRUE)}
    prob[is.na(prob)] = -Inf
    s[i] = sample(sum(S), 1, prob = exp(prob) / sum(exp(prob)))
    #m[i] = stom[s[i]]
  }
  return(list(s=s))
}


rearrange_ms <- function(m, s, S){
  m_temp = m; S_temp = S; s_temp = s; stom = NULL
  
  m_unique = unique(m); k = length(m_unique); h_m = rep(0, k)
  old_m = m_unique[order(m_unique)][order(table(m), decreasing = TRUE)]
  
  count_m = 1; count_s = 1
  
  for (x in old_m){
    m_temp[m==x] = count_m
    S_temp[count_m] = S[x]
    s_unique = unique(s[m==x]); h_m[count_m] = length(s_unique)
    old_s = s_unique[order(s_unique)][order(table(s[m==x]), decreasing = TRUE)]
    for (y in old_s){
      s_temp[s==y] = count_s
      stom[count_s] = count_m
      count_s = count_s + 1
    }
    if(S_temp[count_m] > h_m[count_m]){
      for (y in seq(1, (S_temp[count_m] - h_m[count_m]))){
        stom[count_s] = count_m 
        count_s = count_s + 1
      }}
    count_m = count_m + 1
  }
  
  for (x in setdiff(seq(1,M),m_unique)){
    S_temp[count_m] = S[x]
    for (y in (1:S_temp[count_m])){
      stom[count_s] = count_m
      count_s = count_s + 1
    }
    count_m = count_m + 1
  }
  return(list(m=m_temp, s=s_temp, stom=stom, S=S_temp, k=k, h_m=h_m))
}

phi_dir<-function(u, hyper){
  return ((u+1)^(-hyper[10]))
}

sample_u_low <-function(M, m, q_un){
  n_m = c(as.vector(table(m)[unique(m)]))
  un_m = unique(m)
  k = length(un_m)
  sumq = rep(0,k)
  for (mm in un_m){
    if(mm == 1){
      a = 0; b = S[1]
    }else{
      a = sum(S[1:(mm-1)]) ; b = a + S[mm]
    }
    sumq[mm] = sum(q_un[(a+1):b])
  }
  u_low = c(rgamma(k, n_m, sumq)) #auxiliary variables down
  return(u_low)
}

sampleM <-function(k, u, x, hyper){
  prob = log(x  + k) - log(factorial(x)) + x * log( phi_dir(u, hyper) * hyper[8])
  prob = exp(prob) 
  prob = prob / sum(prob)
  x = sample(x, 1, prob = prob)
  return (k + x)
}

sampleS<-function(h_m, u_low, x, hyper){
  M = length(h_m)
  new = rep(0, M)
  for (mm in 1:M){
    if (h_m[mm]>0){
      prob = log(x  + h_m[mm]) - log(factorial(x)) + x * log( phi_dir(u_low[mm], hyper) * hyper[9])
      prob = exp(prob) 
      prob = prob / sum(prob)
      new[mm] = sample(x, 1, prob = prob)
    }else{
      new[mm] = rpois(1, hyper[9]) + 1
    }
  }
  return (h_m + new)
}

sample_w_un<-function(M, m, u, hyper){
  empty = M - length( unique(m)) 
  n = c(as.vector(table(m)), rep(0,empty))
  w_un = rgamma(length(n), n + hyper[10], 1 + u)
  return(w_un)
}

sample_q_un<-function(M, S, m, s, u_low, hyper){
  q_un = rep(0, sum(S))
  for (mm in (1:M)){
    if(mm == 1){
      a = 0; b = S[1]
    }else{
      a = sum(S[(1:(mm-1))]); b = a + S[mm]
    }
    if(length(s[m == mm])>0){
      q_un[(a+1):b] = sample_w_un(S[mm], s[m==mm], u_low[mm], hyper)
    }else{ 
      q_un[(a+1):b] = rgamma(S[mm], hyper[10], 1)
    }
  }
  return(q_un)
}

sample_p<-function(ybin_count, hyper, M, d, n, m){
  p = matrix(rep(0, d*M), nrow = d,ncol = M)
  for (mm in 1:M){
    n1 = colSums(ybin_count[m==mm,,drop = FALSE])
    n2 = Tmax * sum(m==mm) - n1
    p[,mm] = rbeta(d, hyper[3] + n1, hyper[4] + n2)
  }
  return(p)
}

sampler <- function(y, s, S, hyper){
  r = matrix(, nrow = d, ncol = sum(S))
  for (jj in seq(1,d)){
    for (ss in unique(s)){
      x = y[s==ss,jj,] - 1; n = length(x[x>=0]); x = x[x>=0]
      out = lbeta(hyper[6] + sum(x) , hyper[7] + n * seq(1:length(geodist)) ) - 
        lbeta(hyper[6], hyper[7])
      for (xx in x[x>0]){
        temp = lgamma(xx + seq(1:length(geodist))) -
          lgamma(seq(1:length(geodist))) - 
          lgamma(xx + 1)
        out = out + temp
      }
      out = out + geodist
      if(sum(exp(out)) == 0){out = out - min(out)}
      if(sum(exp(out)) == Inf){out = out - max(out)}
      prob = exp(out)
      r[jj, ss] = sample(length(geodist), 1, prob = prob/sum(prob))
    }
    for (ss in setdiff(seq(1,sum(S)),unique(s))){
      r[jj, ss] = rgeom(1, hyper[5]) + 1
    }
  }
  return(r)
}

sampletheta<-function(yminus1, ypos_count, s, S, hyper, r, d){
  theta = matrix(, nrow = d, ncol = sum(S))
  for (ss in unique(s)){
    eta_new = hyper[6] + colSums(ypos_count[s==ss, , drop = FALSE]) * r
    lambda_new = hyper[7] + colSums(yminus1[s==ss, , drop = FALSE])     
    theta[, ss] = rbeta(d, eta_new , lambda_new)
  }
  for (ss in setdiff(seq(1,sum(S)),unique(s))){
    theta[, ss] = rbeta(d, hyper[6], hyper[7])
  }
  return(theta)
}