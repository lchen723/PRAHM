library(MASS)
library(pracma)
library (rootSolve)
library(MultiRNG)
library(ahaz)
corrected_ridge_AH_with_bootstrap <- function(time, status, W, Z, Sigma_e, rep.measurement.times,
                                              max_outer_iter = 10, tol = 1e-6,
                                              optim_ctrl = list(maxit = 200, reltol = 1e-8),
                                              verbose = TRUE, bootstrap_times, bootstrap_property,
                                              critical.value = NULL, trim = 0) {
  
  n0 = length(time)
  beta_hat_M = matrix(NA,bootstrap_times,ncol(W)+ncol(Z))
  omega_M = matrix(NA,bootstrap_times,ncol(W)+ncol(Z))
  index_M = matrix(NA,ceiling(n0*bootstrap_property),bootstrap_times)
  ni = rep.measurement.times
  kappa = Inf
  
  status = status[order(time)]
  W = W[order(time),]
  Z = Z[order(time),]
  time = time[order(time)]
  
  status0 = status
  W0 = W
  Z0 = Z
  time0 = time
  tau0 = max(time)
  event_times0 <- sort(unique(time[status == 1]))
  m0 <- length(event_times0)
  if (m0 == 0) stop("No events found (status==1).")
  R_mat0 <- outer(time0, event_times0, FUN = function(y, t) as.numeric(y >= t))
  dN_mat0 <- outer(time0, event_times0, FUN = function(y, t) as.numeric(y == t)) * status
  for (b in seq(bootstrap_times)) {
    boot_index = sample(1:n0,ceiling(n0*bootstrap_property),replace = T)
    boot_index = boot_index[order(boot_index)]
    index_M[,b] = boot_index
    status = status0[boot_index]
    W = W0[boot_index,]
    Z = Z0[boot_index,]
    time = time0[boot_index]+runif(length(boot_index),0,0.000001)
    
    n = length(time)
    tau = max(time)
    W = as.matrix(W)
    Z = as.matrix(Z)
    V_star = cbind(W, Z)
    p_x = ncol(W)
    p_z = ncol(Z)
    p = p_x + p_z
    if (is.null(kappa)) kappa = 1/p
    
    # Sigma_e0
    Sigma_e0 <- as.matrix(Matrix::bdiag(Sigma_e, matrix(0, p_z, p_z)))
    
    # event times and indicator mats
    event_times <- sort(unique(time[status == 1]))
    m <- length(event_times)
    if (m == 0) stop("No events found (status==1).")
    R_mat <- outer(time, event_times, FUN = function(y, t) as.numeric(y >= t))
    dN_mat <- outer(time, event_times, FUN = function(y, t) as.numeric(y == t)) * status
    
    # --- USE ahaz to get A_base ---
    surv_obj <- survival::Surv(time, status)
    # call ahaz: second argument is covariate matrix (without intercept)
    fit_surv <- ahaz::ahaz(surv_obj, V_star)
    A_base <- fit_surv$D
    
    int2 = function(ts) {(1-1/sum(Y>ts))*sum((Y>ts)/ni)}
    corr_total = quad(int2,0,(tau-tol))*Sigma_e0
    A <- A_base - corr_total
    B <- fit_surv$d
    
    # --- helper: £N0 at each subject's Y_i ---
    Lambda0_subjects <- function(beta_hat) {
      dLambda <- numeric(m)
      for (j in seq_len(m)) {
        R_j <- R_mat[, j]; dN_j <- dN_mat[, j]; denom <- sum(R_j)
        if (denom == 0) { dLambda[j] <- 0; next }
        num <- sum(dN_j) - sum(R_j * as.numeric(V_star %*% beta_hat))
        dLambda[j] <- num / denom
      }
      Lambda0_vec <- numeric(n)
      for (i in seq_len(n)) {
        idx <- which(event_times <= time[i])
        if (length(idx) == 0) {
          Lambda0_vec[i] <- 0
        } else {
          Lambda0_vec[i] <- sum(dLambda[idx])
        }
      }
      Lambda0_vec
    }
    
    # helper: compute D given omega
    compute_D_given_omega <- function(omega_vec) {
      omega_vec <- pmax(omega_vec, 1e-6)
      Omega_inv <- diag(1 / omega_vec, p, p)
      S_mat <- A + Omega_inv + diag(1e-8, p)
      beta <- tryCatch(as.numeric(solve(S_mat, B)),
                       error = function(e) as.numeric(MASS::ginv(S_mat) %*% B))
      
      Lambda_subj <- Lambda0_subjects(beta)
      
      # residuals
      quad_e <- as.numeric(t(beta) %*% Sigma_e0 %*% beta)
      r_i <- numeric(n)
      for (i in seq_len(n)) {
        r_i[i] <- Lambda_subj[i] + time[i] * as.numeric(V_star[i, ] %*% beta) + (time[i]^2) * quad_e
      }
      rM <- status - r_i
      inner_log_arg <- pmax(1e-8, status - rM)
      tmp <- -2 * (rM + status * log(inner_log_arg))
      tmp <- pmax(tmp, 0)
      rD <- sign(rM) * sqrt(tmp)
      D_val <- sum(rD^2)
      
      list(D = D_val, beta = beta, Lambda_subjects = Lambda_subj, rD = rD, rM = rM)
    }
    
    # outer loop
    omega <- rep(0.01, p)
    omega <- pmax(omega, 1e-6)
    prev_D <- Inf
    beta_hat <- NULL
    Lambda0_subj <- NULL
    
    for (outer in seq_len(max_outer_iter)) {
      if (verbose) message(sprintf("Outer iter %d: sum(omega)=%.6g", outer, sum(omega)))
      obj_for_optim <- function(omega_try) {
        pen <- 0
        ssum <- sum(omega_try)
        if (ssum > kappa) pen <- 1e6 * (ssum - kappa)^2
        res <- compute_D_given_omega(omega_try)
        res$D + pen
      }
      lower <- rep(1e-6, p); upper <- rep(kappa, p)
      start_omega <- pmax(omega, 0.01)
      opt <- tryCatch({
        optim(par = start_omega, fn = obj_for_optim, method = "L-BFGS-B",
              lower = lower, upper = upper, control = optim_ctrl)
      }, error = function(e) {
        if (verbose) message("optim failed; keep current omega")
        list(par = start_omega, value = obj_for_optim(start_omega), convergence = 1)
      })
      info_new <- compute_D_given_omega(opt$par)
      D_new <- info_new$D
      beta_new <- info_new$beta
      Lambda_subj_new <- info_new$Lambda_subjects
      if (verbose) message(sprintf("  After optim: D=%.6g, sum(omega)=%.6g", D_new, sum(opt$par)))
      if (abs(prev_D - D_new) < tol) {
        omega <- opt$par
        beta_hat <- beta_new
        Lambda0_subj <- Lambda_subj_new
        prev_D <- D_new
        if (verbose) message("Converged by D change.")
        break
      }
      omega <- opt$par
      beta_hat <- beta_new
      Lambda0 <- Lambda_subj_new
      prev_D <- D_new
    }
    
    beta_hat_M[b,] = beta_hat
    omega_M[b,] = omega
  }
  mean2 = function(x){
    x = sort(x)
    n0 = length(x)
    l = round(n0*trim)+1
    u = n0-l+1
    mean(x[l:u])
  }
  sd2 = function(x){
    x = sort(x)
    n0 = length(x)
    l = round(n0*trim)+1
    u = n0-l+1
    sd(x[l:u])
  }
  beta = apply(beta_hat_M,2,mean2)
  sd = apply(beta_hat_M,2,sd2)
  t.value = beta/sd
  p.value = 2*(1-pt(abs(t.value),n-p-1))
  
  ##refit part
  if(is.null(critical.value)){
    sel = which(p.value<0.05) 
  }else{
    sel = which((abs(beta))>critical.value) 
  }
  if(length(sel)==0){
    re.beta_M = matrix(0,bootstrap_times,p)
  }else{
    V_star02 = as.matrix(cbind(W0,Z0)[,sel])
    Sigma_e02 = Sigma_e0[sel,sel]
    re.beta_M = matrix(NA,bootstrap_times,p)
    for (b in seq(bootstrap_times)) {
      index2 = index_M[,b]
      status = status0[index2]
      V_star2 = V_star02[index2,]
      time = time0[index2]+runif(length(index2),0,0.000001)
      surv_obj2 = survival::Surv(time, status)
      tau = max(time)
      # call ahaz: second argument is covariate matrix (without intercept)
      fit_surv = ahaz::ahaz(surv_obj2, V_star2)
      A_base2 = fit_surv$D
      
      int2 = function(ts) {(1-1/sum(Y>ts))*sum((Y>ts)/ni)}
      corr_total = quad(int2,0,(tau-tol))*Sigma_e02
      A2 = A_base2 - corr_total
      B2 <- fit_surv$d
      
      re.beta = rep(0,p)
      re.beta[sel] = solve(A2)%*%B2
      re.beta_M[b,] = re.beta
    }
  }
  re.beta = apply(re.beta_M,2,mean2)
  #re.sd = apply(re.beta_M,2,sd)
  re.t.value = re.beta/sd
  re.p.value = 2*(1-pt(abs(re.t.value),n-length(sel)-1))
  
  Lambda0_subjects_2 <- function(beta_hat) {
    dLambda <- numeric(m0)
    for (j in seq_len(m0)) {
      R_j = R_mat0[, j]
      dN_j = dN_mat0[, j]
      denom = sum(R_j)
      if (denom == 0) { dLambda[j] <- 0; next }
      num <- sum(dN_j) - sum(R_j * as.numeric(cbind(W0,Z0) %*% beta_hat))
      dLambda[j] <- num / denom
    }
    Lambda0_vec <- numeric(n0)
    for (i in seq_len(n0)) {
      idx <- which(event_times0 <= time0[i])
      if (length(idx) == 0) {
        Lambda0_vec[i] <- 0
      } else {
        Lambda0_vec[i] <- sum(dLambda[idx])
      }
    }
    Lambda0_vec
  }
  
  return(list(
    coef = cbind(beta,sd,t.value,p.value),
    refit.coef = cbind(re.beta,sd,re.t.value,re.p.value),
    Lambda0 = cbind(time0[order(time0)],Lambda0_subjects_2(beta)),
    re.Lambda0 = cbind(time0[order(time0)],Lambda0_subjects_2(re.beta))
  ))
}
#########################
n = 200 #sample size
ni = 2  #replicate size

beta_X = c(1,1,0,0)  
px = length(beta_X) #dimension of X
ax=matrix(1:px,px,px)
tempSigmax=(0.5^abs(ax-t(ax))) 

beta_Z = c(1,0,0,0) 
pz = length(beta_Z) #dimension of Z
az=matrix(1:pz,pz,pz)
tempSigmaz=(0.5^abs(az-t(az))) 
p = px+pz

#measurement error setting
SE = 0.25
ME_CM = diag(px)*SE
ME_CM_0 = diag(c(diag(ME_CM),rep(0,pz)))

beta1 = c()
beta2 = c()
beta3 = c()
beta4 = c()
nai_beta = c()
nai_beta_ME = c()
nai_beta_ME_adj = c()
nai_beta_ME_adj = c()
for (i in seq(100)) {
  if(px==1){
    X = runif(300,0,1)
  }else{X = draw.d.variate.uniform(300,px,tempSigmax)}
  #X = draw.d.variate.uniform(n,px,tempSigmax) #real value
  W1 = X+mvrnorm(300,rep(0,px),ME_CM)           #replicate 1
  W2 = X+mvrnorm(300,rep(0,px),ME_CM)           #replicate 2
  if(pz==1){
    Z = runif(300,0,1)
  }else{Z = cbind(rbinom(300,1,0.5),draw.d.variate.uniform(300,pz-1,tempSigmaz[-1,-1]))}
  #Z = draw.d.variate.uniform(n,pz,tempSigmaz) 
  W = (W1+W2)/ni                               #W_bar
  
  V = cbind(X,Z)
  V_s = cbind(W,Z)
  beta = c(beta_X,beta_Z)
  
  U = runif(300,0,1)
  Ts = -log(1-U)/(1+V%*%beta)
  C = rexp(300,1)
  sigma_0_hat = (t(W1-W)%*%(W1-W) + t(W2-W)%*%(W2-W))/300
  sigma_1_hat = diag(c(diag(sigma_0_hat),rep(0,pz)))
  
  T1 = cbind(Ts,C)
  
  Y = apply(T1,1,min)
  delta = (Ts<C)*1
  
  rei = which(((V%*%beta+Ts>0)*(V_s%*%beta+Ts>0)) == 1)
  Y = Y[rei]
  delta = delta[rei]
  X = X[rei,]
  W = W[rei,]
  Z = Z[rei,]
  V = V[rei,]
  V_s = V_s[rei,]
  
  Y = Y[1:n]
  delta = delta[1:n]
  X = X[1:n,]
  W = W[1:n,]
  W1 = W1[1:n,]
  W2 = W2[1:n,]
  Z = Z[1:n,]
  V = V[1:n,]
  V_s = V_s[1:n,]
  
  
  sigma_0_hat = (t(W1-W)%*%(W1-W) + t(W2-W)%*%(W2-W))/ n
  sigma_1_hat = diag(c(diag(sigma_0_hat),rep(0,pz)))
  
  As2 = corrected_ridge_AH_with_bootstrap(Y,delta,X,Z,ME_CM*0,2,bootstrap_times = 100, bootstrap_property = 1, critical.value = 0.1, trim = 0.16)
  As2_ME = corrected_ridge_AH_with_bootstrap(Y,delta,W,Z,sigma_0_hat,2,,bootstrap_times = 100, bootstrap_property = 1, critical.value = 0.1, trim = 0.16)
  
  surv = Surv(Y,delta)
  fit_surv = ahaz(surv,V)
  fit_surv_ME = ahaz(surv,V_s)
  int2 = function(ts) {(1-1/sum(Y>ts))*sum((Y>ts)/2)}
  corr_total = quad(int2,0,(max(Y)-0.000001))*ME_CM_0
  
  A = fit_surv$D
  B = fit_surv$d
  beta_n = solve(A, B)
  A_ME = fit_surv_ME$D
  B_ME = fit_surv_ME$d
  A_ME2 = A_ME - corr_total
  beta_n = solve(A, B)
  beta_n_ME = solve(A_ME, B_ME)
  beta_n_ME_adj = solve(A_ME2, B_ME)
  
  beta1[[i]] = As2$coef
  beta2[[i]] = As2$refit.coef
  beta3[[i]] = As2_ME$coef
  beta4[[i]] = As2_ME$refit.coef
  nai_beta = cbind(nai_beta,beta_n)
  nai_beta_ME = cbind(nai_beta_ME,beta_n_ME)
  nai_beta_ME_adj = cbind(nai_beta_ME_adj,beta_n_ME_adj)
}
apply(beta1, 1, mean)
apply(beta2, 1, mean)
apply(beta3, 1, mean)
apply(beta4, 1, mean)
matrix_mean = function(M,p,dim){
  a = rep(0,p)
  for (i in seq(length(M))) {
    a = a+M[[i]][,dim]
  }
  a/length(M)
}
sel_bias = function(M){
  a = rep(0,8)
  b = rep(0,8)
  for (i in seq(length(M))) {
    a = a+(M[[i]][,1]>0.1)
    b = b+abs(M[[i]][,1]-c(1,1,0,0,1,0,0,0))
  }
  return(cbind(a,b/length(M)))
}
cbind(matrix_mean(beta4,4,1),matrix_mean(beta4,4,2),matrix_mean(beta4,4,1)/matrix_mean(beta4,4,2),2*(1-pt(abs(matrix_mean(beta4,4,1)/matrix_mean(beta4,4,2)),200-8-1)),sel_bias(beta4))
cbind(apply(nai_beta_ME, 1, mean),apply(nai_beta_ME, 1, sd),apply(nai_beta_ME, 1, mean)/apply(nai_beta_ME, 1, sd),2*(1-pt(abs(apply(nai_beta_ME, 1, mean)/apply(nai_beta_ME, 1, sd)),200-8-1)),apply(abs(nai_beta_ME)>0.1, 1, sum),apply(abs(nai_beta_ME-c(1,1,0,0,1,0,0,0)), 1, mean))