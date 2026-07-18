library(medicaldata)
ridgereg_ME = function(w,y,gam,penalty.factor,ME_CM) {
  H = w%*%diag(sqrt(gam))%*%solve(diag(sqrt(gam))%*%(t(w)%*%w-nrow(w)*ME_CM)%*%diag(sqrt(gam))+diag(penalty.factor))%*%diag(sqrt(gam))%*%t(w)
  coef = diag(sqrt(gam))%*%solve(diag(sqrt(gam))%*%(var(w)-ME_CM)%*%diag(sqrt(gam))+diag(penalty.factor))%*%diag(sqrt(gam))%*%cov(w,y)
  D = diag(sqrt(gam))%*%solve(diag(sqrt(gam))%*%(var(w)-ME_CM)%*%diag(sqrt(gam))+diag(penalty.factor))%*%diag(sqrt(gam))%*%t(w)
  return(list(ssr = sum((y-H%*%y)^2), df = sum(diag(H)), coef = coef, D = D))
}
ridge_operator<-function(w, y, kappa, gaminit=NULL, ME_CM, adaptivel = F, critical.value = 0) {
  
  n = dim(w)[1]
  p = dim(w)[2]
  
  if(adaptivel) {
    penalty.factor = 1/abs(coef(lm(y~w-1)))
  } else {
    penalty.factor = rep(1,p) #adaWeight*0+1
  }
  
  if(is.null(p)){p = 1}
  
  if(is.null(gaminit)) {
    gamcur = rep(kappa/p,p)
  } else {
    gamcur = kappa*gaminit/sum(gaminit)
  }
  
  kkk = 20
  mygrid = (0:kkk)/kkk
  test = 1
  
  count = 1
  
  while(test) {
    
    oldgam = gamcur
    
    for (j in 1:p) {
      gamcandj = rep(0,p)
      gamcandj[j] = 1
      
      gamcandMj = gamcur
      gamcandMj[j] = 0
      
      if (sum(gamcandMj)>0.01*kappa) {
        gamcandMj = gamcandMj/sum(gamcandMj)
        
        fffun = function(ttt){
          return(ridgereg_ME(w,y,(gamcandj*ttt+gamcandMj*(1-ttt))*kappa,penalty.factor,ME_CM)$ssr)
        }
        tttmin = optimize(fffun,c(0,1),tol = 0.0001)
        ttt = tttmin$minimum
        gamcur = (gamcandj*ttt+gamcandMj*(1-ttt))*kappa
      } # end if
    }  # end for over j
    count = count+1
    if(max(abs(gamcur-oldgam))<0.001*kappa) {test=0}
    if (count>20) {
      test = 0
      print("Takes more than 20 loops to converge in modified coordinate descent!!!!")
    }
  } # end while over test
  gamcur[which(gamcur<=2*0.001*kappa)]=0
  lamhat = gamcur
  ridge1 = ridgereg_ME(w,y,lamhat,penalty.factor,ME_CM)
  MSE_ridge = ridge1$ssr/(n-p-1)
  var_ridge = (MSE_ridge*ridge1$D%*%t(ridge1$D))/(n^2)
  Estimate =  ridge1$coef
  Std.Error = sqrt(diag(var_ridge))
  t.value = Estimate/Std.Error
  p.value = 2*(1-pt(abs(t.value),n-p-1))
  ridge_coef = cbind(Estimate, Std.Error, t.value, p.value)
  ##refit part
  sel = which((abs(ridge1$coef))>critical.value)
  w2 = w[ , sel]
  ME_CM2 = ME_CM[sel, sel]
  lamhat2 = lamhat[sel]
  re.Estimate = solve(t(w2)%*%w2-n*ME_CM2)%*%t(w2)%*%y
  re.residual = y-w2%*%re.Estimate
  re.MSE = as.numeric((t(y-w2%*%re.Estimate)%*%(y-w2%*%re.Estimate))/(n-length(sel)-1))
  re.Std.Error = sqrt(diag(re.MSE*solve(t(w2)%*%w2-n*ME_CM2+diag(1/lamhat2))))
  re.Estimate2 = rep(0,p)
  re.Std.Error2 = rep(0,p)
  re.Estimate2[sel] = re.Estimate
  re.Std.Error2[sel] = re.Std.Error
  re.t.value = re.Estimate2/re.Std.Error2
  re.p.value = 2*(1-pt(abs(re.t.value),n-length(sel)-1))
  refit_coef = cbind(re.Estimate2, re.Std.Error2, re.t.value, re.p.value)
  ##
  return(list(lamhat = lamhat, ridge = ridge1, ridge.coef = ridge_coef, selection = sel, refit.coef = refit_coef))
}
ridge_operator_DK<-function(w, y, gaminit=NULL, ME_CM, adaptivel = F, critical.value = 0){
  d = dim(w)
  KD = 1/d[2]
  KS = seq(KD/20,KD,length=20)
  SSE = rep(NA,20)
  for (i in seq(20)) {
    a = ridge_operator(w,y,KS[i], gaminit, ME_CM, adaptivel, critical.value)
    SSE[i] = sum((y-w%*%a$refit.coef[,1])^2)
  }
  kappa = KS[min(which(SSE == min(SSE)))]
  a = ridge_operator(w,y,kappa, gaminit, ME_CM, adaptivel, critical.value)
  return(list(lamhat = a$lamhat, ridge = a$ridge, ridge.coef = a$ridge.coef, selection = a$selection, refit.coef = a$refit.coef, kappa = kappa))
}
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
    
    # --- helper: Λ0 at each subject's Y_i ---
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
#################
library(BART)
data("ACTG175")
naive_method = c()
rep_naive_method = c()
Yi_method = c()
our_method = c()

n = 1000
for (i in 1:100) {
  boot = sample(1:2139,n,T)
  ACTG = ACTG175[boot,]
  #the data have the tied survival times, so I add a small random number to break it
  Y = as.numeric(ACTG[,'days'])+runif(n,0,0.00000001)
  delta = as.numeric(ACTG[,'cens'])[order(Y)]
  ##one hard encoding
  #strat
  strat1 = (ACTG[,'strat']==1)*1
  strat2 = (ACTG[,'strat']==2)*1
  #arms
  arms1 = (ACTG[,'arms']==1)*1
  arms2 = (ACTG[,'arms']==2)*1
  arms3 = (ACTG[,'arms']==3)*1
  
  #Adjusting the variables required for X and Z selection. Which variables use in X, need to put in W1 and W2 at the same time.
  #Z is the category variables and X is the continuous variables.
  Z = cbind(ACTG[,'treat'],ACTG[,'oprior'],ACTG[,'z30'],ACTG[,'race'],ACTG[,'gender'],ACTG[,'hemo'],ACTG[,'homo'],
            ACTG[,'drugs'],ACTG[,'str2'],ACTG[,'symptom'],ACTG[,'offtrt'],ACTG[,'r'])
  Z = Z[order(Y),]
  X = cbind(log(ACTG[,'cd40']+1),log(ACTG[,'cd80']+1),log(ACTG[,'preanti']+1),ACTG[,'wtkg'],ACTG[,'age'],ACTG[,'karnof'])
  X = X[order(Y),]
  Y = Y[order(Y)]
  W1 = cbind(log(ACTG[,'cd40']+1),log(ACTG[,'cd80']+1),log(ACTG[,'preanti']+1),ACTG[,'wtkg'],ACTG[,'age'],ACTG[,'karnof'])
  W2 = cbind(log(ACTG[,'cd420']+1),log(ACTG[,'cd820']+1),log(ACTG[,'preanti']+1),ACTG[,'wtkg'],ACTG[,'age'],ACTG[,'karnof'])
  W = (W1+W2)/2
  sigma_0_hat = (t(W1-W)%*%(W1-W) + t(W2-W)%*%(W2-W))/ n
  V = cbind(X,Z)
  V_s = cbind(W,Z)
  
  sigma_0_hat = (t(W1-W)%*%(W1-W) + t(W2-W)%*%(W2-W))/ n
  sigma_1_hat = diag(c(diag(sigma_0_hat),rep(0,12)))
  data = as.data.frame(cbind(Y,delta,X,Z,W1,W2))
  
  #just using cd40 and cd80
  A = A1_fun(Y,delta,V)
  B = B_fun(Y,delta,V)
  naive_method = cbind(naive_method,ginv(A)%*%B)
  
  #replicate measurements
  A_s = A1_fun(Y,delta,V_s)
  B_s = B_fun(Y,delta,V_s)
  rep_naive_method = cbind(rep_naive_method,ginv(A_s)%*%B_s)
  
  #correction with method by Yan & Yi
  A1 = A1_fun(Y,delta,V_s)
  int2 = function(ts) {(1-1/sum(Y>ts))*sum((Y>ts)/2)}
  A2t = quad(int2,0,tau-0.00001)*sigma_1_hat
  A_ss = A1-A2t
  B_s = B_fun(Y,delta,V_s)
  Yi_method = cbind(Yi_method,ginv(A_ss)%*%B_s)
  
  #our method
  ALL_ACTG_model_ME = corrected_ridge_AH_with_bootstrap(Y,delta,W,Z,sigma_0_hat,2,bootstrap_times = 2, bootstrap_property = 1, critical.value = 0.1, trim = 0.16)
  our_method = cbind(our_method,ALL_ACTG_model_ME$coef) 
}

cbind(apply(naive_method, 1, mean),apply(naive_method, 1, sd),
      apply(naive_method, 1, mean)/apply(naive_method, 1, sd),
      pnorm(abs(apply(naive_method, 1, mean)/apply(naive_method, 1, sd)),lower.tail = F))
cbind(apply(rep_naive_method, 1, mean),apply(rep_naive_method, 1, sd),
      apply(rep_naive_method, 1, mean)/apply(rep_naive_method, 1, sd),
      pnorm(abs(apply(rep_naive_method, 1, mean)/apply(rep_naive_method, 1, sd)),lower.tail = F))
cbind(apply(Yi_method, 1, mean),apply(Yi_method, 1, sd),
      apply(Yi_method, 1, mean)/apply(Yi_method, 1, sd),
      pnorm(abs(apply(Yi_method, 1, mean)/apply(Yi_method, 1, sd)),lower.tail = F))
cbind(apply(our_method, 1, mean),apply(our_method, 1, sd),
      apply(our_method, 1, mean)/apply(our_method, 1, sd),
      pnorm(abs(apply(our_method, 1, mean)/apply(our_method, 1, sd)),lower.tail = F))
#If no comparison is needed, simply using our method allows for 100 bootstrap iterations using the following function.
ALL_ACTG_model_ME = corrected_ridge_AH_with_bootstrap(Y,delta,W,Z,sigma_0_hat,2,bootstrap_times = 100, bootstrap_property = 1, critical.value = 0.1, trim = 0.16)