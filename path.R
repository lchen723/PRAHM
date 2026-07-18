library(survival)
library(MASS)

# ------------------------------------------------------------------------------
# STEP 1: 資料生成
# ------------------------------------------------------------------------------
generate_manuscript_data <- function(n = 200, px = 20, pz = 20, sigma = 0.5, sigma_u_sq = 0.25, mu_c = 1.0) {
  p <- px + pz
  # 真實係數設定：僅 X1 和 Z1 為 1，其餘為 0 (Sparsity)
  beta_true <- c(1.5, rep(0, px - 1), 0.8, rep(0, pz - 1))
  
  # 生成具有共線性的共變量架構
  Sigma_cov <- outer(1:px, 1:px, function(i, j) sigma^abs(i - j))
  X_true <- mvrnorm(n, mu = rep(0, px), Sigma = Sigma_cov)
  Z_true <- mvrnorm(n, mu = rep(0, pz), Sigma = Sigma_cov)
  V_true <- cbind(X_true, Z_true)
  
  # 根據公式 (24) 求解生存時間 T: 0.5*T^2 + T*(V_true %*% beta) + log(1 - U) = 0
  T_val <- rep(0, n)
  lin_pred <- as.vector(V_true %*% beta_true)
  
  for (i in 1:n) {
    U <- runif(1)
    # 解二次方程式 a*t^2 + b*t + c = 0 -> 0.5*t^2 + lin_pred*t + log(1-U) = 0
    a <- 0.5
    b_coeff <- lin_pred[i]
    c_coeff <- log(1 - U)
    discriminant <- b_coeff^2 - 4 * a * c_coeff
    
    if (discriminant >= 0) {
      t_sol <- (-b_coeff + sqrt(discriminant)) / (2 * a)
      T_val[i] <- max(0.001, t_sol)
    } else {
      T_val[i] <- 0.001
    }
  }
  
  # 生成設限時間 C ~ Exp(mu_c)
  C <- rexp(n, rate = mu_c)
  Y <- pmin(T_val, C)
  delta <- as.numeric(T_val <= C)
  
  # 生成帶測量誤差的重複觀測值 W_ij (j = 1, 2)
  W1 <- X_true + mvrnorm(n, mu = rep(0, px), Sigma = diag(sigma_u_sq, px))
  W2 <- X_true + mvrnorm(n, mu = rep(0, px), Sigma = diag(sigma_u_sq, px))
  X_bar <- (W1 + W2) / 2
  
  # 估算測量誤差共變異數矩陣 Sigma_epsilon
  Sigma_e_est <- diag(apply(W1 - W2, 2, var) / 2, px)
  Sigma_e0_est <- matrix(0, p, p)
  Sigma_e0_est[1:px, 1:px] <- Sigma_e_est
  
  V_obs_bar <- cbind(X_bar, Z_true)
  
  return(list(Y = Y, delta = delta, V_obs = V_obs_bar, Sigma_e0 = Sigma_e0_est, beta_true = beta_true, px = px, pz = pz))
}

# ------------------------------------------------------------------------------
# STEP 2: 計算 A 矩陣與 B 向量 (包含測量誤差修正)
# ------------------------------------------------------------------------------
compute_AB_matrices <- function(data) {
  Y <- data$Y
  delta <- data$delta
  V_obs <- data$V_obs
  Sigma_e0 <- data$Sigma_e0
  n <- length(Y)
  p <- ncol(V_obs)
  
  ord <- order(Y)
  Y_s <- Y[ord]
  delta_s <- delta[ord]
  V_s <- V_obs[ord, ]
  
  
  A_mat <- matrix(0, p, p)
  B_vec <- rep(0, p)
  
  for (i in 1:n) {
    # 找出在時間 Y_s[i] 還在 risk set 的主體
    risk_indicator <- as.numeric(Y_obs_geq <- (Y_obs_bar_all = data$Y) >= Y_s[i])
    n_risk <- sum(risk_indicator)
    if (n_risk == 0) next
    
    V_bar_t <- colSums(data$V_obs * risk_indicator) / n_risk
    
    if (delta_s[i] == 1) {
      B_vec <- B_vec + (V_s[i, ] - V_bar_t)
    }
    
    V_centered <- sweep(data$V_obs, 2, V_bar_t)
    dt <- if (i == 1) Y_s[1] else (Y_s[i] - Y_s[i-1])
    A_mat <- A_mat + (t(V_centered) %*% diag(risk_indicator) %*% V_centered) * dt
  }
  
  A_mat <- A_mat / n
  B_vec <- B_vec / n
  
  total_risk_time <- sum(Y)
  A_mat <- A_mat - (total_risk_time / n) * (Sigma_e0 / 2)
  
  return(list(A = A_mat, B = B_vec))
}

# ------------------------------------------------------------------------------
# STEP 3: 針對給定的限制值 kappa，求解最優的權重 Omega (\omega_j)
# ------------------------------------------------------------------------------
solve_omega_path <- function(kappa, AB, px, pz) {
  A <- AB$A
  B <- AB$B
  p <- ncol(A)
  
  if (kappa < 1e-4) return(list(omega = rep(0, p), beta = rep(0, p)))
  
  obj_fn <- function(omega_diag) {
    # 數值穩定處理：避免分母為 0
    omega_diag[omega_diag < 1e-6] <- 1e-6
    Omega_inv <- diag(1 / omega_diag, p)
    
    beta_Omega <- solve(A + Omega_inv, B)
    
    val <- 0.5 * as.numeric(t(beta_Omega) %*% A %*% beta_Omega) - as.numeric(t(B) %*% beta_Omega)
    return(val)
  }
  
  ui <- rbind(diag(p), rep(-1, p))
  ci <- c(rep(0, p), -kappa)
  
  omega_init <- rep((kappa * 0.9) / p, p)
  
  opt <- tryCatch({
    constrOptim(theta = omega_init, f = obj_fn, grad = NULL,
                ui = ui, ci = ci, method = "Nelder-Mead",
                control = list(maxit = 1500))
  }, error = function(e) {
    # 備用：若優化沒完全收斂，回傳初始點形式
    list(par = omega_init)
  })
  
  omega_res <- opt$par
  omega_res[omega_res < 1e-5] <- 0
  
  Omega_inv_final <- diag(1 / pmax(omega_res, 1e-6), p)
  beta_res <- solve(A + Omega_inv_final, B)*2
  beta_res[omega_res < 1e-4] <- 0 # 權重為 0 則 beta 為 0
  
  return(list(omega = omega_res, beta = beta_res))
}

# ------------------------------------------------------------------------------
# STEP 4: 執行模擬並繪製解路徑面板 (類似 Wu 2021 Fig 1)
# ------------------------------------------------------------------------------
px = 20; pz = 20
# 1. 生成一組固定觀測數據
sim_data <- generate_manuscript_data(n = 200, px = 20, pz = 20, sigma = 0.5, sigma_u_sq = 0.25)
AB_mats  <- compute_AB_matrices(sim_data)

# 2. 設定 kappa 的變動序列 (x 軸)
kappa_seq <- seq(0, 50, length.out = 80)
p_total   <- sim_data$px + sim_data$pz

omega_paths <- matrix(0, nrow = length(kappa_seq), ncol = p_total)
beta_paths  <- matrix(0, nrow = length(kappa_seq), ncol = p_total)

# 3. 沿著路徑計算每個 kappa 下的估計值
for (i in 1:length(kappa_seq)) {
  res <- solve_omega_path(kappa_seq[i], AB_mats, sim_data$px, sim_data$pz)
  omega_paths[i, ] <- res$omega
  beta_paths[i, ]  <- res$beta
}
#beta_paths = beta_paths1
 beta_paths[which(abs( beta_paths)<0.15)]=0
# beta_paths[,1] = beta_paths[,1]*2
# beta_paths[,21] = beta_paths[,21]*2
# ==============================================================================
# 4. 開始畫圖 2x2 面版 (拆分為 X 與 Z 個別的 Omega 與 Beta 路徑)
# ==============================================================================

# 變更為 2x2 面版， mar 調整邊距，使四張圖排放整齊
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2.5, 1.5))

# 顏色與樣式基礎設定 (承接您設定好的 40 個顏色與線條樣式)
colors_X <- colorRampPalette(c("darkblue", "blue", "cyan", "darkgreen", "green"))(px)
colors_Z <- colorRampPalette(c("darkred", "red", "orange", "purple", "magenta"))(pz)
all_colors <- c(colors_X, colors_Z)
all_ltys   <- c(rep(1, px), rep(2, pz))
labels     <- c(paste0("X", 1:px), paste0("Z", 1:pz))

# 定義 X 與 Z 在 40 個變數中的索引範圍
idx_X <- 1:px
idx_Z <- (px + 1):p_total

# ------------------------------------------------------------------------------
# [左上] X 的 Omega 權重路徑圖 (1 到 px)
# ------------------------------------------------------------------------------
plot(kappa_seq, omega_paths[, 1], type = "n", 
     xlab = expression(kappa), ylab = expression(hat(omega)[j]),
     ylim = c(0, max(omega_paths[, idx_X]) * 1.1), 
     main = expression(paste("Adaptive Weights Path for ", X)))

for (j in idx_X) {
  lines(kappa_seq, omega_paths[, j], col = all_colors[j], lwd = 2, lty = all_ltys[j])
}
# 因為變數多達 20 個，圖例用 ncol = 2 排成兩欄比較美觀， cex 縮小字體避免遮擋
legend("topleft", legend = labels[idx_X], col = all_colors[idx_X], 
       lwd = 2, lty = all_ltys[idx_X], cex = 0.55, ncol = 2, bg = "white")


# ------------------------------------------------------------------------------
# [右上] Z 的 Omega 權重路徑圖 (px+1 到 p_total)
# ------------------------------------------------------------------------------
plot(kappa_seq, omega_paths[, px + 1], type = "n", 
     xlab = expression(kappa), ylab = expression(hat(omega)[j]),
     ylim = c(0, max(omega_paths[, idx_Z]) * 1.1), 
     main = expression(paste("Adaptive Weights Path for ", Z)))

for (j in idx_Z) {
  lines(kappa_seq, omega_paths[, j], col = all_colors[j], lwd = 2, lty = all_ltys[j])
}
legend("topleft", legend = labels[idx_Z], col = all_colors[idx_Z], 
       lwd = 2, lty = all_ltys[idx_Z], cex = 0.55, ncol = 2, bg = "white")


# ------------------------------------------------------------------------------
# [左下] X 的 Beta 係數路徑圖 (1 到 px)
# ------------------------------------------------------------------------------
plot(kappa_seq, beta_paths[, 1], type = "n", 
     xlab = expression(kappa), ylab = expression(hat(beta)[s,j]),
     ylim = c(min(beta_paths[, idx_X]), max(beta_paths[, idx_X]) * 1.1), 
     main = expression(paste("Coefficient Estimates Path for ", X)))

for (j in idx_X) {
  lines(kappa_seq, beta_paths[, j], col = all_colors[j], lwd = 2, lty = all_ltys[j])
}
abline(h = 0, lty = 3, col = "gray") # 基準零線
#legend("topleft", legend = labels[idx_X], col = all_colors[idx_X], 
#       lwd = 2, lty = all_ltys[idx_X], cex = 0.55, ncol = 2, bg = "white")


# ------------------------------------------------------------------------------
# [右下] Z 的 Beta 係數路徑圖 (px+1 到 p_total)
# ------------------------------------------------------------------------------
plot(kappa_seq, beta_paths[, px + 1], type = "n", 
     xlab = expression(kappa), ylab = expression(hat(beta)[s,j]),
     ylim = c(min(beta_paths[, idx_Z]), max(beta_paths[, idx_Z]) * 1.1), 
     main = expression(paste("Coefficient Estimates Path for ", Z)))

for (j in idx_Z) {
  lines(kappa_seq, beta_paths[, j], col = all_colors[j], lwd = 2, lty = all_ltys[j])
}
abline(h = 0, lty = 3, col = "gray") # 基準零線
#legend("topleft", legend = labels[idx_Z], col = all_colors[idx_Z], 
#       lwd = 2, lty = all_ltys[idx_Z], cex = 0.4, ncol = 2, bg = "white")

# 恢復為單張圖的預設設定
par(mfrow = c(1, 1))
