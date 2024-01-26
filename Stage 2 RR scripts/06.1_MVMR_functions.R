#06.1_MVMR_functions.R

#make function for depression MVMR across outcomes
perform_mvmr_analysis <- function(outcome_filename) {
  #read in the outcome data
  outcome_data <- read_outcome_data(
    snps = Exposures_H4$SNP,
    filename = outcome_filename,
    sep = " ",
    snp_col = "ID",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "ALLELE1",
    other_allele_col = "ALLELE0",
    pval_col = "P"
  )
  
  #harmonize the exposures and outcome
  mv_data <- harmonise_data(Exposures_H4, outcome_data, action = 1)
  
  #keep only mr_keep = TRUE
  mv_data <- mv_data[mv_data$mr_keep == TRUE, ]
  
  #format for MVMR package and run analysis
  bX1 <- c(mv_data$beta.exposure[mv_data$id.exposure == 1])
  bX2 <- c(mv_data$beta.exposure[mv_data$id.exposure == 2])
  bY <- c(mv_data$beta.outcome[mv_data$id.exposure == 1])
  bYse <- c(mv_data$se.outcome[mv_data$id.exposure == 1])
  bXse1 <- c(mv_data$se.exposure[mv_data$id.exposure == 1])
  bXse2 <- c(mv_data$se.exposure[mv_data$id.exposure == 2])
  
  df <- data.frame(bX1, bXse1, bX2, bXse2, bY, bYse)
  df_mvmr <- format_mvmr(df[, c(1, 3)], df[, 5], df[, c(2, 4)], df[, 6])
  
  res <- ivw_mvmr(df_mvmr) %>%
    as.data.frame() %>%
    mutate(lci = Estimate - `Std. Error` * 1.96,
           uci = Estimate + `Std. Error` * 1.96)
  
  #include F-statistic and pleiotropy tests
  Fstat <- strength_mvmr(df_mvmr, gencov=0)
  ptr <- pleiotropy_mvmr(df_mvmr, gencov=0)
  
  #reformat data for MVMR-Egger sensitivity analysis
  mvmr_egger_input <- mr_mvinput(
    bx = cbind(bX1, bX2),
    by = bY,
    bxse = cbind(bXse1, bXse2),
    byse = bYse)
  
  #MVMR-Egger
  mvmr_egger <- mr_mvegger(
    mvmr_egger_input,
    orientate = 1)
  
  #reformat data for further MVMR sensitivity analyses
  bx = as.matrix(df_mvmr[, c("betaX1", "betaX2")])  # Exposure effect sizes
  sebx = as.matrix(df_mvmr[, c("sebetaX1", "sebetaX2")])  # Exposure standard errors
  by = df_mvmr$betaYG  # Outcome effect sizes
  seby = df_mvmr$sebetaYG  # Outcome standard errors
  
  #use mvmr_median function to obtain results
  mvmr_median_results = mvmr_median(bx, sebx, by, seby, boot_it = 1000)
  
  #format results and calculate CIs
  mvmr_median <- mvmr_median_results %>%
    as.data.frame() %>%
    mutate(lci = coefficients - se*1.96,
           uci = coefficients + se*1.96)
  
  #use mvmr_lasso function to run lasso
  mvmr_lasso_results = mvmr_lasso(bx, by, seby)
  
  #format results (estimates and standard errors)
  mvmr_lasso_ests <- as.data.frame(mvmr_lasso_results$th_post)
  colnames(mvmr_lasso_ests)[1] <- "coefficients"
  
  mvmr_lasso_ses <- as.data.frame(mvmr_lasso_results$se_post)
  colnames(mvmr_lasso_ses)[1] <- "se"
  
  #merge estimates and SEs, and calculate CIs
  mvmr_lasso <- mvmr_lasso_ests %>%
    cbind(mvmr_lasso_ses) %>%
    mutate(lci = coefficients - se*1.96,
           uci = coefficients + se*1.96)
  
  #return a list of results
  return(list(mvmr_results = res, F_statistic = Fstat, pleiotropy_test = ptr,
              mvmr_egger = mvmr_egger, mvmr_median = mvmr_median, mvmr_lasso = mvmr_lasso))
}

##MVMR sensitivity functions

#MVMR-Median
mvmr_med_boot = function(bx, sebx, by, seby, N){
  est = sapply(1:N, function(i){
    p = length(by)
    k = dim(bx)[2]
    Sx = lapply(1:p, function(j){diag(sebx[j, ]^2)})
    bxboot = sapply(1:p, function(j){mvrnorm(1, bx[j, ], Sx[[j]])})
    bxboot = t(bxboot)
    byboot = rnorm(p, by, seby)
    rq(byboot ~ bxboot - 1, weights = seby^-2)$coefficients
  })
  apply(est, 1, sd)
}

mvmr_median = function(bx, sebx, by, seby, boot_it = 1000){
  qr_mod = rq(by ~ bx - 1, weights = seby^-2)
  boot_se = mvmr_med_boot(bx, sebx, by, seby, boot_it)
  return(list("coefficients" = qr_mod$coefficients, "se" = boot_se))
}

#MVMR-Lasso
cv.mvmr_lasso = function(bx, by, seby){
  p = dim(bx)[1]
  k = dim(bx)[2]
  S = diag(seby^-2)
  b = S^(1/2) %*% bx
  Pb = b %*% solve(t(b) %*% b, t(b))
  xlas = (diag(p) - Pb) %*% S^(1/2)
  ylas = (diag(p) - Pb) %*% S^(1/2) %*% by
  alas = glmnet(xlas, ylas, intercept = FALSE)
  lamseq = sort(alas$lambda)
  lamlen = length(lamseq)
  rse = sapply(1:lamlen, function(j){
    av = which(alas$beta[, (lamlen - j + 1)] == 0)
    mod = lm.fit(as.matrix(S[av, av]^(1/2) %*% bx[av, ]), S[av, av]^(1/2) %*% by[av])
    c(sqrt(t(mod$residuals) %*% (mod$residuals) / (mod$df.residual)), length(av))
  })
  rse_inc = rse[1, 2:lamlen] - rse[1, 1:(lamlen-1)]
  het = which(rse[1, 2:lamlen] > 1 & rse_inc > ((qchisq(0.95, 1) / rse[2, 2:lamlen])))
  if (length(het) == 0){
    lam_pos = 1
  } else {
    lam_pos = min(het)
  }
  num_valid = rev(sapply(1:lamlen, function(j){sum(alas$beta[, j]==0)}))
  min_lam_pos = min(which(num_valid > k))
  if (lam_pos < min_lam_pos){lam_pos = min_lam_pos}
  return(list(fit = alas$beta[, (lamlen - lam_pos + 1)], lambda = lamseq[lam_pos]))
}

mvmr_lasso = function(bx, by, seby){
  p = dim(as.matrix(bx))[1]
  k = dim(as.matrix(bx))[2]
  S = diag(seby^-2)
  sn = sign(bx[, 1])
  bx_or = bx * sn
  by_or = by * sn
  cv.alas = cv.mvmr_lasso(bx_or, by_or, seby)
  a1 = cv.alas$fit
  e = by_or - a1
  thest = solve(t(bx_or) %*% S %*% bx_or, t(bx_or) %*% S %*% e)
  v = which(a1==0)
  mvmr_mod = mr_mvivw(mr_mvinput(bx = bx_or[v, ], bxse = bx_or[v, ],
                                 by = by_or[v], byse = seby[v]))
  th_post = mvmr_mod$Estimate
  se_post = mvmr_mod$StdError
  return(list(thest = thest, a = a1, lambda = cv.alas$lambda,
              th_post = th_post, se_post = se_post))
}

#MVMR-Median2
mvmr_med_boot2 = function(bx2, sebx2, by2, seby2, N){
  est = sapply(1:N, function(i){
    p = length(by2)
    k = dim(bx2)[2]
    Sx = lapply(1:p, function(j){diag(sebx2[j, ]^2)})
    bxboot = sapply(1:p, function(j){mvrnorm(1, bx2[j, ], Sx[[j]])})
    bxboot = t(bxboot)
    byboot = rnorm(p, by2, seby2)
    rq(byboot ~ bxboot - 1, weights = seby2^-2)$coefficients
  })
  apply(est, 1, sd)
}


mvmr_median2 = function(bx2, sebx2, by2, seby2, boot_it = 1000){
  qr_mod = rq(by2 ~ bx2 - 1, weights = seby2^-2)
  boot_se = mvmr_med_boot2(bx2, sebx2, by2, seby2, boot_it)
  return(list("coefficients" = qr_mod$coefficients, "se" = boot_se))
}

#MVMR-Lasso2
cv.mvmr_lasso2 = function(bx2, by2, seby2){
  p = dim(bx2)[1]
  k = dim(bx2)[2]
  S = diag(seby2^-2)
  b = S^(1/2) %*% bx2
  Pb = b %*% solve(t(b) %*% b, t(b))
  xlas = (diag(p) - Pb) %*% S^(1/2)
  ylas = (diag(p) - Pb) %*% S^(1/2) %*% by
  alas = glmnet(xlas, ylas, intercept = FALSE)
  lamseq = sort(alas$lambda)
  lamlen = length(lamseq)
  rse = sapply(1:lamlen, function(j){
    av = which(alas$beta[, (lamlen - j + 1)] == 0)
    mod = lm.fit(as.matrix(S[av, av]^(1/2) %*% bx[av, ]), S[av, av]^(1/2) %*% by[av])
    c(sqrt(t(mod$residuals) %*% (mod$residuals) / (mod$df.residual)), length(av))
  })
  rse_inc = rse[1, 2:lamlen] - rse[1, 1:(lamlen-1)]
  het = which(rse[1, 2:lamlen] > 1 & rse_inc > ((qchisq(0.95, 1) / rse[2, 2:lamlen])))
  if (length(het) == 0){
    lam_pos = 1
  } else {
    lam_pos = min(het)
  }
  num_valid = rev(sapply(1:lamlen, function(j){sum(alas$beta[, j]==0)}))
  min_lam_pos = min(which(num_valid > k))
  if (lam_pos < min_lam_pos){lam_pos = min_lam_pos}
  return(list(fit = alas$beta[, (lamlen - lam_pos + 1)], lambda = lamseq[lam_pos]))
}

mvmr_lasso2 = function(bx2, by2, seby2){
  p = dim(as.matrix(bx2))[1]
  k = dim(as.matrix(bx2))[2]
  S = diag(seby2^-2)
  sn = sign(bx2[, 1])
  bx_or = bx2 * sn
  by_or = by2 * sn
  cv.alas = cv.mvmr_lasso(bx_or, by_or, seby2)
  a1 = cv.alas$fit
  e = by_or - a1
  thest = solve(t(bx_or) %*% S %*% bx_or, t(bx_or) %*% S %*% e)
  v = which(a1==0)
  mvmr_mod = mr_mvivw(mr_mvinput(bx = bx_or[v, ], bxse = bx_or[v, ],
                                 by = by_or[v], byse = seby2[v]))
  th_post = mvmr_mod$Estimate
  se_post = mvmr_mod$StdError
  return(list(thest = thest, a = a1, lambda = cv.alas$lambda,
              th_post = th_post, se_post = se_post))
}
