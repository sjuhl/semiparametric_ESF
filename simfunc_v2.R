#################################
# Simulation Function (Version 2)
#################################

# NOTE:
# This script contains the revised simulation functions used
# for the Monte Carlo study.

### DGPs:
# 1) Non-spatial
# 2) SAR
# 3) SEM
# 4) SLX
# 5) Spatial Heterogeneity (SLX)

### Models:
# 1) Non-spatial OLS
# 2) Parametric spatial regression models
#     - SAR
#     - SEM
#     - SLX
# 3) ESF (different objective functions - see 'spfilteR' package documentation)
#     - R2
#     - MI
#     - p
#     - pMI



# general simulation function
sim_func <- function(spmultiplier, W, x, beta, theta, dgp_type, ideal.setsize = FALSE){
  ### SIMULATE TRUE DGP ###
  n <- length(x)
  
  # iid errors
  e <- rnorm(n, 0, 1)
  
  # DGP type
  if(dgp_type == "OLS"){
    y <- x * beta + e
  } else if(dgp_type == "SAR"){
    y <- spmultiplier %*% (x * beta + e)
  } else if(dgp_type == "SEM"){
    y <- x * beta + spmultiplier %*% e
  } else if(dgp_type == "SLX"){
    y <- x * beta + theta * W %*% x + e
  } else if(dgp_type == 'HET'){
    theta1 <- theta - .5
    theta2 <- theta

    thetaW <- W
    
    thetaW[(nrow(thetaW)/2 + 1) : nrow(thetaW), (nrow(thetaW)/2 + 1) : nrow(thetaW)] <- theta1 * thetaW[(nrow(thetaW)/2 + 1) : nrow(thetaW), (nrow(thetaW)/2 + 1) : nrow(thetaW)]
    thetaW[1 : (nrow(thetaW)/2), 1 : (nrow(thetaW)/2)] <- theta1 * thetaW[1 : (nrow(thetaW)/2), 1 : (nrow(thetaW)/2)]

    thetaW[1 : (nrow(thetaW)/2), (nrow(thetaW)/2 + 1) : nrow(thetaW)] <- theta2 * thetaW[1 : (nrow(thetaW)/2), (nrow(thetaW)/2 + 1) : nrow(thetaW)]
    thetaW[(nrow(thetaW)/2 + 1) : nrow(thetaW), 1 : (nrow(thetaW)/2)] <- theta2 * thetaW[(nrow(thetaW)/2 + 1) : nrow(thetaW), 1 : (nrow(thetaW)/2)]

    y <- x * beta + thetaW %*% x + e
  }
  
  ### ESTIMATION ###
  # 1) OLS
  ols <- lm(y ~ x)
  lm_tests <- lm.LMtests(model = ols, listw = mat2listw(W, style = 'M')
                        ,zero.policy = TRUE, test = c('RLMerr', 'RLMlag'))
  names(lm_tests$RLMerr)
  lm_tests$RLMerr$statistic
  lm_tests$RLMerr$p.value
  
  # 2) SAR
  sar <- lagsarlm(y ~ x, listw = mat2listw(W, style = 'M'), zero.policy = TRUE)

  # 3) SEM
  sem <- errorsarlm(y ~ x, listw = mat2listw(W, style = 'M'), zero.policy = TRUE)
   
  # 4) SLX
  WX <- W %*% x
  slx <- lm(y ~ x + WX)
  
  # 5) ESF
  esf_R2 <- lmFilter(y = y, x = x, W = W, positive = TRUE, objfn = "R2", alpha = .1 
                     ,ideal.setsize = ideal.setsize)
#  esf_MI <- lmFilter(y = y, x = x, W = W, positive = TRUE, objfn = "MI", tol=.1
  esf_MI <- lmFilter(y = y, x = x, W = W, positive = TRUE, objfn = "MI", tol = .05, alpha = .1
                      ,ideal.setsize = ideal.setsize)
#  esf_p <- lmFilter(y = y, x = x, W = W, positive = TRUE, objfn = "p", sig = .05, bonferroni = TRUE
  esf_p <- lmFilter(y = y, x = x, W = W, positive = TRUE, objfn = "p", sig = .05, alpha = .1, bonferroni = FALSE
                    ,ideal.setsize = ideal.setsize)
  esf_pMI <- lmFilter(y = y, x = x, W = W, positive = TRUE, objfn = "pMI", sig = .05, alpha = .1
                    ,ideal.setsize = ideal.setsize)

  # 6) Moran's I
  m.ols <- MI.resid(resid = residuals(ols), x = x, W = W, alternative = "greater")
  m.sar <- MI.resid(resid = residuals(sar), x = x, W = W, alternative = "greater")
  m.sem <- MI.resid(resid = residuals(sem), x = x, W = W, alternative = "greater")
  m.slx <- MI.resid(resid = residuals(slx), x = x, W = W, alternative = "greater")
  
  ### Output ###
  out <- data.frame(ols = coef(ols)["x"]
                    ,sar = coef(sar)["x"]
                    ,sem = coef(sem)["x"]
                    ,slx = coef(slx)["x"]
                    ,rho = coef(sar)["rho"]
                    ,lambda = coef(sem)["lambda"]
                    ,theta = coef(slx)["WX"]
                    ,filtered_R2 = coef(esf_R2)["beta_1"]
                    ,filtered_p = coef(esf_p)["beta_1"]
                    ,filtered_MI = coef(esf_MI)["beta_1"]
                    ,filtered_pMI = coef(esf_pMI)["beta_1"]
                    ,se_ols = summary(ols)$coefficients["x", "Std. Error"]
                    ,se_sar = sar$rest.se["x"]
                    ,se_sem = sem$rest.se["I(x - lambda * WX)x"]
                    ,se_slx = summary(slx)$coefficients["x", "Std. Error"]
                    ,se_rho = sar$rho.se
                    ,se_lambda = sem$lambda.se
                    ,se_theta = summary(slx)$coefficients["WX", "Std. Error"]
                    ,se_filtered_R2 = esf_R2$estimates["beta_1", "SE"]
                    ,se_filtered_p = esf_p$estimates["beta_1", "SE"]
                    ,se_filtered_MI = esf_MI$estimates["beta_1", "SE"]
                    ,se_filtered_pMI = esf_pMI$estimates["beta_1", "SE"]
                    ,RLMlag_test = lm_tests$RLMlag$statistic
                    ,RLMlag_test_p = lm_tests$RLMlag$p.value
                    ,RLMerr_test = lm_tests$RLMerr$statistic
                    ,RLMerr_test_p = lm_tests$RLMerr$p.value
                    ,fit_init = esf_R2$fit["Initial"]
                    ,fit_filtered_R2 = esf_R2$fit["Filtered"]
                    ,fit_filtered_p = esf_p$fit["Filtered"]
                    ,fit_filtered_MI = esf_MI$fit["Filtered"]
                    ,fit_filtered_pMI = esf_pMI$fit["Filtered"]
                    ,fit_ols = summary(ols)$adj.r.squared
                    ,fit_slx = summary(slx)$adj.r.squared
                    ,moran_ols = m.ols$zI
                    ,moran_sar = m.sar$zI
                    ,moran_sem = m.sem$zI
                    ,moran_slx = m.slx$zI
                    ,moran_init = esf_R2$moran["Initial", "z"]
                    ,moran_filtered_R2 = esf_R2$moran["Filtered", "z"]
                    ,moran_filtered_p = esf_p$moran["Filtered", "z"]
                    ,moran_filtered_MI = esf_MI$moran["Filtered", "z"]
                    ,moran_filtered_pMI = esf_pMI$moran["Filtered", "z"]
                    ,nev_R2 = esf_R2$other$nev
                    ,nev_p = esf_p$other$nev
                    ,nev_MI = esf_MI$other$nev
                    ,nev_pMI = esf_pMI$other$nev
                    ,ncandidates = esf_MI$other$ncandidates
  )
  # return
  return(out)
}
