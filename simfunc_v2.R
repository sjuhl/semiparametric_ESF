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
sim_func <- function(spmultiplier, W, x, beta, dgp_type, ideal.setsize=F){
  ### SIMULATE TRUE DGP ###
  n <- length(x)
  
  # iid errors
  e <- rnorm(n,0,1)
  
  # DGP type
  if(dgp_type == "OLS"){
    y <- x*beta + e
  } else if(dgp_type == "SAR"){
    y <- spmultiplier %*% (x*b + e)
  } else if(dgp_type == "SEM"){
    y <- x*b + spmultiplier %*% e
  } else if(dgp_type == "SLX"){
    y <- x*b + W %*% x + e
  }
  
  ### ESTIMATION ###
  # 1) OLS
  ols <- lm(y ~ x)
  
  # 2) SAR
  sar <- lagsarlm(y ~ x, listw = mat2listw(W))
  
  # 3) SEM
  sem <- errorsarlm(y ~ x, listw = mat2listw(W))
  
  # 4) SLX
  WX <- W %*% x
  slx <- lm(y ~ x + WX)
  
  # 5) ESF
  esf_R2 <- lmFilter(y=y,x=x,W=W,positive=T,objfn="R2",ideal.setsize=ideal.setsize)
  esf_MI <- lmFilter(y=y,x=x,W=W,positive=T,objfn="MI",tol=.1,ideal.setsize=ideal.setsize)
  esf_p <- lmFilter(y=y,x=x,W=W,positive=T,objfn="p",sig=.05,bonferroni=T,ideal.setsize=ideal.setsize)
  esf_pMI <- lmFilter(y=y,x=x,W=W,positive=T,objfn="pMI",sig=.05,ideal.setsize=ideal.setsize)
  
  ### Output ###
  out <- data.frame(ols=coef(ols)["x"]
                    ,sar = coef(sar)["x"]
                    ,sem = coef(sem)["x"]
                    ,slx = coef(slx)["x"]
                    ,filtered_R2 = coef(esf_R2)["beta_1"]
                    ,filtered_p = coef(esf_p)["beta_1"]
                    ,filtered_MI = coef(esf_MI)["beta_1"]
                    ,filtered_pMI = coef(esf_pMI)["beta_1"]
                    ,se_ols = summary(ols)$coefficients["x","Std. Error"]
                    ,se_sar = sar$rest.se["x"]
                    ,se_sem = sem$rest.se["I(x - lambda * WX)x"]
                    ,se_slx = summary(slx)$coefficients["x","Std. Error"]
                    ,se_filtered_R2 = esf_R2$estimates["beta_1","SE"]
                    ,se_filtered_p = esf_p$estimates["beta_1","SE"]
                    ,se_filtered_MI = esf_MI$estimates["beta_1","SE"]
                    ,se_filtered_pMI = esf_pMI$estimates["beta_1","SE"]
                    ,fit_init = esf_R2$fit["Initial"]
                    ,fit_filtered_R2 = esf_R2$fit["Filtered"]
                    ,fit_filtered_p = esf_p$fit["Filtered"]
                    ,fit_filtered_MI = esf_MI$fit["Filtered"]
                    ,fit_filtered_pMI = esf_pMI$fit["Filtered"]
                    ,moran_init = esf_R2$moran["Initial","z"]
                    ,moran_filtered_R2 = esf_R2$moran["Filtered","z"]
                    ,moran_filtered_p = esf_p$moran["Filtered","z"]
                    ,moran_filtered_MI = esf_MI$moran["Filtered","z"]
                    ,moran_filtered_pMI = esf_pMI$moran["Filtered","z"]
                    ,nev_R2 = esf_R2$other$nev
                    ,nev_p = esf_p$other$nev
                    ,nev_MI = esf_MI$other$nev
                    ,nev_pMI = esf_pMI$other$nev
                    ,ncandidates = esf_MI$other$ncandidates
  )
  # return
  return(out)
}
