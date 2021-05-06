##########################
# Simulation Function
##########################

sim_func <- function(spmultiplier,x,beta,sd.e,W,SEM=F,ideal.setsize=F){
  ### SIMULATE TRUE DGP ###
  n <- nrow(spmultiplier)
  # iid errors
  e <- rnorm(n,0,sd.e)
  # spatial DGPs
  if(SEM){
    y <- x*b + spmultiplier %*% e
  } else y <- spmultiplier %*% (x*b + e)
  
  ### ESTIMATION ###
  # ols
  ols <- lm(y ~ x)
  # filtering
  esf_R2 <- lmFilter(y=y,x=x,W=W,positive=T,objfn="R2",ideal.setsize=ideal.setsize)
  esf_MI <- lmFilter(y=y,x=x,W=W,positive=T,objfn="MI",tol=.1,ideal.setsize=ideal.setsize)
  esf_p <- lmFilter(y=y,x=x,W=W,positive=T,objfn="p",sig=.05,bonferroni=T,ideal.setsize=ideal.setsize)
  esf_pMI <- lmFilter(y=y,x=x,W=W,positive=T,objfn="p",sig=.05,ideal.setsize=ideal.setsize)
  
  ### Output ###
  out <- data.frame(ols=summary(ols)$coefficients["x","Estimate"]
                    ,filtered_R2=esf_R2$estimates["beta_1","Estimate"]
                    ,filtered_p=esf_p$estimates["beta_1","Estimate"]
                    ,filtered_MI=esf_MI$estimates["beta_1","Estimate"]
                    ,filtered_pMI=esf_pMI$estimates["beta_1","Estimate"]
                    ,se_ols=summary(ols)$coefficients["x","Std. Error"]
                    ,se_filtered_R2=esf_R2$estimates["beta_1","SE"]
                    ,se_filtered_p=esf_p$estimates["beta_1","SE"]
                    ,se_filtered_MI=esf_MI$estimates["beta_1","SE"]
                    ,se_filtered_pMI=esf_pMI$estimates["beta_1","SE"]
                    ,fit_init=esf_R2$fit["Initial"]
                    ,fit_filtered_R2=esf_R2$fit["Filtered"]
                    ,fit_filtered_p=esf_p$fit["Filtered"]
                    ,fit_filtered_MI=esf_MI$fit["Filtered"]
                    ,fit_filtered_pMI=esf_pMI$fit["Filtered"]
                    ,moran_init=esf_R2$moran["Initial","z"]
                    ,moran_filtered_R2=esf_R2$moran["Filtered","z"]
                    ,moran_filtered_p=esf_p$moran["Filtered","z"]
                    ,moran_filtered_MI=esf_MI$moran["Filtered","z"]
                    ,moran_filtered_pMI=esf_pMI$moran["Filtered","z"]
                    ,nev_R2=esf_R2$other$nev
                    ,nev_p=esf_p$other$nev
                    ,nev_MI=esf_MI$other$nev
                    ,nev_pMI=esf_pMI$other$nev
                    ,ncandidates=esf_MI$other$ncandidates
  )
  # return
  return(out)
}


