##########################
# Monte Carlo Simulations
##########################

# clean environment
rm(list=ls())

### load packages
# spatial packages
library(spdep) # for W
library(spfilteR) # for filtering

# parallel computing
library(foreach)
library(doParallel)
library(doRNG)

# read simulation function
source("./simfunc.R")

# load matrices
load("./W.RData")

# spectral normalization
z.W_states <- W_states/as.numeric(eigen(W_states)$values[1])
z.W_cities <- W_cities/as.numeric(eigen(W_cities)$values[1])
z.W_cd116 <- W_cd116/as.numeric(eigen(W_cd116)$values[1])

# check
all.equal(target=1, current=as.numeric(eigen(z.W_states)$values[1]))
all.equal(target=1, current=as.numeric(eigen(z.W_cities)$values[1]))
all.equal(target=1, current=as.numeric(eigen(z.W_cd116)$values[1]))



###################
# MONTE CARLO SIMULATION

# parameters
b <- 1
sigma2 <- c(1,2) # variance of disturbances
SEM <- c(TRUE,FALSE)

# calculate spatial multipliers
p <- c(0,.25,.5,.75)
W <- list(z.W_states,z.W_cities,z.W_cd116)
W_id <- c(1,2,3)

grid <- matrix(as.matrix(cbind(expand.grid(p,W_id,sigma2,SEM),seq_len(length(p)*length(W)))),ncol=5
               ,dimnames=list(NULL, c("p", "W_id","sigma2","SEM","multi_id")))
#unique(grid[grid[,"multi_id"]==6,c("p","W_id")])
multipliers <- apply(grid,1,function(x) solve(diag(1,nrow(W[[x[2]]]))-x[1]*W[[x[2]]]))

# cross-sections
n <- sapply(W, nrow)

# covariates
set.seed(123)
covars <- list()
for(i in seq_along(n)){
  covars[[i]] <- rnorm(n[i],0,1)
}

ninput <- nrow(grid)
nsim <- 1000

# specify function inputs
input <- data.frame(do.call(rbind,replicate(nsim,grid,simplify=F)))
input <- input[order(input$W_id,input$p,input$multi_id,input$sigma2,input$SEM),]
input$SEM <- input$SEM==1
unique(input)
rownames(input) <- 1:nrow(input)

# add n to input matrix
input$n <- n[input$W_id]

# check
nrow(input)==ninput*nsim


# test
#sim_func(spmultiplier=multipliers[[input$multi_id[1]]],x=covars[[input$W_id[1]]]
#         ,beta=b,sigma2=input$sigma2[1],W=W[[input$W_id[1]]],SEM=input$SEM[1]
#         ,ideal.setsize=F)


# parallel computing
(ncores <- detectCores())
nworkers <- ncores - 1
cl <- makeCluster(nworkers)
registerDoParallel(cl)
registerDoRNG(12345)

# start timer
start.time <- Sys.time()
# simulate
sim_out <- foreach(i=1:nrow(input), .combine=rbind,
                   .packages=c("spdep","spfilteR")
                   ) %dopar% {
                     sim_func(spmultiplier=multipliers[[input$multi_id[i]]],x=covars[[input$W_id[i]]]
                              ,beta=b,sigma2=input$sigma2[i],W=W[[input$W_id[i]]],SEM=input$SEM[i]
                              ,ideal.setsize=F)
                     }
stopCluster(cl)

(time.taken <- Sys.time() - start.time)
session_info <- sessionInfo()


# create empty folder (if not existing already)
ifelse(!dir.exists(file.path("./SimOut"))
       ,dir.create(file.path("./SimOut")), FALSE)

# save output
save(input,sim_out,covars,nsim,time.taken,session_info,file="./SimOut/MC_Out.RData")

