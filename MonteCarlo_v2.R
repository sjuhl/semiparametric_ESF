##########################
# Monte Carlo Simulations
##########################

# clean environment
rm(list=ls())

### load packages
# spatial packages
library(spdep) # for W
library(spatialreg)
library(spfilteR) # for filtering

# parallel computing
library(foreach)
library(doParallel)
library(doRNG)
library(doSNOW) # for the progress bar

# read simulation function
source("./simfunc_v2.R")

# load matrices
load("./W.RData")

# spectral normalization
z.W_states <- W_states / as.numeric(eigen(W_states)$values[1])
z.W_cities <- W_cities / as.numeric(eigen(W_cities)$values[1])
z.W_cd116 <- W_cd116 / as.numeric(eigen(W_cd116)$values[1])

# check
all.equal(target = 1, current = as.numeric(eigen(z.W_states)$values[1]))
all.equal(target = 1, current = as.numeric(eigen(z.W_cities)$values[1]))
all.equal(target = 1, current = as.numeric(eigen(z.W_cd116)$values[1]))


###################
# MONTE CARLO SIMULATION

# parameters
b <- 10

# DGPs
dgp_type <- c("OLS", "SAR", "SEM", "SLX", 'HET')
dgp_id <- seq_along(dgp_type)

# calculate spatial multipliers
p <- c(0, .25, .5, .75)
W <- list(z.W_states, z.W_cities, z.W_cd116)
W_id <- c(1, 2, 3)

grid <- matrix(as.matrix(cbind(expand.grid(p, W_id, dgp_id), seq_len(length(p) * length(W)))), ncol = 4
               ,dimnames = list(NULL, c("p", "W_id", "dgp_type", "multi_id")))
#multipliers <- apply(grid, 1 ,function(x) solve(diag(1, nrow(W[[x[2]]])) - x[1] * W[[x[2]]]))

multipliers <- list()
for(i in seq_len(max(grid[,'multi_id']))){
  x <- unique(grid[grid[,'multi_id'] == i, c('p', 'W_id')])
  multipliers[[i]] <- solve(diag(1, nrow(W[[x[2]]])) - x[1] * W[[x[2]]])
}

# cross-sections
n <- sapply(W, nrow)

# covariates
set.seed(123)
covars <- list()
for(i in seq_along(n)){
  covars[[i]] <- rnorm(n[i], 0, 1)
}

ninput <- nrow(grid)
nsim <- 1000

# specify function inputs
input <- data.frame(do.call(rbind, replicate(nsim, grid, simplify = FALSE)))
input <- input[order(input$dgp_type, input$W_id, input$p, input$multi_id),]

# put correct dgp_type into dataframe
input$dgp_type <- dgp_type[input[,'dgp_type']]

# add n to input matrix
input$n <- n[input$W_id]

# check
nrow(input) == ninput*nsim

# remove rho > 0 when DGP type == 'OLS'
input <- input[!(input$dgp_type == 'OLS' & input$p > 0),]

# remove rho == 0 when DGP != 'OLS'
input <- input[!(input$dgp_type != 'OLS' & input$p == 0),]

# add correct rownames
unique(input)
rownames(input) <- seq_len(nrow(input))

# test
#sim_func(spmultiplier = multipliers[[input$multi_id[1]]], W = W[[input$W_id[1]]]
#                              ,x = covars[[input$W_id[1]]], beta = b, theta = input$p[1]
#                              ,dgp_type = input$dgp_type[1], ideal.setsize = FALSE)

# parallel computing
(ncores <- detectCores())
nworkers <- ncores - 1
cl <- makeCluster(nworkers)
registerDoSNOW(cl)
registerDoRNG(12345)

# start timer
start.time <- Sys.time()

# initialize progress bar
pb <- txtProgressBar(max = nrow(input), style = 3, char = '*')
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# simulate
sim_out <- foreach(i = seq_len(nrow(input)), .combine = rbind, .options.snow = opts,
                   .packages = c("spdep", "spatialreg", "spfilteR")
                   ) %dopar% {
                     # simulation function
                     sim_func(spmultiplier = multipliers[[input$multi_id[i]]], W = W[[input$W_id[i]]]
                              ,x = covars[[input$W_id[i]]], beta = b, theta = input$p[i]
                              ,dgp_type = input$dgp_type[i], ideal.setsize = FALSE)
                   }
close(pb)
stopCluster(cl)

(time.taken <- Sys.time() - start.time)
session_info <- sessionInfo()

# create empty folder (if not existing already)
ifelse(!dir.exists(file.path("./SimOut"))
       ,dir.create(file.path("./SimOut")), FALSE)

# save output
save(input, sim_out, covars, nsim, time.taken, session_info
    ,file="./SimOut/MC_Out_v2.RData")
