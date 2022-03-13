setwd("~/Downloads/GitHub/stat520-project/analysis")
library(here)
library(rjags)
library(coda)
library(purrr)
rm(list = ls())

A_lcc <- readRDS(here("data", "A_lcc_16.rds"))

get_jags_model <- function(A_lcc, K, n.iter, n.burnin){
  ################################################################################
  # Set up data file
  N <- nrow(A_lcc)
  alpha <- rep(1,K)
  
  lcc_dat <- list(A = A_lcc,
                  N = N,
                  K = K,
                  alpha = alpha)
  
  ################################################################################
  # Setup initial parameters
  Omega = matrix(0.5, K,K)
  g = sample(1:K, N, replace = TRUE)
  
  inits = list(Omega = Omega, g = g)
  inits_all = list(inits, inits, inits)
  bayes.mod.params <- c("Omega", "g")
  ################################################################################
  # Compile the model
  
  
  model.name <- "1_sbm.bugs"
  
  #### if you want to use the R2jags package
  # Set seed
  set.seed(123+K)
  jags.m <- R2jags::jags(model.file = model.name, 
                         data = lcc_dat, inits = inits_all, parameters.to.save = bayes.mod.params,
                         n.chains = 3, n.iter = n.iter, n.burnin = n.burnin)
  
  return(jags.m)
}


# test_run <- get_jags_model(A_lcc = A_lcc, K = 5, n.iter = 90, n.burnin = 10)
safe_get_jags_model <- safely(.f = get_jags_model)
res_models <- purrr::map(seq(4,20,4), ~safe_get_jags_model(A_lcc = A_lcc, K = .x, n.iter = 9000, n.burnin = 1000))


# saveRDS(res_models, "res_models.rds")
################################################################################
