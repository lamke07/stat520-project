setwd("~/Downloads/GitHub/stat520-project/analysis")
library(here)
library(rjags)
library(coda)
rm(list = ls())

A_lcc_raw <- readRDS(here("data", "A_lcc.rds"))

A_lcc <- A_lcc_raw[1:15,1:15]
# A_lcc <- A_lcc_raw

N <- nrow(A_lcc)
K <- 5
# g_const <- log(choose(N-1, K-1))

lcc_dat <- list(A = A_lcc,
                N = N,
                K = K,
                # g_const = g_const,
                alpha = rep(1,K))

model.name <- "sbm.bugs"
sink(model.name)

cat("
model{
  # Priors
  for(i in 1:K){
    for(j in 1:K){
      Omega[i,j] ~ dunif(0,1)
    }
  }
  
  for(i in 1:N){
    g[i] ~ dcat(pi)
  }
  
  pi[1:K] ~ ddirch(alpha[1:K])
  # Likelihood
  for(i in 1:N){
    for(j in (i+1):N){
      A[i,j] ~ dbern(Omega[g[i],g[j]])
    }
  }
  
  
}
", fill = TRUE)

sink()

Omega = matrix(0.5, K,K)
g = sample(1:K, N, replace = TRUE)

inits = list(Omega = Omega, g = g)

jags.m <- jags.model(file = model.name, data = lcc_dat, inits = inits, n.chains = 1, n.adapt = 0)

params <- c("Omega", "g")

samps <- coda.samples(jags.m, params, n.iter = 1000)
