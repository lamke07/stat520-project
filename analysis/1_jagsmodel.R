setwd("~/Downloads/GitHub/stat520-project/analysis")
library(here)
library(rjags)
library(coda)
rm(list = ls())

# Setup of Stochastic Block Model

model.name <- "1_sbm.bugs"
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
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      A[i,j] ~ dbern(Omega[g[i],g[j]])
    }
  }
  
  
}
", fill = TRUE)
sink()
