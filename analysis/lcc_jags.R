library(here)
library(rjags)
library(coda)

A_lcc <- readRDS(here("data", "A_lcc.rds"))

lcc_dat <- list(A = A_lcc,
                n = nrow(A_lcc),
                K = 30)

fit <- stan(file = 'lcc_stan.stan', data = lcc_dat)
