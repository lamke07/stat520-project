library(here)
library("rstan") # observe startup messages
options(mc.cores = parallel::detectCores()-1)
rstan_options(auto_write = TRUE)

A_lcc <- readRDS(here("data", "A_lcc.rds"))

lcc_dat <- list(A = A_lcc,
                n = nrow(A_lcc),
                K = 30)

fit <- stan(file = 'lcc_stan.stan', data = lcc_dat)
