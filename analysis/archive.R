# A_network <- network(A_edges_extended, directed = FALSE)
# 
# # A_network <- network(A_edges_extended, vertex.attr = A_nodes, matrix.type = "edgelist", ignore.eval = FALSE)
# A_adj <- A_edges_extended %>%
#   dplyr::select(from, to, weight) %>%
#   pivot_wider(names_from = to, values_from = weight)

#### if you want to use the rjags package 
# system.time({
#   jags.m <- jags.model(file = model.name, data = lcc_dat, inits = inits, n.chains = 1, n.adapt = 0)
# })
# params <- c("Omega", "g")
# 
# samps <- coda.samples(jags.m, params, n.iter = 1000)

# Old R2jags
# # Set up data file
# N <- nrow(A_lcc)
# K <- 5
# alpha <- rep(1,K)
# 
# lcc_dat <- list(A = A_lcc,
#                 N = N,
#                 K = K,
#                 alpha = alpha)
# 
# ################################################################################
# # Setup initial parameters
# Omega = matrix(0.5, K,K)
# g = sample(1:K, N, replace = TRUE)
# 
# inits = list(Omega = Omega, g = g)
# inits_all = list(inits, inits, inits)
# bayes.mod.params <- c("Omega", "g")
# ################################################################################
# # Compile the model
# 
# 
# model.name <- "1_sbm.bugs"
# 
# #### if you want to use the R2jags package
# # Set seed
# set.seed(123)
# jags.m <- R2jags::jags(model.file = model.name, 
#                        data = lcc_dat, inits = inits_all, parameters.to.save = bayes.mod.params,
#                        n.chains = 3, n.iter = 9000, n.burnin = 1000)

# for(i in 1:97){
#   print(table(clusters[,i]))
# }

