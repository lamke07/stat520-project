library(here)
library(tidyverse)
library(stringi)
library(MCMCprecision)
rm(list = ls())
res_output <- readRDS(here::here("analysis", "res_models.rds"))
A_lcc <- readRDS(here("data", "A_lcc_16.rds"))

rownames(A_lcc) <- NULL
colnames(A_lcc) <- NULL

res_output_results <- purrr::map(res_output, "result")

output1 <- res_output_results[[1]]

# Use this
names(output1$BUGSoutput$sims.list)
nrow(output1$BUGSoutput$sims.list$g)
names(output1$BUGSoutput)
rownames(output1$BUGSoutput$summary)

################################################################################
# Model Summaries

get_jags_summary <- function(bugs_summary){
  bugs_summary %>%
    as_tibble() %>%
    mutate(parameter = rownames(bugs_summary),
           par_type = gsub("\\[.*","",parameter),
           par_ind1 = stri_extract_first_regex(parameter, "[0-9]+"),
           par_ind2 = stri_extract_last_regex(parameter, "[0-9]+")) %>%
    return()
}

res_output_summary <- purrr::map_df(1:5, ~get_jags_summary(res_output_results[[.x]]$BUGSoutput$summary) %>%
                                      mutate(run = .x,
                                             Rhat_conv = as.numeric(Rhat <= 1.05)))

res_output_deviance <- purrr::map_dfc(1:5, ~ as.vector(res_output_results[[.x]]$BUGSoutput$sims.list$deviance))
colnames(res_output_deviance) <- paste0("model", 1:5)

################################################################################
# Traceplots

p <- res_output_deviance %>%
  mutate(sample = row_number()) %>%
  pivot_longer(-c("sample"), names_to = "Model", values_to = "Deviance") %>%
  ggplot() +
  geom_line(aes(x = sample, y = Deviance)) +
  facet_wrap(~Model) +
  labs(x = "Sample", y = "Deviance", title = "Traceplot of Deviance") +
  theme_light()

dim(res_output_results[[1]]$BUGSoutput$sims.list$Omega)

res_output_results[[2]]$BUGSoutput$sims.matrix %>%
  as_tibble() %>%
  dplyr::select(contains("Omega")) %>% 
  mutate(sample_nr = row_number()) %>%
  pivot_longer(-c("sample_nr"), names_to = "Parameter", values_to = "par_value") %>%
  ggplot() +
  geom_line(aes(x = sample_nr, y = par_value)) +
  facet_wrap(~Parameter)

res_output_results[[2]]$BUGSoutput$sims.matrix %>%
  as_tibble() %>%
  dplyr::select(starts_with("g")) %>% 
  mutate(sample_nr = row_number()) %>%
  pivot_longer(-c("sample_nr"), names_to = "Parameter", values_to = "par_value") %>%
  mutate(g_node = as.numeric(stri_extract_first_regex(Parameter, "[0-9]+"))) %>%
  filter(g_node <= 25) %>%
  ggplot() +
  geom_line(aes(x = sample_nr, y = par_value)) +
  facet_wrap(~g_node) +
  theme_light()


################################################################################
# DIC
res_output_DIC <- purrr::map_dbl(1:5, ~res_output_results[[.x]]$BUGSoutput$DIC)

################################################################################
# Bayes factor
K = 5
N = 97


sample_prior <- function(N, K){
  Omega <- matrix(runif(K*K, 0, 1), nrow = K, ncol = K)
  pi <- MCMCprecision::rdirichlet(1, rep(1,K))
  g <- purrr::map_int(1:N, ~sample(1:K, size = 1, replace = TRUE, prob = pi))
  
  return(list(Omega = Omega, pi = pi, g = g))
}

generate_sample_prior_lik <- function(A, N, K){
  prior_pars <- sample_prior(N, K)
  res <- sample_lik(A = A , Omega = prior_pars$Omega, g = prior_pars$g)
}

sample_lik <- function(A, Omega, g){
  res = 0
  N = nrow(A)
  
  for(i in (1:(N-1))){
    for(j in (i+1):N){
      res = res + dbinom(A[i,j], size = 1, prob = Omega[g[i], g[j]], log = TRUE)
    }
  }
  return(res)
}

model_lik <- purrr::map_dbl(1:500, ~generate_sample_prior_lik(A = A_lcc, N = nrow(A_lcc), K = 4))

################################################################################
clusters <- res_output_results[[5]]$BUGSoutput$sims.list$g

compute_psm <- function(cluster_output){
  # Compute posterior similarity matrix for a cluster output
  p_psm <- ncol(cluster_output)
  
  psm <- matrix(0, nrow = p_psm, ncol = p_psm)
  
  for(i in 1:p_psm){
    for(j in i:p_psm){
      psm[i,j] <- sum(cluster_output[,i] == cluster_output[,j])
    }
  }
  return(psm/nrow(cluster_output))
}

plot_psm_heatmap <- function(psm_matrix){
  df_plot <- as_tibble(psm_matrix)
  colnames(df_plot) <- seq(1,nrow(psm_matrix))
  
  df_plot <- df_plot %>%
    mutate(x = row_number()) %>%
    pivot_longer(-x, names_to = "y", values_to = "psm_val") %>%
    mutate(y = as.integer(y))
  
  p <- ggplot(df_plot, aes(x, y, fill= psm_val)) + 
    geom_tile() +
    scale_fill_viridis_c() +
    labs(fill = "PSM") +
    theme_light()
  
  return(p)
}

compute_dahl_loss <- function(psm, proposed_clusters){
  # Compute Dahl score for a proposed clustering
  g <- as.vector(proposed_clusters)
  p_psm <- nrow(psm)
  
  res <- 0
  for(i in 1:p_psm){
    for(j in i:p_psm){
      res = res + (as.numeric(g[i] == g[j]) - psm[i,j])^2
    }
  }
  return(res)
}

# Get final clustering
psm <- compute_psm(clusters)
p <- plot_psm_heatmap(psm)
x <- purrr::map_dbl(1:nrow(clusters), ~compute_dahl_loss(psm = psm, proposed_clusters = clusters[.x,]))

final_clusters <- as.vector(clusters[which.min(x),])

################################################################################

