library(here)
res_output <- readRDS(here::here("analysis", "res_models.rds"))

res_output_results <- purrr::map(res_output, "result")

output1 <- res_output_results[[1]]
