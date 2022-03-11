library(here)
library(readr)
library(tidyverse)
# library(network)
library(igraph)

rm(list = ls())

raw_edges <- read_csv(here("data", "kelch_IBD_graph.csv"))
raw_nodes <- read_csv(here("data", "meta_data.csv"))

################################################################################
A_nodes <- raw_nodes %>%
  mutate(id = row_number()) %>%
  # mutate(id = paste0("V", row_number())) %>%
  rename(label = Sample) %>%
  relocate(id)

A_edges <- raw_edges %>%
  rename(source = V1, target = V2, weight = Probability_of_edge) %>%
  left_join(A_nodes %>%
              dplyr::select(label, from = id),
            by = c("source" = "label")) %>%
  left_join(A_nodes %>%
              dplyr::select(label, to = id),
            by = c("target" = "label"))

edges <- A_edges %>%
  dplyr::select(from, to, weight)

nodes <- A_nodes %>%
  dplyr::select(id, label)
################################################################################

get_LCC_size <- function(nodes, edges, threshold){
  edge_thresh <- edges %>%
    mutate(weight = as.integer(weight > threshold)) %>%
    filter(weight == 1)
  
  graph_thresh <- graph_from_data_frame(d = edge_thresh, vertices = nodes, directed = FALSE)

  return(max(components(graph_thresh, mode = "weak")$csize))
}

lcc_sizes <- purrr::map_dbl(1-2^(-seq(1:20)), ~get_LCC_size(nodes, edges, threshold = .x))

lcc_sizes_df <- tibble(threshold = 1-2^(-seq(1:20)), power = seq(1:20), lcc_sizes = lcc_sizes)

p <- ggplot(lcc_sizes_df) +
  geom_line(aes(x = power, y = lcc_sizes)) +
  labs(x = "Threshold (1 - 2^(-x))", y = "LCC Size", title = "Thresholding the Edge Weights") +
  theme_light()

ggsave(p, filename = here("fig", "threshold.pdf"), width = 6, height = 3, dpi = "retina")
################################################################################
edges_select <- edges %>%
  mutate(weight = as.integer(weight > 1-2^(-10))) %>%
  filter(weight == 1)

graph_select <- graph_from_data_frame(d = edges_select, vertices = nodes, directed = FALSE)

graph_select_comps <- components(graph_select, mode = "weak")
max_lcc <- which.max(graph_select_comps$csize)

lcc_ids <- V(graph_select)[graph_select_comps$membership == max_lcc]

# subgraph
graph_select_lcc <- igraph::induced_subgraph(graph_select, lcc_ids)

gsize(graph_select)
gorder(graph_select)
gsize(graph_select_lcc)
gorder(graph_select_lcc)

A_lcc <- as_adjacency_matrix(graph_select_lcc, type = "both", sparse = FALSE)

saveRDS(A_lcc, file = here("data", "A_lcc.rds"))
################################################################################
################################################################################
################################################################################

