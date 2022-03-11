# A_network <- network(A_edges_extended, directed = FALSE)
# 
# # A_network <- network(A_edges_extended, vertex.attr = A_nodes, matrix.type = "edgelist", ignore.eval = FALSE)
# A_adj <- A_edges_extended %>%
#   dplyr::select(from, to, weight) %>%
#   pivot_wider(names_from = to, values_from = weight)