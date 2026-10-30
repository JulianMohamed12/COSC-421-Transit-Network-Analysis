library(tidytransit)
library(igraph)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(hms)

# ==========================================
# STEP 0: Set Working Directory Dynamically
# ==========================================

# Ensures paths resolve relative to this script's directory (code/)
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# Path to the zip file inside the 'code' folder
zip_path <- "transit_data.zip"

# ==========================================
# STEP 1: Load and Clean GTFS Data
# ==========================================

cat("Loading GTFS zip archive from:", zip_path, "\n")
gtfs <- read_gtfs(zip_path)

# Process stop sequences and convert time strings to seconds
stop_times <- gtfs$stop_times %>%
  select(trip_id, arrival_time, departure_time, stop_id, stop_sequence) %>%
  mutate(
    dep_sec = as.numeric(as_hms(departure_time)),
    arr_sec = as.numeric(as_hms(arrival_time))
  ) %>%
  arrange(trip_id, stop_sequence)

# Join trips with routes to classify transit modes
trips_routes <- gtfs$trips %>%
  select(trip_id, route_id) %>%
  inner_join(gtfs$routes %>% select(route_id, route_type, route_short_name), by = "route_id")

# Create directed pairs of sequential stops (Edges)
edges_raw <- stop_times %>%
  group_by(trip_id) %>%
  mutate(
    next_stop_id = lead(stop_id),
    next_arr_sec = lead(arr_sec)
  ) %>%
  filter(!is.na(next_stop_id)) %>%
  ungroup() %>%
  inner_join(trips_routes, by = "trip_id")

# Calculate scheduled segment travel times
edges_raw <- edges_raw %>%
  mutate(
    travel_time_sec = next_arr_sec - dep_sec
  ) %>%
  filter(travel_time_sec > 0) # Remove schedule anomalies

# Aggregate network edges and assign mode groupings
edges <- edges_raw %>%
  group_by(from_stop_id = stop_id, to_stop_id = next_stop_id) %>%
  summarise(
    weight = mean(travel_time_sec, na.rm = TRUE),
    route_count = n_distinct(route_id),
    # GTFS modes: 1 = SkyTrain/Subway, 2 = Rail (WCE), 3 = Bus, 4 = SeaBus/Ferry
    mode = ifelse(any(route_type %in% c(1, 2, 4)), "Train/Ferry", "Bus"),
    .groups = "drop"
  )

# Extract unique Network Nodes
nodes <- gtfs$stops %>%
  filter(stop_id %in% c(edges$from_stop_id, edges$to_stop_id)) %>%
  select(stop_id, stop_name, stop_lat, stop_lon)

# Construct full directed, weighted network graph
g_full <- graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)


# ==========================================
# RESEARCH QUESTION 1:
# Bottlenecks and Points of Failure
# ==========================================

V(g_full)$degree_in <- degree(g_full, mode = "in")
V(g_full)$degree_out <- degree(g_full, mode = "out")
V(g_full)$degree_total <- degree(g_full, mode = "all")
V(g_full)$betweenness <- betweenness(g_full, directed = TRUE, weights = E(g_full)$weight)
V(g_full)$closeness <- closeness(g_full, mode = "all", weights = E(g_full)$weight)

metrics_df <- data.frame(
  stop_id = V(g_full)$name,
  stop_name = V(g_full)$stop_name,
  degree = V(g_full)$degree_total,
  betweenness = V(g_full)$betweenness,
  closeness = V(g_full)$closeness
)

top_bottlenecks <- metrics_df %>%
  arrange(desc(betweenness)) %>%
  head(10)

cat("\n--- TOP 10 BOTTLENECKS (BETWEENNESS CENTRALITY) ---\n")
print(top_bottlenecks)


# ==========================================
# RESEARCH QUESTION 2:
# Disruption & Network Removal Simulations
# ==========================================

simulate_removal <- function(graph, top_n = 20) {
  targets <- head(V(graph)$name[order(betweenness(graph, weights = E(graph)$weight), decreasing = TRUE)], top_n)
  
  results <- data.frame(
    nodes_removed = 0,
    avg_path_length = mean_distance(graph, directed = TRUE, unconnected = TRUE),
    num_components = count_components(graph, mode = "weak")
  )
  
  g_temp <- graph
  for (i in 1:top_n) {
    node_to_remove <- targets[i]
    g_temp <- delete_vertices(g_temp, node_to_remove)
    
    results <- rbind(results, data.frame(
      nodes_removed = i,
      avg_path_length = mean_distance(g_temp, directed = TRUE, unconnected = TRUE),
      num_components = count_components(g_temp, mode = "weak")
    ))
  }
  return(results)
}

removal_results <- simulate_removal(g_full, top_n = 20)

# Plot network robustness decay
p <- ggplot(removal_results, aes(x = nodes_removed, y = avg_path_length)) +
  geom_line(color = "#1f78b4", linewidth = 1) +
  geom_point(color = "#e31a1c", size = 2) +
  theme_minimal() +
  labs(
    title = "Network Vulnerability: Impact of Sequential Hub Removal",
    x = "Number of Top Bottlenecks Removed",
    y = "Average Shortest Path Length"
  )

print(p)


# ==========================================
# RESEARCH QUESTION 3:
# Hidden Bottlenecks (High Routes, Low Accessibility)
# ==========================================

metrics_df <- metrics_df %>%
  mutate(
    z_degree = scale(degree),
    z_closeness = scale(closeness),
    inefficiency_score = z_degree - z_closeness
  )

hidden_bottlenecks <- metrics_df %>%
  arrange(desc(inefficiency_score)) %>%
  head(10)

cat("\n--- TOP 10 HIDDEN BOTTLENECKS (HIGH CONNECTIVITY, POOR ACCESSIBILITY) ---\n")
print(hidden_bottlenecks %>% select(stop_id, stop_name, degree, closeness, inefficiency_score))


# ==========================================
# RESEARCH QUESTION 4:
# Structural Comparison: Bus vs. Rapid Transit
# ==========================================

# Subgraph: Bus
edges_bus <- edges %>% filter(mode == "Bus")
g_bus <- graph_from_data_frame(
  d = edges_bus, 
  vertices = nodes %>% filter(stop_id %in% c(edges_bus$from_stop_id, edges_bus$to_stop_id)), 
  directed = TRUE
)

# Subgraph: Rail & Ferry
edges_train <- edges %>% filter(mode == "Train/Ferry")
g_train <- graph_from_data_frame(
  d = edges_train, 
  vertices = nodes %>% filter(stop_id %in% c(edges_train$from_stop_id, edges_train$to_stop_id)), 
  directed = TRUE
)

# Modal comparison summary
mode_comparison <- data.frame(
  Metric = c("Node Count", "Edge Count", "Network Density", "Average Path Length", "Global Clustering Coeff"),
  Bus_Network = c(
    vcount(g_bus),
    ecount(g_bus),
    edge_density(g_bus),
    mean_distance(g_bus, directed = TRUE, unconnected = TRUE),
    transitivity(g_bus, type = "global")
  ),
  Train_Ferry_Network = c(
    vcount(g_train),
    ecount(g_train),
    edge_density(g_train),
    mean_distance(g_train, directed = TRUE, unconnected = TRUE),
    transitivity(g_train, type = "global")
  )
)

cat("\n--- MODAL NETWORK COMPARISON ---\n")
print(mode_comparison)