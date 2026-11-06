# 02_centrality.R
# Network-level summary + node centrality.
#   RQ1: Which stations are the biggest bottlenecks / points of failure?
#   RQ3: Which stations have many routes but are poor for getting across town?

if (!exists("ROOT")) source(if (file.exists("00_setup.R")) "00_setup.R" else file.path("code", "00_setup.R"))
g <- load_graph()

# Network-level summary
ride_g <- subgraph_from_edges(g, E(g)[mode != "Transfer"], delete.vertices = FALSE)
deg_all <- degree(g, mode = "all")
scc <- components(g, mode = "strong")
wcc <- components(g, mode = "weak")

network_summary <- data.table(
  metric = c("Nodes", "Edges (all)", "Edges (rides only)", "Density",
             "Average in-degree (= out-degree)", "Average total degree",
             "Max total degree", "Weak components", "Largest weak component (% nodes)",
             "Strong components", "Largest strong component (% nodes)",
             "Average shortest path (hops)", "Average shortest path (minutes)",
             "Global clustering coefficient"),
  value = c(vcount(g), ecount(g), ecount(ride_g), edge_density(g),
            mean(degree(g, mode = "in")), mean(deg_all), max(deg_all),
            wcc$no, 100 * max(wcc$csize) / vcount(g),
            scc$no, 100 * max(scc$csize) / vcount(g),
            mean_distance(g, weights = NA, directed = TRUE, unconnected = TRUE),  # NA = ignore weights (hops)
            mean_distance(g, weights = E(g)$weight, directed = TRUE, unconnected = TRUE),
            transitivity(as_undirected(g, mode = "collapse"), type = "global"))
)
cat("\n--- NETWORK SUMMARY ---\n"); print(network_summary)
save_table(network_summary, "network_summary")

# Node centrality
cat("\nComputing betweenness and harmonic closeness (weighted by travel time)...\n")
t0 <- Sys.time()
V(g)$degree_in    <- degree(g, mode = "in")
V(g)$degree_out   <- degree(g, mode = "out")
V(g)$degree       <- deg_all
V(g)$ride_degree  <- degree(ride_g, mode = "all")
V(g)$betweenness  <- betweenness(g, directed = TRUE, weights = E(g)$weight, normalized = TRUE)
# Harmonic closeness instead of classic closeness: the directed graph is not
# strongly connected, and classic closeness is undefined/misleading there.
V(g)$closeness    <- harmonic_centrality(g, mode = "out", weights = E(g)$weight, normalized = TRUE)
cat(sprintf("Done in %.1f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# Articulation points (on the undirected version): stops whose loss splits the network
art <- articulation_points(as_undirected(g, mode = "collapse"))
V(g)$articulation <- V(g)$name %in% art$name

metrics <- as.data.table(as_data_frame(g, what = "vertices"))
setorder(metrics, -betweenness)
save_table(metrics, "node_metrics")
saveRDS(g, GRAPH_RDS)   # keep the centrality scores on the graph for later scripts

# RQ1: bottlenecks
top_betweenness <- metrics[1:15, .(stop_name, node_mode, routes_served, degree,
                                   betweenness = round(betweenness, 4), articulation)]
cat("\n--- RQ1: TOP 15 BOTTLENECKS (weighted betweenness) ---\n"); print(top_betweenness)
save_table(top_betweenness, "rq1_top_betweenness")

top_degree <- metrics[order(-routes_served)][1:15, .(stop_name, node_mode, routes_served, degree,
                                                     betweenness = round(betweenness, 4))]
cat("\n--- RQ1: TOP 15 STOPS BY ROUTES SERVED ---\n"); print(top_degree)
save_table(top_degree, "rq1_top_routes_served")

cat(sprintf("\nArticulation points (single stops whose loss splits the network): %d of %d\n",
            length(art), vcount(g)))
save_table(metrics[articulation == TRUE][1:min(.N, 25),
                   .(stop_name, node_mode, routes_served, betweenness)], "rq1_top_articulation_points")

# RQ3: hidden bottlenecks
# Degree in a stop-level GTFS graph is almost always 2-4, so "number of routes
# serving the stop" is the better measure of how connected a stop looks.
# A hidden bottleneck is a stop in the top 5% for routes served whose
# closeness is still below the network median. Score = z(routes) - z(closeness).
# Stops with closeness 0 are sinks (end-of-line stops with no outgoing ride);
# they are excluded here and reported separately as a data artefact.
sinks <- metrics[closeness == 0]
cat(sprintf("\nSink stops (no outgoing edges, closeness = 0): %d\n", nrow(sinks)))
save_table(sinks[, .(stop_name, node_mode, routes_served)], "sink_stops")

metrics[, z_routes := as.numeric(scale(routes_served))]
metrics[, z_close  := as.numeric(scale(closeness))]
metrics[, hidden_score := z_routes - z_close]
route_cut <- quantile(metrics$routes_served, 0.95)
close_median <- median(metrics$closeness)

hidden <- metrics[routes_served >= route_cut & closeness > 0 & closeness < close_median][order(-hidden_score)][1:min(.N, 15),
            .(stop_name, node_mode, routes_served, closeness = round(closeness, 4),
              closeness_pct = round(100 * ecdf(metrics$closeness)(closeness)),
              hidden_score = round(hidden_score, 2))]
cat(sprintf("\n--- RQ3: HIDDEN BOTTLENECKS (top 5%% by routes served, i.e. >= %d routes, closeness below median) ---\n", route_cut))
print(hidden)
save_table(hidden, "rq3_hidden_bottlenecks")

rho <- cor(metrics$routes_served, metrics$closeness, method = "spearman")
cat(sprintf("Spearman correlation, routes served vs closeness: %.3f\n", rho))

# Figures
p_deg <- ggplot(metrics, aes(degree)) +
  geom_bar(fill = "#2b6cb0") +
  scale_y_log10() +
  theme_minimal() +
  labs(title = "Degree distribution (all edges)", x = "Total degree", y = "Number of stops (log scale)")
save_plot(p_deg, "degree_distribution", 7, 4.5)

p_rq3 <- ggplot(metrics, aes(routes_served, closeness, colour = node_mode)) +
  geom_jitter(alpha = 0.35, size = 0.8, width = 0.25, height = 0) +
  geom_point(data = metrics[stop_name %in% hidden$stop_name], colour = "black", shape = 21, size = 2.5) +
  scale_colour_manual(values = c("Bus" = "#8aa1b1", "Rapid transit" = "#d53f3f", "Transfer only" = "grey70")) +
  theme_minimal() +
  labs(title = "RQ3: Routes served vs. harmonic closeness",
       subtitle = "Circled = top hidden bottlenecks (many routes, low closeness)",
       x = "Routes serving the stop", y = "Harmonic closeness (travel-time weighted)", colour = NULL)
save_plot(p_rq3, "rq3_routes_vs_closeness")

# Geographic network map, nodes sized by betweenness
el <- as.data.table(as_data_frame(g, what = "edges"))
coords <- metrics[, .(name, lat, lon)]
el <- merge(el, coords[, .(from = name, x = lon, y = lat)], by = "from")
el <- merge(el, coords[, .(to = name, xend = lon, yend = lat)], by = "to")
el <- el[mode != "Transfer"]
el[, mode_group := factor(mode_group, levels = c("Bus", "Rapid transit"))]
setorder(el, mode_group)

p_map <- ggplot() +
  geom_segment(data = el, aes(x, y, xend = xend, yend = yend, colour = mode_group,
                              linewidth = mode_group), lineend = "round") +
  geom_point(data = metrics[1:15], aes(lon, lat, size = betweenness), colour = "#111111", alpha = 0.8) +
  geom_text(data = metrics[1:10], aes(lon, lat, label = stop_name), size = 2.4, vjust = -1, check_overlap = TRUE) +
  scale_colour_manual(values = c("Bus" = "#9fb3c8", "Rapid transit" = "#d53f3f")) +
  scale_linewidth_manual(values = c("Bus" = 0.25, "Rapid transit" = 0.9), guide = "none") +
  coord_quickmap() +
  theme_void() +
  labs(title = "Metro Vancouver transit network (TransLink GTFS)",
       subtitle = "Black points = top 15 stops by travel-time betweenness", colour = NULL, size = "Betweenness")
save_plot(p_map, "network_map_betweenness", 10, 8)
