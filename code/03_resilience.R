# ==========================================================
# 03_resilience.R
#   RQ2: How well does the network hold up when key stations are removed?
#
# Three removal strategies are compared:
#   - targeted by betweenness (bottlenecks from RQ1)
#   - targeted by routes served (the "busiest-looking" stops)
#   - random removal (baseline, averaged over several runs)
# After each step we measure:
#   - share of stops still in the largest weakly connected component
#   - network efficiency (mean of 1 / travel time), from a fixed sample of sources
# ==========================================================

if (!exists("ROOT")) source(if (file.exists("00_setup.R")) "00_setup.R" else file.path("code", "00_setup.R"))
g <- load_graph()
if (is.null(V(g)$betweenness)) {
  source(file.path(ROOT, "code", "02_centrality.R")); g <- readRDS(GRAPH_RDS)
}

N0       <- vcount(g)
SOURCES  <- sample(V(g)$name, N_EFFICIENCY_SOURCES)
MAX_REMOVE <- 50
STEPS    <- c(0, 1, 2, 3, 5, 10, 15, 20, 30, 40, 50)
N_RANDOM <- 5

E0 <- sampled_efficiency(g, SOURCES, N0)

evaluate <- function(gr, k, strategy, run = 1L) {
  data.table(strategy = strategy, run = run, removed = k,
             lcc_share  = largest_component_share(gr, N0, "weak"),
             efficiency = sampled_efficiency(gr, SOURCES, N0) / E0)
}

run_removal <- function(order_names, strategy, run = 1L) {
  rbindlist(lapply(STEPS, function(k) {
    gr <- if (k == 0) g else delete_vertices(g, order_names[seq_len(k)])
    evaluate(gr, k, strategy, run)
  }))
}

by_betweenness <- V(g)$name[order(V(g)$betweenness, decreasing = TRUE)]
by_routes      <- V(g)$name[order(V(g)$routes_served, V(g)$betweenness, decreasing = TRUE)]

cat("Running removal simulations...\n")
res <- rbind(
  run_removal(by_betweenness, "Targeted: betweenness"),
  run_removal(by_routes,      "Targeted: routes served"),
  rbindlist(lapply(seq_len(N_RANDOM), function(r) run_removal(sample(V(g)$name, MAX_REMOVE), "Random", r)))
)
summary_res <- res[, .(lcc_share = mean(lcc_share), efficiency = mean(efficiency)), by = .(strategy, removed)]
cat("\n--- RQ2: EFFICIENCY / LARGEST COMPONENT AFTER REMOVAL ---\n"); print(dcast(summary_res, removed ~ strategy, value.var = "efficiency"))
save_table(summary_res, "rq2_removal_results")

# Which single station hurts the most on its own? (top 25 by betweenness)
base_lcc <- max(components(g, mode = "weak")$csize)
single <- rbindlist(lapply(by_betweenness[1:25], function(v) {
  gr <- delete_vertices(g, v)
  data.table(stop_name = V(g)[v]$stop_name, node_mode = V(g)[v]$node_mode,
             efficiency_loss_pct = round(100 * (1 - sampled_efficiency(gr, SOURCES, N0) / E0), 2),
             stops_cut_off = base_lcc - 1 - max(components(gr, mode = "weak")$csize))
}))
setorder(single, -efficiency_loss_pct)
cat("\n--- RQ2: SINGLE-STATION SHUTDOWNS (top 25 by betweenness) ---\n"); print(single)
save_table(single, "rq2_single_station_removal")

p_eff <- ggplot(summary_res, aes(removed, efficiency, colour = strategy)) +
  geom_line(linewidth = 1) + geom_point(size = 1.8) +
  scale_colour_manual(values = c("Targeted: betweenness" = "#d53f3f",
                                 "Targeted: routes served" = "#dd8a1f", "Random" = "#2b6cb0")) +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  labs(title = "RQ2: Network efficiency as stations are removed",
       x = "Stations removed", y = "Efficiency (% of intact network)", colour = NULL)
save_plot(p_eff, "rq2_efficiency_decay", 8, 5)

p_lcc <- ggplot(summary_res, aes(removed, lcc_share, colour = strategy)) +
  geom_line(linewidth = 1) + geom_point(size = 1.8) +
  scale_colour_manual(values = c("Targeted: betweenness" = "#d53f3f",
                                 "Targeted: routes served" = "#dd8a1f", "Random" = "#2b6cb0")) +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  labs(title = "RQ2: Largest connected component as stations are removed",
       x = "Stations removed", y = "Share of stops still connected", colour = NULL)
save_plot(p_lcc, "rq2_largest_component", 8, 5)
