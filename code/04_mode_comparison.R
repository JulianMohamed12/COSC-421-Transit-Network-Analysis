# 04_mode_comparison.R
#   RQ4: Do buses and rapid transit (SkyTrain, WCE, SeaBus) play
#        different roles, or do they duplicate each other?
#
# Part A - structure of each mode's own subnetwork
# Part B - what happens to the integrated network when a mode
#          (or the transfers between modes) is removed
# Part C - for every rapid transit link, how long the same trip
#          takes using only buses + walking (substitution ratio)

if (!exists("ROOT")) source(if (file.exists("00_setup.R")) "00_setup.R" else file.path("code", "00_setup.R"))
g <- load_graph()
N0 <- vcount(g)
SOURCES <- sample(V(g)$name, N_EFFICIENCY_SOURCES)
E0 <- sampled_efficiency(g, SOURCES, N0)

# Part A: subnetwork structure
mode_subgraph <- function(groups) {
  subgraph_from_edges(g, E(g)[mode_group %in% groups], delete.vertices = TRUE)
}
describe <- function(gr, label) {
  data.table(network = label, nodes = vcount(gr), edges = ecount(gr),
             density = edge_density(gr),
             avg_degree = mean(degree(gr, mode = "all")),
             avg_path_hops = mean_distance(gr, weights = NA, directed = TRUE),
             avg_travel_min = mean_distance(gr, weights = E(gr)$weight, directed = TRUE),
             avg_edge_min = mean(E(gr)$weight),
             clustering = transitivity(as_undirected(gr, mode = "collapse"), type = "global"),
             weak_components = components(gr, mode = "weak")$no)
}
g_bus   <- mode_subgraph("Bus")
g_rapid <- mode_subgraph("Rapid transit")
structure_tab <- rbind(describe(g_bus, "Bus only"), describe(g_rapid, "Rapid transit only"),
                       describe(g, "Integrated (all modes + transfers)"))
cat("\n--- RQ4A: SUBNETWORK STRUCTURE ---\n"); print(structure_tab)
save_table(structure_tab, "rq4_mode_structure")

# Part B: remove a whole mode
scenario <- function(drop_groups, label) {
  gr <- delete_edges(g, E(g)[mode_group %in% drop_groups])
  data.table(scenario = label,
             efficiency_pct = round(100 * sampled_efficiency(gr, SOURCES, N0) / E0, 1),
             largest_component_pct = round(100 * largest_component_share(gr, N0, "weak"), 1))
}
mode_removal <- rbind(
  scenario(character(0),    "Integrated network"),
  scenario("Rapid transit", "No rapid transit"),
  scenario("Transfer",      "No walking transfers (modes not integrated)"),
  scenario("Bus",           "No buses")
)
cat("\n--- RQ4B: MODE REMOVAL ---\n"); print(mode_removal)
save_table(mode_removal, "rq4_mode_removal")

# Part C: can buses substitute for each rapid transit link?
g_no_rapid <- delete_edges(g, E(g)[mode_group == "Rapid transit"])
rapid_el <- as.data.table(as_data_frame(g, what = "edges"))[mode_group == "Rapid transit"]
from_ids <- unique(rapid_el$from)
D <- distances(g_no_rapid, v = from_ids, to = unique(rapid_el$to), mode = "out",
               weights = E(g_no_rapid)$weight)
rapid_el[, bus_min := D[cbind(from, to)]]
rapid_el[, substitution_ratio := bus_min / weight]
rapid_el[, from_name := V(g)[from]$stop_name]
rapid_el[, to_name   := V(g)[to]$stop_name]

subst <- rapid_el[, .(mode, routes, from_name, to_name, rail_min = round(weight, 1),
                      bus_min = round(bus_min, 1), substitution_ratio = round(substitution_ratio, 1))]
setorder(subst, -substitution_ratio)
cat("\n--- RQ4C: RAPID TRANSIT LINKS WITH THE WORST BUS ALTERNATIVE ---\n"); print(head(subst, 15))
cat(sprintf("\nMedian substitution ratio: %.1fx | links with no bus alternative: %d of %d\n",
            median(subst$substitution_ratio[is.finite(subst$substitution_ratio)]),
            sum(!is.finite(subst$substitution_ratio)), nrow(subst)))
save_table(subst, "rq4_substitution_ratio")

p_sub <- ggplot(subst[is.finite(substitution_ratio)], aes(substitution_ratio, fill = mode)) +
  geom_histogram(bins = 30, colour = "white") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  theme_minimal() +
  labs(title = "RQ4: How much slower is the best bus/walk alternative to each rapid transit link?",
       subtitle = "Ratio = fastest bus + walk time / scheduled rail or ferry time (log scale; 1 = same speed)",
       x = "Substitution ratio", y = "Number of rapid transit links", fill = NULL)
save_plot(p_sub, "rq4_substitution_ratio", 9, 5)
