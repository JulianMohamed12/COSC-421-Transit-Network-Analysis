# 00_setup.R
# Shared packages, paths and helper functions.
# Every other script sources this file first.

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
  library(ggplot2)
})

# Locate the repository root
# Works whether the script is run from the repo root, from code/, or in RStudio.
find_repo_root <- function() {
  start <- getwd()
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(p)) start <- dirname(p)
  }
  d <- normalizePath(start)
  for (i in 1:5) {
    if (file.exists(file.path(d, "code", "transit_data.zip"))) return(d)
    d <- dirname(d)
  }
  stop("Could not find the repo root (looking for code/transit_data.zip). ",
       "Run this from inside the COSC-421-Transit-Network-Analysis folder.")
}

ROOT     <- find_repo_root()
ZIP_PATH <- file.path(ROOT, "code", "transit_data.zip")
OUT_DIR  <- file.path(ROOT, "output")
FIG_DIR  <- file.path(OUT_DIR, "figures")
TAB_DIR  <- file.path(OUT_DIR, "tables")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

GRAPH_RDS <- file.path(OUT_DIR, "transit_graph.rds")

# Analysis settings
# A normal weekday inside the feed's calendar (Sept 7 2026 - Jan 3 2027),
# chosen to avoid holidays. Using one day avoids mixing weekday and weekend
# schedules into the same edge weights.
SERVICE_DATE <- as.Date("2026-11-18")

MIN_EDGE_MIN          <- 0.25   # floor on edge travel time (15 s); GTFS rounds to 30 s so some hops show 0
DEFAULT_TRANSFER_MIN  <- 2      # walking time used when transfers.txt has no min_transfer_time
WALK_RADIUS_M         <- 400    # extra walking links between rapid transit stations and nearby stops
WALK_SPEED_M_PER_MIN  <- 80     # about 4.8 km/h
N_EFFICIENCY_SOURCES  <- 300    # sampled source stops for the efficiency estimate
set.seed(421)

# Helpers

# GTFS times look like " 5:05:00" or "25:10:30" (after midnight). Convert to seconds.
gtfs_time_to_sec <- function(x) {
  x <- trimws(x)
  p <- tstrsplit(x, ":", fixed = TRUE)
  as.integer(p[[1]]) * 3600L + as.integer(p[[2]]) * 60L + as.integer(p[[3]])
}

# GTFS route_type codes used in the TransLink feed
mode_from_route_type <- function(rt) {
  fcase(rt == 1L, "SkyTrain",
        rt == 2L, "West Coast Express",
        rt == 3L, "Bus",
        rt == 4L, "SeaBus",
        default = "Other")
}

# Network efficiency (Latora & Marchiori), estimated from a fixed sample of
# source stops so it stays fast on ~9,000 nodes. Uses travel time (minutes).
# Normalised by the ORIGINAL node count so removing nodes lowers the score.
sampled_efficiency <- function(g, sources, n_original) {
  src <- intersect(sources, V(g)$name)
  if (length(src) == 0) return(0)
  D <- distances(g, v = src, mode = "out", weights = E(g)$weight)
  inv <- 1 / D
  inv[!is.finite(inv)] <- 0          # unreachable pairs and self-distance (1/0)
  sum(inv) / (length(sources) * (n_original - 1))
}

# Share of the original nodes still in the largest connected piece
largest_component_share <- function(g, n_original, mode = "weak") {
  if (vcount(g) == 0) return(0)
  max(components(g, mode = mode)$csize) / n_original
}

# Great-circle distance in metres
haversine_m <- function(lat1, lon1, lat2, lon2) {
  r <- pi / 180
  a <- sin((lat2 - lat1) * r / 2)^2 + cos(lat1 * r) * cos(lat2 * r) * sin((lon2 - lon1) * r / 2)^2
  2 * 6371000 * asin(pmin(1, sqrt(a)))
}

save_table <- function(df, name) {
  fwrite(df, file.path(TAB_DIR, paste0(name, ".csv")))
  invisible(df)
}

save_plot <- function(p, name, width = 8, height = 6) {
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), p, width = width, height = height, dpi = 200, bg = "white")
  invisible(p)
}

load_graph <- function() {
  if (!file.exists(GRAPH_RDS)) {
    message("Graph not built yet - running 01_build_network.R")
    source(file.path(ROOT, "code", "01_build_network.R"))
  }
  readRDS(GRAPH_RDS)
}
