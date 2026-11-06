# 01_build_network.R
# Load the TransLink GTFS feed, clean it, and build the
# directed, weighted, multimodal transit graph.
#
# Nodes : stops. SkyTrain / WCE platforms are merged into
#         their parent station so a station is one node.
# Edges : a direct ride between two consecutive stops on any
#         trip running on SERVICE_DATE. Weight = median
#         scheduled travel time (minutes). Edges are tagged by
#         mode. Walking transfers from transfers.txt are added
#         as "Transfer" edges so bus and rail connect.

if (!exists("ROOT")) source(if (file.exists("00_setup.R")) "00_setup.R" else file.path("code", "00_setup.R"))

# 1. Read GTFS
gtfs_dir <- file.path(tempdir(), "gtfs_translink")
dir.create(gtfs_dir, showWarnings = FALSE)
unzip(ZIP_PATH, exdir = gtfs_dir)
read_gtfs_file <- function(f) fread(file.path(gtfs_dir, f), colClasses = "character", na.strings = "")

stops      <- read_gtfs_file("stops.txt")
routes     <- read_gtfs_file("routes.txt")
trips      <- read_gtfs_file("trips.txt")
stop_times <- read_gtfs_file("stop_times.txt")
calendar   <- read_gtfs_file("calendar.txt")
cal_dates  <- read_gtfs_file("calendar_dates.txt")
transfers  <- read_gtfs_file("transfers.txt")

cat(sprintf("Raw feed: %d stops, %d routes, %d trips, %d stop_times\n",
            nrow(stops), nrow(routes), nrow(trips), nrow(stop_times)))

# 2. Keep only service running on SERVICE_DATE
day_col  <- tolower(weekdays(SERVICE_DATE))
date_str <- format(SERVICE_DATE, "%Y%m%d")

active <- calendar[get(day_col) == "1" & start_date <= date_str & end_date >= date_str, service_id]
added   <- cal_dates[date == date_str & exception_type == "1", service_id]
removed <- cal_dates[date == date_str & exception_type == "2", service_id]
active  <- setdiff(union(active, added), removed)

routes[, route_type := as.integer(route_type)]
routes[, mode := mode_from_route_type(route_type)]

# Drop HandyDART (route_type 715): it is on-demand, not a fixed route network
trips <- trips[service_id %in% active]
trips <- merge(trips[, .(trip_id, route_id)],
               routes[mode != "Other", .(route_id, mode, route_short_name, route_long_name)],
               by = "route_id")
trips[, route_label := fifelse(is.na(route_short_name), route_long_name, route_short_name)]

cat(sprintf("Service on %s (%s): %d trips on %d routes\n",
            SERVICE_DATE, day_col, nrow(trips), uniqueN(trips$route_id)))

# 3. Map platforms to stations (node ids)
stops[, node_id := fifelse(!is.na(parent_station), parent_station, stop_id)]

# West Coast Express platforms ("X Station @ West Coast Express") have no parent
# in the feed. Attach them to the SkyTrain station of the same name if one exists.
stations <- stops[location_type == "1", .(station_name = stop_name, station_id = stop_id)]
wce <- stops[is.na(parent_station) & grepl(" @ West Coast Express$", stop_name)]
wce[, station_name := sub(" @ West Coast Express$", "", stop_name)]
wce <- merge(wce, stations, by = "station_name")
stops[wce, node_id := i.station_id, on = "stop_id"]
cat(sprintf("Merged %d West Coast Express platforms into SkyTrain stations\n", nrow(wce)))
stop_to_node <- stops[, .(stop_id, node_id)]

# 4. Build ride edges from consecutive stop_times
st <- stop_times[trip_id %in% trips$trip_id,
                 .(trip_id, stop_id, stop_sequence = as.integer(stop_sequence),
                   arr = gtfs_time_to_sec(arrival_time),
                   dep = gtfs_time_to_sec(departure_time))]
setorder(st, trip_id, stop_sequence)

st[, `:=`(next_trip = shift(trip_id, -1L),
          next_stop = shift(stop_id, -1L),
          next_arr  = shift(arr, -1L))]
hops <- st[next_trip == trip_id, .(trip_id, from_stop = stop_id, to_stop = next_stop,
                                   tt_min = (next_arr - dep) / 60)]
rm(st, stop_times); invisible(gc())

hops <- merge(hops, trips[, .(trip_id, route_id, mode, route_label)], by = "trip_id")
hops[, from := stop_to_node[.(hops$from_stop), on = "stop_id", node_id]]
hops[, to   := stop_to_node[.(hops$to_stop),   on = "stop_id", node_id]]
hops <- hops[from != to]                     # platform-to-platform inside one station
n_negative <- hops[tt_min < 0, .N]
hops <- hops[tt_min >= 0]

ride_edges <- hops[, .(
  weight        = max(median(tt_min), MIN_EDGE_MIN),
  trips_per_day = .N,
  n_routes      = uniqueN(route_id),
  routes        = paste(sort(unique(route_label)), collapse = ";"),
  mode          = fifelse(any(mode != "Bus"), mode[mode != "Bus"][1], "Bus")
), by = .(from, to)]

cat(sprintf("Ride edges: %d (dropped %d hops with negative times)\n", nrow(ride_edges), n_negative))

# 5. Walking transfer edges
# transfer_type 4 rows are trip-specific (in-seat) transfers, not walks between stops.
walk <- transfers[is.na(from_trip_id) & from_stop_id != to_stop_id & transfer_type != "3",
                  .(from_stop = from_stop_id, to_stop = to_stop_id,
                    w = as.numeric(min_transfer_time) / 60)]
walk[is.na(w), w := DEFAULT_TRANSFER_MIN]
walk[, from := stop_to_node[.(walk$from_stop), on = "stop_id", node_id]]
walk[, to   := stop_to_node[.(walk$to_stop),   on = "stop_id", node_id]]
walk <- walk[from != to & !is.na(from) & !is.na(to)]
transfer_edges <- walk[, .(weight = max(min(w), MIN_EDGE_MIN), trips_per_day = NA_integer_,
                           n_routes = 0L, routes = "", mode = "Transfer"), by = .(from, to)]
# Only keep a transfer if no ride edge already joins the same pair
transfer_edges <- transfer_edges[!ride_edges, on = .(from, to)]

# transfers.txt misses some obvious walks (e.g. Templeton Station to the bus
# stop outside it). Add walking links from each rapid transit station to any
# stop within WALK_RADIUS_M.
coords <- stops[, .(node_id = stop_id, lat = as.numeric(stop_lat), lon = as.numeric(stop_lon))]
rapid_nodes <- unique(hops[mode != "Bus", c(from, to)])
all_nodes   <- coords[node_id %in% unique(c(hops$from, hops$to))]
near <- rbindlist(lapply(rapid_nodes, function(r) {
  rc <- coords[node_id == r]
  d  <- haversine_m(rc$lat, rc$lon, all_nodes$lat, all_nodes$lon)
  keep <- d <= WALK_RADIUS_M & all_nodes$node_id != r
  data.table(from = r, to = all_nodes$node_id[keep], weight = d[keep] / WALK_SPEED_M_PER_MIN)
}))
near <- rbind(near, near[, .(from = to, to = from, weight)])
near[, weight := pmax(weight, MIN_EDGE_MIN)]
near <- unique(near, by = c("from", "to"))
near <- near[!ride_edges, on = .(from, to)][!transfer_edges, on = .(from, to)]
proximity_edges <- near[, .(from, to, weight, trips_per_day = NA_integer_, n_routes = 0L,
                            routes = "", mode = "Transfer")]
transfer_edges <- rbind(transfer_edges, proximity_edges)

edges <- rbind(ride_edges, transfer_edges)
cat(sprintf("Transfer edges: %d (incl. %d added by %d m proximity rule)\n",
            nrow(transfer_edges), nrow(proximity_edges), WALK_RADIUS_M))

# 6. Node table
used_nodes <- unique(c(edges$from, edges$to))
node_info <- stops[stop_id %in% used_nodes,
                   .(name = stop_id, stop_name,
                     lat = as.numeric(stop_lat), lon = as.numeric(stop_lon))]

# Number of distinct routes serving each node, and whether rapid transit stops there
served <- rbind(hops[, .(node = from, route_id, mode)], hops[, .(node = to, route_id, mode)])
node_routes <- served[, .(routes_served = uniqueN(route_id),
                          node_mode = fifelse(any(mode != "Bus"), "Rapid transit", "Bus")),
                      by = .(name = node)]
node_info <- merge(node_info, node_routes, by = "name", all.x = TRUE)
node_info[is.na(routes_served), `:=`(routes_served = 0L, node_mode = "Transfer only")]
setcolorder(node_info, "name")

# 7. Graph
edges[, mode_group := fifelse(mode == "Bus", "Bus",
                       fifelse(mode == "Transfer", "Transfer", "Rapid transit"))]
g <- graph_from_data_frame(edges, directed = TRUE, vertices = node_info)

cat(sprintf("\nGraph built: %d nodes, %d edges (%s)\n", vcount(g), ecount(g),
            paste(names(table(E(g)$mode)), table(E(g)$mode), sep = "=", collapse = ", ")))

saveRDS(g, GRAPH_RDS)
save_table(node_info, "nodes")
save_table(edges, "edges")
