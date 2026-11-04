# Code

All analysis is in R (igraph, data.table, ggplot2).

Run everything from the repo root:

```
Rscript code/run_all.R
```

| Script | What it does |
|---|---|
| `00_setup.R` | Packages, paths, settings (service date, walking radius) and helper functions |
| `01_build_network.R` | Reads `transit_data.zip` (TransLink GTFS), keeps one weekday of service, builds the directed travel-time weighted graph |
| `02_centrality.R` | Network summary, degree / betweenness / harmonic closeness, RQ1 (bottlenecks) and RQ3 (hidden bottlenecks) |
| `03_resilience.R` | RQ2: targeted vs random station removal, single-station shutdowns |
| `04_mode_comparison.R` | RQ4: bus vs rapid transit structure, mode removal, bus substitution ratio for each rail link |

Outputs go to `output/tables/` (CSV) and `output/figures/` (PNG). The whole run takes about a minute.

`main.R` is the earlier single-file draft and needs the `tidytransit` package; the numbered scripts replace it.
