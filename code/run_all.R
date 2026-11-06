# run_all.R - runs the full milestone pipeline in order.
# From the repo root:  Rscript code/run_all.R
# Outputs: output/tables/*.csv, output/figures/*.png

source(if (file.exists("00_setup.R")) "00_setup.R" else file.path("code", "00_setup.R"))

for (s in c("01_build_network.R", "02_centrality.R", "03_resilience.R", "04_mode_comparison.R")) {
  cat("\n==========", s, "==========\n")
  set.seed(421)
  source(file.path(ROOT, "code", s), echo = FALSE)
}
cat("\nAll done. Tables in", TAB_DIR, "| figures in", FIG_DIR, "\n")
