# Optional convenience runner. Execute from the repository root.

analyses <- c(
  "01_breed_stratified_DE",
  "02_basal_control_DE",
  "03_interaction_model",
  "04_performance_stratified_DE"
)

scripts <- c(
  "scripts/02_run_breed_deseq2.R",
  "scripts/03_run_basal_control_DE.R",
  "scripts/01_run_B_three_way_interaction_only.R",
  "scripts/03_run_performance_stratified_condition_DE.R"
)

root <- getwd()
for (i in seq_along(analyses)) {
  setwd(file.path(root, analyses[i]))
  source(scripts[i])
  setwd(root)
}
