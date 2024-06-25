#!/usr/bin/env -S Rscript --vanilla

utils <- new.env()
sys.source("cmdstan_utils.R", envir = utils)
utils$load_packages(plot_libs = TRUE)

fit <- readRDS(file.path(
  "Rout", "pop_fit_free_weights_cmdstan_max_delay_9_job3.rds"
))

parameters <- c(
  "sigma_CBC", "sigma_haemocue", "sigma_retic", "CBC_correction",
  "logit_alpha", "beta", "h", "log_k", "Hb_star", "T_E_star",
  "alpha_diff1", "alpha_delta1", "alpha_diff2", "alpha_delta2"
)
draws <- utils$get_fit_draws_long(fit, parameters)

inv_logit <- function(x) {
  exp(x)/(1+exp(x))
}

draws_alpha <- draws |>
  filter(name == "logit_alpha") |>
  mutate(
    value = 100 * inv_logit(value),
    name = "alpha"
  )

all_draws <- bind_rows(draws, draws_alpha)

draw_summary <- all_draws |>
  group_by(name) |>
  summarise(
    mean = mean(value),
    lower = quantile(value, probs = 0.05),
    upper = quantile(value, probs = 0.95),
    .groups = "drop"
  )

labels <- c(
  "alpha" = "Max RBC lifespan reduction (%)",
  "beta" = "Half-maximal effect dose (mg/kg)",
  "h" = "Dose-response slope",
  "Hb_star" = "Steady-state haemoglobin (g/dL)",
  "T_E_star" = "Steady-state RBC lifespan (days)"
)

df_table <- draw_summary |>
  filter(name %in% names(labels)) |>
  mutate(name = labels[name])

print(
  df_table |>
    knitr::kable(
      format = "latex",
      booktabs = TRUE,
      digits = 2,
      col.names = c("Parameter", "Mean", "5%", "95%")
    )
)
