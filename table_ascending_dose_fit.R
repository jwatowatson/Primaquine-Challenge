#!/usr/bin/env -S Rscript --vanilla

suppressPackageStartupMessages(library(tidyverse))

utils <- new.env()
sys.source("cmdstan_utils.R", envir = utils)
utils$load_packages(plot_libs = TRUE)

fit <- readRDS(file.path(
  "Rout", "pop_fit_free_weights_cmdstan_max_delay_9_job3.rds"
))

parameters <- c("logit_alpha", "beta", "h", "Hb_star", "T_E_star")
draws <- utils$get_fit_draws_long(fit, parameters)

# Sample the individual random effects, noting that `h` has no random effect.
calc <- new.env()
sys.source(
  file.path("optimal-dose-regimens", "calculate-max-daily-Hb-drops.R"),
  envir = calc
)
individual_effects <- calc$sample_individual_effects(fit, seed = 12345) |>
  mutate(.draw = row_number()) |>
  pivot_longer(! .draw) |>
  filter(name %in% parameters) |>
  rename(effect = value)

draws_and_effects <- draws |>
  select(.draw, name, value) |>
  left_join(individual_effects, by = c(".draw", "name")) |>
  # NOTE: effects combine with the base values in parameter-specific ways.
  mutate(
    effect = case_when(
      name == "beta" ~ value * exp(effect),
      name == "h" ~ value,
      name == "Hb_star" ~ value + effect,
      name == "logit_alpha" ~ value + effect,
      name == "T_E_star" ~ value + effect,
    )
  )

inv_logit <- function(x) {
  exp(x)/(1+exp(x))
}

draws_alpha <- draws_and_effects |>
  filter(name == "logit_alpha") |>
  mutate(
    value = 100 * inv_logit(value),
    effect = 100 * inv_logit(effect),
    name = "alpha"
  )

all_draws <- bind_rows(draws_and_effects, draws_alpha) |>
  filter(name != "logit_alpha") |>
  arrange(.draw, name)

popn_intervals_width <- 95
popn_prob_lower <- 0.5 * (1 - popn_intervals_width / 100)
popn_prob_upper <- 1 - popn_prob_lower

pred_intervals_width <- 80
pred_prob_lower <- 0.5 * (1 - pred_intervals_width / 100)
pred_prob_upper <- 1 - pred_prob_lower

# Define the number of decimal places for each parameter.
n_digits <- c("Hb_star" = 1, "T_E_star" = 0, "alpha" = 0, "beta" = 2, "h" = 2)

draw_summary <- all_draws |>
  group_by(name) |>
  summarise(
    mean = mean(value),
    lower = quantile(value, probs = popn_prob_lower),
    upper = quantile(value, probs = popn_prob_upper),
    pred_lower = quantile(effect, probs = pred_prob_lower),
    pred_upper = quantile(effect, probs = pred_prob_upper),
    .groups = "drop"
  ) |>
  mutate(
    format = paste0("%0.", n_digits[name], "f"),
    popn_95 = paste0(
      sprintf(format, mean),
      " (",
      sprintf(format, lower),
      ", ",
      sprintf(format, upper),
      ")"
    ),
    pred_95 = case_when(
      # Do not show a predictive interval for the dose-response slope, because
      # it does not have any inter-individual variation.
      name == "h" ~ "---",
      .default = paste0(
        sprintf(format, pred_lower),
        ", ",
        sprintf(format, pred_upper)
      )
    ),
    format = NULL,
    mean = NULL,
    lower = NULL,
    upper = NULL,
    pred_lower = NULL,
    pred_upper = NULL
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

tbl_latex <- df_table |>
  knitr::kable(
    format = "latex",
    booktabs = TRUE,
    digits = 2,
    col.names = c(
      "Parameter", "Population mean (95% CI)", "Predictive interval"
    )
  )

print(tbl_latex)
cat(tbl_latex, file = "table-ascending-dose-fit.tex")

# Print the 80% CI for RBC lifespan.
tbl_rbc_lifespan_80pcnt_ints <- all_draws |>
  filter(name == "T_E_star") |>
  group_by(name) |>
  summarise(
    popn_lower = round(quantile(value, 0.1), digits = 2),
    popn_upper = round(quantile(value, 0.9), digits = 2),
    indiv_lower = round(quantile(effect, 0.1), digits = 2),
    indiv_upper = round(quantile(effect, 0.9), digits = 2),
    .groups = "drop"
  )

cat("\n")
print(tbl_rbc_lifespan_80pcnt_ints)
write_csv(
  tbl_rbc_lifespan_80pcnt_ints,
  "table-ascending-dose-rbc-lifespan.csv"
)
