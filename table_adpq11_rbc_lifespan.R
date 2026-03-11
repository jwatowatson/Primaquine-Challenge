#!/usr/bin/env -S Rscript --vanilla
#
# Print the estimated maximum drug effect for individual ADPQ 11 in the
# ascending dose study, and the impact this has on their RBC lifespan.
#

utils <- new.env()
sys.source("cmdstan_utils.R", envir = utils)
utils$load_packages(plot_libs = TRUE)
suppressPackageStartupMessages(library(tidyverse))


# Retrieve the dose_response() function from the Stan model.
model_file <- file.path(
  "Stan_models", "RBC_model_master_pop_free_weights_cmdstan.stan"
)
model <- utils$compile_model_with_exposed_functions(model_file)
model_fns <- model$functions


# Calculate the maximum reduction in RBC lifespan due to primaquine, for a
# specific individual.
#
# NOTE: RBC max lifespan is:
#       T_E_star + theta_rand[id[j], 8]
#
# NOTE: Max % decrease in RBC lifespan is:
#       inv_logit(logit_alpha + theta_rand[id[j], 6])
#
calculate_rbc_lifespan_reduction <- function(fit, patient_number) {
  T_E_effect <- paste0("theta_rand[", patient_number, ",8]")
  alpha_effect <- paste0("theta_rand[", patient_number, ",6]")
  variables <- c(
    "T_E_star", T_E_effect,
    "logit_alpha", alpha_effect
  )

  draws <- fit$draws(variables) |>
    as_draws_df() |>
    as.data.frame() |>
    pivot_longer(! c(.chain, .iteration, .draw))

  T_E_star <- filter(draws, name == "T_E_star")
  T_E_rand <- filter(draws, name == T_E_effect)
  T_E_values <- T_E_star$value + T_E_rand$value

  inv_logit <- function(x) {
    exp(x) / (1 + exp(x))
  }

  alpha_values <- inv_logit(
    (draws |>
       filter(name == "logit_alpha") |>
       pull(value))
    + (draws |>
         filter(name == alpha_effect) |>
         pull(value))
  )

  T_E_reduced <- T_E_values * (1 - alpha_values)

  tibble(
    T_E_initial = T_E_values,
    T_E_reduced = T_E_reduced,
    alpha = alpha_values
  )
}


subject_response <- function(fit, subject_id, dose_mgkg = 1) {
  alpha_effect_var <- paste0("theta_rand[", subject_id, ",6]")
  beta_effect_var <- paste0("theta_rand[", subject_id, ",7]")
  T_E_effect_var <- paste0("theta_rand[", subject_id, ",8]")

  subject_response <- cross_join(
    data.frame(effective_dose = dose_mgkg),
    utils$get_fit_draws_wide(
      fit, c("logit_alpha", "beta", "h", alpha_effect_var, beta_effect_var,
             "T_E_star", T_E_effect_var)
    ) |>
      rename(
        alpha_effect = !!alpha_effect_var,
        beta_effect = !!beta_effect_var,
        T_E_effect = !!T_E_effect_var
      ) |>
      mutate(
        logit_alpha = logit_alpha + alpha_effect,
        beta = beta * exp(beta_effect)
      )
  ) |>
    mutate(
      alpha = purrr:::pmap_dbl(
        list(effective_dose, logit_alpha, h, beta),
        model_fns$dose_response
      ),
      response = 100 * alpha,
      T_E_init = T_E_star + T_E_effect,
      # NOTE: the dose response `alpha` is the reduction in lifespan.
      T_E_reduced = T_E_init * (1 - alpha)
    )

  subject_response
}



# Load the model fit.
# NOTE: Job #3 is the fit against the ascending-dose study.
results_file <- file.path(
  "Rout", "pop_fit_free_weights_cmdstan_max_delay_9_job3.rds"
)
fit <- readRDS(results_file)

#
# NOTE: the `id` array comes from the job data, it is 1:23.
# We want individual 11, but their index is 3 because the IDs are sorted
# alphabetically rather than numerically (i.e., 1, 10, 11, 12, ...).
patient_number <- 3

tbl_redn <- calculate_rbc_lifespan_reduction(fit, patient_number)

tbl_redn_intervals <- tbl_redn |>
  pivot_longer(everything()) |>
  group_by(name) |>
  summarise(
    as_tibble_row(
      c(
        mean(value),
        quantile(
          value,
          c(0.025, 0.975, 0.050, 0.950, 0.100, 0.900, 0.500)
        )
      ),
      .name_repair = function(x) {
        names <- paste0("q", parse_number(x))
        names[is.na(parse_number(x))] <- "mean"
        names
      }
    )
  )

print(tbl_redn_intervals)


df_subj <- subject_response(fit, patient_number)

tbl_subj <- df_subj |>
  select(T_E_init, T_E_reduced, alpha) |>
  pivot_longer(everything()) |>
  group_by(name) |>
  summarise(
    as_tibble_row(
      c(
        mean(value),
        quantile(
          value,
          c(0.025, 0.975, 0.050, 0.950, 0.100, 0.900, 0.500)
        )
      ),
      .name_repair = function(x) {
        names <- paste0("q", parse_number(x))
        names[is.na(parse_number(x))] <- "mean"
        names
      }
    )
  )
print(tbl_subj)
