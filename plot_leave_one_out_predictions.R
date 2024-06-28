#!/usr/bin/env -S Rscript --vanilla
#
# Plot leave-one-out predictions for haemoglobin (g/dL) and reticulocytes (%),
# using each patient's initial haemoglobin measurement for the initial state,
# without adding random effects to this measurement.
#

main <- function() {
  library(tidyverse)
  library(cmdstanr)
  library(posterior)

  true_data <- collect_ground_truth()
  predictions <- collect_leave_one_out_predictions()

  p_hb <- plot_haemoglobin(true_data, predictions)
  p_retic <- plot_reticulocytes(true_data, predictions)

  # Save plots.
  retic_file <- "leave-one-out-retic-percent.png"
  cat("Writing", retic_file, "...")
  png(retic_file, width = 8, height = 12, units = "in", res = 150)
  print(p_retic)
  invisible(dev.off())
  cat("\n")

  hb_file <- "leave-one-out-haemoglobin.png"
  cat("Writing", hb_file, "...")
  png(hb_file, width = 8, height = 12, units = "in", res = 150)
  print(p_hb)
  invisible(dev.off())
  cat("\n")
}


plot_haemoglobin <- function(true_data, predictions) {
  blues <- scales::brewer_pal(palette = "Blues")(3)

  df_hb <- predictions |>
    filter(name == "hb")

  ggplot(df_hb) +
    geom_rect(
      aes(xmin = Start_Day, xmax = Final_Day, ymin = -Inf, ymax = Inf),
      true_data$regimen,
      fill = "#9f9f9f",
      alpha = 0.3
    ) +
    geom_ribbon(
      aes(day, ymin = lower, ymax = upper),
      fill = blues[2]
    ) +
    geom_vline(
      aes(xintercept = Start_Day),
      true_data$regimen,
      linetype = "dashed"
    ) +
    geom_vline(
      aes(xintercept = Final_Day),
      true_data$regimen,
      linetype = "dashed"
    ) +
    geom_line(
      aes(day, mean),
      colour = blues[3]
    ) +
    geom_point(
      aes(Study_Day, value, colour = name),
      true_data$hb
    ) +
    scale_x_continuous(
      "Days since start of primaquine",
      breaks = scales::breaks_width(7)
    ) +
    scale_y_continuous(
      "Haemoglobin (g/dL)",
      breaks = scales::breaks_width(3),
      minor_breaks = NULL
    ) +
    scale_colour_brewer(NULL, palette = "Dark2") +
    facet_wrap(
      ~ ID2,
      ncol = 4
    ) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "inside",
      legend.position.inside = c(1, 0),
      legend.justification = c(1, 0)
    )
}


plot_reticulocytes <- function(true_data, predictions) {
  blues <- scales::brewer_pal(palette = "Blues")(3)

  df_retic <- predictions |>
    filter(name == "retic")

  ggplot(df_retic) +
    geom_rect(
      aes(xmin = Start_Day, xmax = Final_Day, ymin = -Inf, ymax = Inf),
      true_data$regimen,
      fill = "#9f9f9f",
      alpha = 0.3
    ) +
    geom_ribbon(
      aes(day, ymin = lower, ymax = upper),
      fill = blues[2]
    ) +
    geom_vline(
      aes(xintercept = Start_Day),
      true_data$regimen,
      linetype = "dashed"
    ) +
    geom_vline(
      aes(xintercept = Final_Day),
      true_data$regimen,
      linetype = "dashed"
    ) +
    geom_line(
      aes(day, mean),
      colour = blues[3]
    ) +
    geom_point(
      aes(Study_Day, value, colour = name),
      true_data$retic
    ) +
    scale_x_continuous(
      "Days since start of primaquine",
      breaks = scales::breaks_width(7)
    ) +
    scale_y_continuous(
      "Reticulocytes (%)",
      limits = c(0, NA),
      breaks = scales::breaks_width(5),
      minor_breaks = NULL
    ) +
    scale_colour_brewer(NULL, palette = "Dark2") +
    facet_wrap(
      ~ ID2,
      ncol = 4
    ) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "inside",
      legend.position.inside = c(1, 0),
      legend.justification = c(1, 0)
    )
}


get_dose_weights <- function(fit) {
  dose_weights <- fit$draws("dose_weights", format = "df") |>
    as_tibble() |>
    select(starts_with("dose_weights["))

  dose_weights_mat <- dose_weights |>
    as.matrix() |>
    unname()

  dose_weights_vec <- purrr::map(
    seq_len(nrow(dose_weights_mat)),
    function(i) dose_weights_mat[i, ]
  )

  tibble(dose_weights = dose_weights_vec)
}


sample_individual_effects <- function(fit, seed) {
  omega <- fit$draws("L_Omega", format = "df") |>
    as_tibble() |>
    select(starts_with("L_Omega[")) |>
    as.matrix() |>
    unname()

  omega_list <- purrr::map(
    seq_len(nrow(omega)), function(i) matrix(omega[i, ], ncol = 8)
  )

  sigmasq <- fit$draws("sigmasq_u", format = "df") |>
    as_tibble() |>
    select(starts_with("sigmasq_u[")) |>
    as.matrix() |>
    unname()

  sigmasq_list <- purrr::map(
    seq_len(nrow(sigmasq)), function(i) sigmasq[i, ]
  )

  effect_names <- c(
    "Hb_star", "alpha_diff1", "alpha_delta1", "alpha_diff2", "alpha_delta2",
    "logit_alpha", "beta", "T_E_star"
  )

  prev_random_state <- .Random.seed
  set.seed(seed)

  effects_list <- purrr::map(
    seq_len(length(omega_list)),
    function(i) {
      sig <- sigmasq_list[[i]]
      omg <- omega_list[[i]]

      # NOTE: the following was provided by James Watson.
      L <- diag(sig) %*% omg
      Epsilon <- L %*% t(L)
      random_effects <- mvtnorm::rmvnorm(n = 1, sigma = Epsilon) |>
        as.numeric()

      names(random_effects) <- effect_names
      as_tibble_row(random_effects)
    }
  )

  .Random.seed <- prev_random_state

  do.call(bind_rows, effects_list)
}


get_draws <- function(fit, individual_effects = TRUE) {
  dose_weights_tbl <- get_dose_weights(fit)

  param_names <- c(
    "Hb_star", "alpha_diff1", "alpha_delta1", "alpha_diff2", "alpha_delta2",
    "logit_alpha", "h", "beta", "T_E_star", "log_k"
  )

  draws <- fit$draws(variables = param_names, format = "df") |>
    as_tibble() |>
    mutate(dose_weights = dose_weights_tbl$dose_weights)

  # Add random effects to each draw.
  if (individual_effects) {
    # Sample the random effects for each draw.
    effects <- sample_individual_effects(fit, seed = 12345)

    # Apply the random effects to each parameter.
    draws$Hb_star <- draws$Hb_star + effects$Hb_star
    draws$alpha_diff1 <- draws$alpha_diff1 * exp(effects$alpha_diff1)
    draws$alpha_delta1 <- draws$alpha_delta1 * exp(effects$alpha_delta1)
    draws$alpha_diff2 <- draws$alpha_diff2 * exp(effects$alpha_diff2)
    draws$alpha_delta2 <- draws$alpha_delta2 * exp(effects$alpha_delta2)
    draws$logit_alpha <- draws$logit_alpha + effects$logit_alpha
    draws$beta <- draws$beta * exp(effects$beta)
    draws$T_E_star <- draws$T_E_star + effects$T_E_star
  }

  draws
}


build_forward_sim_args <- function(draws, job_data) {
  data_args <- c(
    "T_nmblast", "T_retic", "T_RBC_max", "T_transit_steady_state", "K_weights",
    "sigma"
  )
  args <- job_data[data_args]
  args$nComp_sim <- 29

  for (name in names(args)) {
    draws[[name]] <- args[[name]]
  }

  arg_names <- c(
    "Hb_star", "alpha_diff1", "alpha_delta1", "alpha_diff2", "alpha_delta2",
    "logit_alpha", "h", "beta", "T_E_star", "log_k", "nComp_sim", "T_nmblast",
    "T_retic", "T_RBC_max", "T_transit_steady_state", "dose_weights",
    "K_weights", "sigma"
  )

  draws |>
    select(all_of(arg_names))
}


get_data_inputs <- function(job_data) {
  # NOTE: need to adjust the working directory for `create_job_data()`.
  curr_wd <- getwd()
  setwd("..")

  utils <- new.env()
  sys.source(file.path("cmdstan_utils.R"), envir = utils)
  utils$load_packages(plot_libs = TRUE)

  # NOTE: the following arguments are data:
  data_args <- c(
    "T_nmblast", "T_retic", "T_RBC_max", "T_transit_steady_state", "K_weights",
    "sigma"
  )
  job_data <- utils$create_job_data(
    job_number = 3,
    max_dose_delay = 9,
    quiet = TRUE
  )
  args <- job_data[data_args]
  args$nComp_sim <- 29

  # NOTE: restore the original working directory.
  setwd(curr_wd)

  args
}


simulate_model <- function(forward_sim, job_data, draws) {
  sim_args <- build_forward_sim_args(draws, job_data)

  result_tibbles <- lapply(
    seq_len(nrow(sim_args)),
    function(args_ix) {
      args <- sim_args[args_ix, ] |> as.list()
      args$dose_weights <- args$dose_weights[[1]]
      args$drug_regimen <- job_data$drug_regimen_pred
      # NOTE: this is where we override the random effects and keep the
      # starting haemoglobin fixed at the initial measurement.
      args$Hb_star <- job_data$Y_true_haemocue
      out <- do.call(forward_sim, args)
      tibble(
        sample = args_ix,
        day = seq_len(dim(out)[2]) - 1,
        hb = out[1, ],
        retic = out[2, ],
        eff_dose = out[3, ]
      )
    }
  )

  results <- bind_rows(result_tibbles)

  # Return mean and 90% credible intervals.
  results |>
    pivot_longer(c(hb, retic, eff_dose)) |>
    group_by(day, name) |>
    summarise(
      mean = mean(value),
      lower = quantile(value, 0.05),
      upper = quantile(value, 0.95),
      .groups = "drop"
    )

}


factor_patients_by_number <- function(patient_ids) {
  unique_ids <- unique(patient_ids)
  # Strip the "ADPQ " prefix and convert to integers.
  patient_numbers <- as.integer(substring(unique_ids, 5))
  # Sort the patients by number, rather than alphabetically.
  patient_order <- unique_ids[order(patient_numbers)]
  # Return an ordered factor.
  factor(patient_ids, levels = patient_order, ordered = TRUE)
}


collect_ground_truth <- function() {
  data_env <- new.env()
  load("Data/RBC_model_data.RData", envir = data_env)
  study_data <- data_env$PQdat |> filter(study == "Part1")

  # Extract the data collected from each individual.
  df_true_retic <- study_data |>
    select(ID2, Study_Day, CBC_retic, Manual_retic) |>
    rename_with(~ gsub("_retic$", "", .x)) |>
    pivot_longer(! c(ID2, Study_Day)) |>
    mutate(ID2 = factor_patients_by_number(ID2)) |>
    filter(! is.na(value))

  df_true_hb <- study_data |>
    select(ID2, Study_Day, Haemocue_hb, CBC_hb) |>
    rename_with(~ gsub("_hb$", "", .x)) |>
    pivot_longer(! c(ID2, Study_Day)) |>
    mutate(ID2 = factor_patients_by_number(ID2)) |>
    filter(! is.na(value))

  # Extract the first and last day of treatment.
  df_regimen <- study_data |>
    filter(dosemgkg > 0) |>
    mutate(ID2 = factor_patients_by_number(ID2)) |>
    group_by(ID2) |>
    summarise(
      Start_Day = min(Study_Day),
      Final_Day = max(Study_Day),
    )

  list(retic = df_true_retic, hb = df_true_hb, regimen = df_regimen)
}


collect_leave_one_out_predictions <- function(forward_sim = NULL) {
  utils <- new.env()
  sys.source("cmdstan_utils.R", envir = utils)

  if (is.null(forward_sim)) {
    model_file <- file.path(
      "Stan_models", "RBC_model_master_pop_free_weights_cmdstan.stan"
    )
    model <- utils$compile_model_with_exposed_functions(model_file)
    forward_sim <- model$functions$forwardsim
  }

  max_dose_delay <- 9
  prefix <- paste0(
    "pop_fit_mechanistic_G6PD_cmdstan_max_delay_",
    max_dose_delay
  )

  jobs <- utils$create_job_data_ascending_dose_leave_one_out(
    max_dose_delay, quiet = TRUE
  )

  data_env <- new.env()
  load("Data/RBC_model_data.RData", envir = data_env)
  PQdat <- data_env$PQdat
  patient_ids <- PQdat |> filter(study == "Part1") |> pull(ID2) |> unique()

  lapply(
    seq_along(jobs),
    function(job_ix) {
      results_file <- file.path(
        "Rout", paste0(prefix, "_leave_one_out_", job_ix, ".rds")
      )
      fit <- readRDS(results_file)

      job_data <- jobs[[job_ix]]
      draws <- get_draws(fit)

      simulate_model(forward_sim, job_data, draws) |>
        mutate(ID2 = patient_ids[job_ix])
    }
  ) |>
    bind_rows() |>
    # Order patients numerically, noting that ADPQ 18 was excluded.
    mutate(ID2 = factor_patients_by_number(ID2))
}


call_main <- function(script_name, main) {
  if (interactive()) {
    return()
  }

  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=.+", args, value = TRUE)
  if (length(file_arg) != 1) {
    warning("Found ", length(file_arg), " '--file=' arguments")
    return()
  }

  file_path <- substring(file_arg, 8)
  file_name <- basename(file_path)
  if (file_name == script_name) {
    main_args <- commandArgs(trailingOnly = TRUE)
    status <- main(main_args)
    if (! is.numeric(status)) {
      status <- 0
    }
    quit(status = status)
  }
}

call_main("plot_leave_one_out_predictions.R", main)
