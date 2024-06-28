#!/usr/bin/env -S Rscript --vanilla
#
# Fit the RBC model with independent dose-weighting parameters.
#

main <- function(args) {
  utils <- new.env()
  sys.source("cmdstan_utils.R", envir = utils)
  utils$load_packages(plot_libs = TRUE)

  output_dir <- "Rout"
  results_regex <- "^.*_leave_one_out_(.+)\\.rds$"
  results_files <- list.files(
    path = output_dir,
    pattern = results_regex,
    full.names = TRUE
  )
  left_out_ixs <- as.integer(sub(results_regex, "\\1",
                                 basename(results_files)))

  truth_dfs <- collect_ground_truth()
  ids <- unique(truth_dfs$hb$ID2)[left_out_ixs]

  # Retrieve the dose_response() function from the Stan model.
  model_file <- file.path(
    "Stan_models", "RBC_model_master_pop_free_weights_cmdstan.stan"
  )
  model <- utils$compile_model_with_exposed_functions(model_file)
  dose_response_fn <- model$functions$dose_response

  # Calculate the dose responses for each fit.
  df_loo_responses <- dose_responses(results_files, ids, dose_response_fn)

  # Calculate the dose responses when fitting to the entire study.
  whole_study_file <- file.path(
    "Rout",
    "pop_fit_free_weights_cmdstan_max_delay_9_job3.rds"
  )
  df_net_response <- dose_responses(whole_study_file, "", dose_response_fn) |>
    mutate(ID2 = NULL)

  # Plot the dose responses.
  p_resp <- plot_dose_responses(df_loo_responses, df_net_response)

  resp_file <- "leave-one-out-dose-responses.png"
  cat("Writing", resp_file, "...")
  png(resp_file, width = 8, height = 12, units = "in", res = 150)
  print(p_resp)
  invisible(dev.off())
  cat("\n")

  invisible(0)
}


plot_dose_responses <- function(df_loo_responses, df_net_response) {
  ggplot() +
    geom_ribbon(
      aes(effective_dose, ymin = Lower, ymax = Upper),
      data = df_loo_responses,
      fill = "#9f9f9f"
    ) +
    geom_line(
      aes(effective_dose, Median),
      data = df_loo_responses
    ) +
    geom_line(
      aes(effective_dose, Median),
      data = df_net_response,
      colour = "red",
      linetype = "dashed"
    ) +
    geom_line(
      aes(effective_dose, Lower),
      data = df_net_response,
      colour = "red",
      linetype = "dashed"
    ) +
    geom_line(
      aes(effective_dose, Upper),
      data = df_net_response,
      colour = "red",
      linetype = "dashed"
    ) +
    xlab("Dose (mg/kg)") +
    ylab("Reduction in RBC lifespan (%)") +
    facet_wrap(~ ID2, ncol = 4) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      panel.grid.minor = element_blank()
    )
}


collect_ground_truth <- function() {
  data_env <- new.env()
  load("Data/RBC_model_data.RData", envir = data_env)
  study_data <- data_env$PQdat |> filter(study == "Part1")

  # Extract the data collected from each individual.
  # TODO: do we need to apply a correction to the reticulocyte data?
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


factor_patients_by_number <- function(patient_ids) {
  unique_ids <- unique(patient_ids)
  # Strip the "ADPQ " prefix and convert to integers.
  patient_numbers <- as.integer(substring(unique_ids, 5))
  # Sort the patients by number, rather than alphabetically.
  patient_order <- unique_ids[order(patient_numbers)]
  # Return an ordered factor.
  factor(patient_ids, levels = patient_order, ordered = TRUE)
}


dose_responses <- function(files, ids, dose_response_fn) {
  df_draws <- dose_response_draws(files, ids)
  dose_response_intervals(df_draws, dose_response_fn)
}


dose_response_draws <- function(results_files, ids) {
  all_draws <- list()

  for (ix in seq_along(results_files)) {
    fit <- readRDS(results_files[ix])

    # Extract the dose-response parameter draws.
    draws <- as_draws_df(fit$draws(c("logit_alpha", "h", "beta"))) |>
      as.data.frame() |>
      mutate(ID2 = ids[ix]) |>
      select(! c(.chain, .iteration))

    all_draws[[length(all_draws) + 1]] <- draws
  }

  bind_rows(all_draws)
}


dose_response_intervals <- function(df_draws, dose_response_fn,
                                    low = 0.05, high = 0.95) {
  dose_response_df <- cross_join(
    # Effective Primaquine dose (mg/kg).
    data.frame(effective_dose = seq(0, 1, length.out = 200)),
    df_draws
  ) |>
    mutate(
      response = 100 * purrr:::pmap_dbl(
        list(effective_dose, logit_alpha, h, beta),
        dose_response_fn
      )
    )

  # Calculate response intervals for each leave-one-out fit.
  dose_response_df |>
    group_by(effective_dose, ID2) |>
    summarise(
      Mean = mean(response),
      Median = median(response),
      Lower = quantile(response, probs = low),
      Upper = quantile(response, probs = high),
      .groups = "drop"
    )
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

call_main("plot_free_weights_leave_one_out.R", main)
