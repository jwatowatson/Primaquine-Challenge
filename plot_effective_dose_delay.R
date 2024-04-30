#!/usr/bin/env -S Rscript --vanilla
#
# Plot the delay in reaching the steady-state effective dose, under a constant
# drug regimen (i.e., the same dose is given each day).
#

main <- function(args) {
  utils <- new.env()
  sys.source("cmdstan_utils.R", envir = utils)
  utils$load_packages(plot_libs = TRUE)

  # Retrieve the compute_effective_dose() function from the Stan model.
  model_file <- file.path(
    "Stan_models", "RBC_model_master_pop_free_weights_cmdstan.stan"
  )
  model <- utils$compile_model_with_exposed_functions(model_file)
  model_fns <- model$functions

  # Load the model fit.
  results_file <- file.path(
    "Rout", "pop_fit_free_weights_cmdstan_max_delay_9_job4.rds"
  )
  fit <- readRDS(results_file)

  # Calculate the effective dose at each day, for a constant drug regimen.
  dose_delay <- get_effective_dose_delay(utils, model_fns, fit)
  p_effective_dose <- plot_effective_dose_delay(dose_delay)

  # Save the plot.
  plot_file <- "ascending-dose-fit-effective-dose-delay.png"
  cat("Writing", plot_file, "...")
  png(plot_file, width = 4.5, height = 4.5, units = "in", res = 150)
  print(p_effective_dose)
  invisible(dev.off())
  cat("\n")

  invisible(0)
}

get_effective_dose_delay <- function(
  utils, model_fns, fit, interval = 0.95, dose_mgkg = 0.5
) {
  # Retrieve the dose-delay weights for each sample.
  weights <- utils$get_fit_dose_delay_weights(fit) |>
    select(delay, value, .draw) |>
    rename(draw = .draw)

  # Define a constant drug regimen.
  n_days <- weights |> pull(delay) |> unique() |> length()
  K_weights <- n_days
  regimen <- rep(dose_mgkg, n_days)

  # Calculate the effective dose at each day, for each sample.
  effective_doses <- list()
  for (draw in unique(weights$draw)) {
    sample_weights <- weights |> filter(draw == !!draw) |> pull(value)
    df_eff <- data.frame(
      draw = rep(draw, n_days),
      day = seq(n_days) - 1,
      effective_dose = c(model_fns$compute_effective_dose(
        regimen, n_days, sample_weights, K_weights
      ))
    )
    effective_doses[[length(effective_doses) + 1]] <- df_eff
  }
  df_effective_doses <- bind_rows(effective_doses)

  # Calculate median and intervals at each day.
  pr_lower <- 0.5 * (1 - interval)
  pr_upper <- 1 - pr_lower

  df_intervals <- df_effective_doses |>
    group_by(day) |>
    summarise(
      Mean = mean(effective_dose),
      Median = median(effective_dose),
      Lower = quantile(effective_dose, probs = pr_lower),
      Upper = quantile(effective_dose, probs = pr_upper),
      .groups = "drop"
    )

  list(
    dose_mgkg = dose_mgkg,
    n_days = n_days,
    samples = df_effective_doses,
    intervals = df_intervals
  )
}

plot_effective_dose_delay <- function(dose_delay) {
  ggplot() +
    geom_line(
      aes(day, effective_dose, group = draw),
      dose_delay$samples,
      colour = "#afafaf",
      alpha = 0.1
    ) +
    geom_errorbar(
      aes(day, ymin = Lower, ymax = Upper),
      dose_delay$intervals,
      linewidth = 0.5,
    ) +
    geom_point(
      aes(day, Median),
      dose_delay$intervals,
      size = 2,
    ) +
    geom_hline(
      yintercept = dose_delay$dose_mgkg,
      linetype = "dashed"
    ) +
    scale_x_continuous(
      "Time (days)",
      breaks = seq(to = dose_delay$n_days, by = 3)
    ) +
    scale_y_continuous(
      "Effective Dose (% of administered dose)",
      breaks = (0:4) * dose_delay$dose_mgkg / 4,
      labels = paste0(100 * (0:4) / 4, "%")
    ) +
    theme_bw()
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

call_main("plot_effective_dose_delay.R", main)
