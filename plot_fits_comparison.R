#!/usr/bin/env -S Rscript --vanilla
#
# Plot fits for the RBC model with independent dose-weighting parameters.
#
# USAGE:
#
#     ./plot_fits_comparison.R <MAX_DELAY>
#
main <- function(args) {
  utils <- new.env()
  sys.source("cmdstan_utils.R", envir = utils)
  utils$load_packages(plot_libs = TRUE)

  if (length(args) == 0) {
    # Default value
    max_delay <- 9
  } else if (length(args) == 1) {
    max_delay <- as.integer(args[1])
  } else {
    stop("Invalid arguments: ", args)
  }

  plot_comparison_of_jobs(utils, max_delay)
}

plot_comparison_of_jobs <- function(utils, max_delay) {
  re_match <- paste0(
    "pop_fit_free_weights_cmdstan_max_delay_",
    max_delay,
    "_job(.*)\\.rds"
  )

  # Provide a description for each job.
  job_names <- c(
    "Both studies,\ncommon individuals",
    "Both studies,\ndifferent individuals",
    "Ascending-dose",
    "Single-dose"
  )

  # Identify the results file for each completed fit.
  results_files <- list.files(
    path = "Rout",
    pattern = re_match,
    full.names = TRUE
  )
  job_numbers <- as.numeric(sub(re_match, "\\1", basename(results_files)))

  params <- c("logit_alpha", "beta", "h", "Hb_star", "T_E_star")
  draws <- bind_rows(
    lapply(
      seq_along(results_files),
      function(i) {
        fit <- readRDS(results_files[i])
        utils$get_fit_draws_long(fit, params) |>
          mutate(
            job = factor(
              job_names[job_numbers[i]],
              levels = job_names,
              ordered = TRUE
            )
          )
      }
    )
  )

  # Replace the "logit_alpha" values with "alpha" values (%).
  all_draws <- bind_rows(
    draws |>
      filter(name != "logit_alpha"),
    draws |>
      filter(name == "logit_alpha") |>
      mutate(
        value = 100 * exp(value) / (1 + exp(value)),
        name = "alpha"
      )
  )

  # Create descriptive labels for each facet.
  #
  # NOTE: use plotmath expressions and `label_parsed` to display the
  # descriptive label and the notation used in the manuscript.
  labels <- c(
    "alpha" = "paste('Max RBC lifespan reduction (', alpha, ', %)')",
    "beta" = "paste('Half-maximal effect dose (', beta, ', mg/kg)')",
    "h" = "paste('Dose-response slope (', h, ')')",
    "Hb_star" = "paste('Steady-state haemoglobin (', 'Hb'*'*', ', g/dL)')",
    "T_E_star" = "paste('Steady-state RBC lifespan (', T[E]*'*', ', days)')"
  )

  labelled_draws <- all_draws |>
    mutate(
      label = factor(
        name,
        levels = names(labels),
        labels = unname(labels),
        ordered = TRUE
      )
    )

  # Adjust the y-axis limits for individual facets.
  expand_limits <- data.frame(
    label = unname(labels[c("alpha", "h")]),
    value = 0,
    job = "Ascending-dose"
  )

  p <- ggplot() +
    geom_violin(aes(job, value), labelled_draws) +
    geom_blank(aes(job, value), expand_limits) +
    facet_wrap(~ label, scales = "free_y", labeller = label_parsed) +
    xlab(NULL) +
    ylab(NULL) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = rel(0.75)),
      strip.background = element_blank()
    )

  plot_file <- "fits_comparison.png"
  cat("Writing", plot_file, "...")
  png(plot_file, width = 12, height = 6, unit = "in", res = 150)
  print(p)
  invisible(dev.off())
  cat(" done\n")

  # Print mean values and 95% credible intervals for the two studies.
  labelled_draws |>
    filter(job %in% c("Ascending-dose", "Single-dose")) |>
    group_by(job, name) |>
    summarise(
      mean = mean(value),
      lower = quantile(value, 0.025),
      upper = quantile(value, 0.975),
      .groups = "drop"
    ) |>
    print()
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

call_main("plot_fits_comparison.R", main)
