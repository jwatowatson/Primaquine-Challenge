#!/usr/bin/env -S Rscript --vanilla
#
# Plot leave-one-out fits for the RBC model with independent dose-weighting
# parameters.
#
# USAGE:
#
#     ./plot_comparison_of_leave_one_out_fits.R
#

main <- function(args) {
  utils <- new.env()
  sys.source("cmdstan_utils.R", envir = utils)
  utils$load_packages(plot_libs = TRUE)

  max_delay <- 9

  plot_comparison_of_leave_one_out_fits(utils, max_delay)
}

plot_comparison_of_leave_one_out_fits <- function(utils, max_delay) {
  re_match <- paste0(
    "pop_fit_mechanistic_G6PD_cmdstan_max_delay_",
    max_delay,
    "_leave_one_out_(.*)\\.rds"
  )

  # Identify the results file for each completed fit.
  results_files <- list.files(
    path = "Rout",
    pattern = re_match,
    full.names = TRUE
  )
  job_numbers <- as.numeric(sub(re_match, "\\1", basename(results_files)))

  # Identify the left-out individual for each completed fit.
  data_env <- new.env()
  load("Data/RBC_model_data.RData", envir = data_env)
  left_out_ids <- data_env$PQdat |>
    filter(study == "Part1") |>
    pull(ID2) |>
    unique()
  job_names <- left_out_ids

  # Highlight the results when specific individuals were left out.
  highlight_ids <- c("ADPQ 11", "ADPQ 14", "ADPQ 17", "ADPQ 23", "ADPQ 24")

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
  ) |>
    mutate(
      name = factor(name, levels = params, ordered = TRUE),
      highlight = case_when(
        job %in% highlight_ids ~ job,
        TRUE ~ ""
      )
    )

  p <- ggplot() +
    geom_density(
      aes(value, group = job),
      colour = "#7f7f7f",
      linewidth = 0.5,
      filter(draws, highlight == "")
    ) +
    geom_density(
      aes(value, group = job, colour = highlight),
      linewidth = 1,
      filter(draws, highlight != "")
    ) +
    scale_colour_brewer(
      name = "Left out",
      palette = "Set1"
    ) +
    facet_wrap(~ name, scales = "free", ncol = 1) +
    xlab(NULL) +
    ylab(NULL)

  png(
    "comparison_of_leave_one_out_fits.png",
    width = 12, height = 6, unit = "in", res = 150
  )
  print(p)
  invisible(dev.off())
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

call_main("plot_comparison_of_leave_one_out_fits.R", main)
