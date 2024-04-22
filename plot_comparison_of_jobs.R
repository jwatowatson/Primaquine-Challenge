#!/usr/bin/env -S Rscript --vanilla
#
# Plot fits for the RBC model with independent dose-weighting parameters.
#
# USAGE:
#
#     ./plot_comparison_of_jobs.R <MAX_DELAY>
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
  ) |>
    mutate(
      name = factor(name, levels = params, ordered = TRUE)
    )

  p <- ggplot() +
    geom_violin(aes(job, value), draws) +
    facet_wrap(~ name, scales = "free_y") +
    xlab(NULL) +
    ylab(NULL) +
    theme(axis.text.x = element_text(size = rel(0.75)))

  png(
    "comparison_of_job_fits.png",
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

call_main("plot_comparison_of_jobs.R", main)
