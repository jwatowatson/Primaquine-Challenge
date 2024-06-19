#!/usr/bin/env -S Rscript --vanilla
#
# Plot the estimated relationships between
#
# (a) Administered dose and effective dose; and
# (b) Effective dose and reduction in RBC lifespan
#
# in a single figure.
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
  # NOTE: Job #3 is the fit against the ascending-dose study.
  results_file <- file.path(
    "Rout", "pop_fit_free_weights_cmdstan_max_delay_9_job3.rds"
  )
  fit <- readRDS(results_file)

  resp <- new.env()
  sys.source("plot_dose_response.R", envir = resp)

  df_mean_resp <- resp$mean_dose_response(fit)
  df_subj_resp <- resp$subject_responses(fit)
  p_response <- utils$plot_subject_dose_responses(
    df_mean_resp, df_subj_resp
  ) +
    theme_bw()

  delay <- new.env()
  sys.source("plot_effective_dose_delay.R", envir = delay)

  dose_delay <- delay$get_effective_dose_delay(utils, model_fns, fit)
  p_effective_dose <- delay$plot_effective_dose_delay(dose_delay)

  ggsave(
    file = "dose-response-and-delay.png",
    gridExtra::arrangeGrob(
      p_effective_dose +
        ggtitle("A: Effective dose") +
        # NOTE: units specify top, right, bottom, left
        theme(plot.margin = unit(c(5.5, 22, 5.5, 5.5), "points")),
      p_response +
        ggtitle("B: Dose response") +
        # NOTE: units specify top, right, bottom, left
        theme(plot.margin = unit(c(5.5, 5.5, 5.5, 22), "points")),
      nrow = 1,
      ncol = 2
    ),
    width = 8,
    height = 4,
    units = "in"
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

call_main("plot_dose_response_and_delay.R", main)
