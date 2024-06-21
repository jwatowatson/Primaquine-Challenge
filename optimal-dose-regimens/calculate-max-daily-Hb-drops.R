#!/usr/bin/env -S Rscript --vanilla

main <- function(args = NULL) {
  library(dplyr, warn.conflicts = FALSE)
  library(future)
  library(furrr)
  library(readr)
  library(tibble)
  library(tidyr)

  # Retrieve dose regimen details from `settings.R`.
  settings <- new.env()
  sys.source("settings.R", envir = settings)

  # Load draws from the posterior and sample individual random effects.
  fit <- settings$load_fit()
  n_draws <- 1000
  draws <- get_draws(fit)
  data_inputs <- settings$get_data_inputs()
  forward_sim_params <- build_forward_sim_args(draws, data_inputs)

  # Load the Stan model simulation function.
  forward_sim <- settings$get_forward_sim_function()
  cat("done\n")

  # NOTE: We could use the `callr` strategy provided by the `future.callr`
  # library (https://future.callr.futureverse.org/). But it seems that this
  # might require compiling the model in each worker process --- otherwise we
  # see the error "NULL value passed as symbol address" --- and it isn't clear
  # how to allow each worker process to compile the model separately and in
  # isolation.
  plan(strategy = "multicore", gc = TRUE)

  for (scenario in settings$scenarios) {
    num_days <- scenario$duration
    num_chunks <- scenario$num_chunks
    evaluation_file_prefix <- scenario$evaluation_file_prefix

    cat(
      "Calculating max daily Hb drops for ", num_days, "-day regimens\n",
      sep = ""
    )

    df_regimens <- readr::read_csv(
      scenario$regimens_file, show_col_types = FALSE
    )

    # Divide the regimens into chunks.
    chunks <- seq(
      from = 1,
      to = nrow(df_regimens) + 1,
      length.out = num_chunks + 1
    ) |>
      round() |>
      as.integer()

    df_regimen_list <- lapply(
      seq_len(num_chunks),
      function(ix) {
        start <- chunks[ix]
        until <- chunks[ix + 1] - 1
        df_regimens[start:until, ]
      }
    )

    for (ix in seq_along(df_regimen_list)) {
      evaluation_file <- paste0(evaluation_file_prefix, ix, ".rds")

      if (file.exists(evaluation_file)) {
        cat("Skipping part", ix, "...", evaluation_file, "exists\n")
        next
      }

      df_reg <- df_regimen_list[[ix]]
      n_regs <- nrow(df_reg)

      cat("Evaluating ", n_regs, " regimens (part ", ix, ")\n", sep = "")

      start <- Sys.time()
      hb_drops <- calculate_hb_drops_parallel(
        forward_sim, df_reg, n_regs, forward_sim_params, n_draws
      )
      finish <- Sys.time()
      print(finish - start)

      cat("Saving", evaluation_file, "...")
      saveRDS(hb_drops, evaluation_file, compress = "xz")
      cat(" done\n")

      # NOTE: running multiple passes of this loop body results in the main R
      # process being killed, so we instead terminate after completing a
      # single chunk, and rely on running this script once for each chunk.
      return()
    }
  }
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


build_forward_sim_args <- function(draws, data_inputs) {
  for (name in names(data_inputs)) {
    draws[[name]] <- data_inputs[[name]]
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


calculate_hb_drops_parallel <- function(
  forward_sim, df_regimens, num_regimens, df_draws, num_draws, nComp_sim = 29
) {
  df_regs <- df_regimens[1:num_regimens, ]

  draw_ixs <- seq(
    from = 1,
    to = nrow(df_draws),
    length.out = num_draws
  ) |>
    round()
  df_in <- df_draws[draw_ixs, ]

  # Retrieve dose regimen details from `settings.R`.
  settings <- new.env()
  sys.source("settings.R", envir = settings)

  # Sample body weights from a normal distribution around a reference weight.
  prev_random_state <- .Random.seed
  set.seed(12345)
  body_weights <- rnorm(
    n = num_draws,
    mean = settings$weight_kg,
    sd = settings$weight_sd
  )
  .Random.seed <- prev_random_state

  cat(
    "Body weight ranges from",
    min(body_weights), "to", max(body_weights),
    "kg\n"
  )

  column_names <- paste0("draw_", seq_len(num_draws))

  furrr::future_map_dfr(
    seq_len(nrow(df_regs)),
    function(regimen_ix) {
      hb_drops <- 0 * seq_len(num_draws)
      names(hb_drops) <- column_names

      # NOTE: pad the drug regimen to cover the entire simulation (29 days).
      regimen_units <- as.numeric(df_regs[regimen_ix, ])
      regimen_zeros <- 0 * seq_len(nComp_sim - length(regimen_units))
      regimen_units <- c(regimen_units, regimen_zeros)

      for (draw_ix in seq_len(nrow(df_in))) {
        # Adjust the dose to account for variations in body weight.
        weight_kg <- body_weights[draw_ix]
        dose_scale <- settings$dose_unit_mg / weight_kg
        regimen_mg_kg <- dose_scale * regimen_units

        draw <- df_in[draw_ix, ] |> as.list()
        draw$dose_weights <- draw$dose_weights[[1]]
        draw$drug_regimen <- regimen_mg_kg

        out <- do.call(forward_sim, draw)

        hb_daily_change <- diff(out[1, ])
        hb_drops[draw_ix] <- min(hb_daily_change)
      }

      as_tibble_row(hb_drops)
    },
    .options = furrr_options(seed = TRUE)
  )
}


call_main <- function(script_name) {
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

call_main("calculate-max-daily-Hb-drops.R")
