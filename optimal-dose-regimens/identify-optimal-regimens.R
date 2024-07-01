#!/usr/bin/env -S Rscript --vanilla

main <- function(args = NULL) {
  library(dplyr, warn.conflicts = FALSE)
  library(ggplot2)
  library(readr)
  library(tibble)
  library(tidyr)

  # Retrieve dose regimen details from `settings.R`.
  settings <- new.env()
  sys.source("settings.R", envir = settings)

  # Load the regimens that were given in the ascending dose study.
  ascending_study <- load_ascending_study_regimens()

  # Collect the results for each regimen.
  results_file <- "evaluate-thresholds-results.rds"
  plot_tbl <- load_or_collect_results(results_file)

  # Identify optimal regimens for a subset of thresholds.
  thresholds <- seq(from = 0.5, to = 1.0, by = 0.1)
  optimal_dfs <- load_optimal_regimens(settings, thresholds)

  # Predict daily haemoglobin and reticulocyte % for optimal regimens.
  forward_sim <- settings$get_forward_sim_function()
  tbl_ints <- predict_outcomes(
    settings,
    optimal_dfs$long |> filter(threshold == 1),
    forward_sim = forward_sim
  )

  # For the optimal dose regimens for each threshold, plot the fraction of
  # individuals whose largest haemoglobin drop exceeded the threshold.
  p_cmp <- plot_frac_exceeding_vs_threshold(plot_tbl)

  # Plot the number of optimal regimens for each threshold.
  p_num_optimal <- plot_number_of_optimal_regimens(plot_tbl)

  # Plot optimal dose regimens for a small number of thresholds.
  p_optimal <- plot_optimal_regimen_daily_doses(settings, optimal_dfs$long)

  # Plot the optimal regimens against all of the dose regimens from the
  # ascending-dose study in a separate facet for each individual.
  p_vs_study <- plot_optimal_daily_doses_vs_ascending(
    settings, optimal_dfs$long, ascending_study
  )

  # Plot the optimal regimens against all of the dose regimens from the
  # ascending-dose study in a separate facet for each individual, in a wider
  # format for, e.g., presentations.
  p_vs_study_wide <- p_vs_study +
    facet_wrap(
      ~ ID_num,
      ncol = 6,
      labeller = as_labeller(function(x) paste("ADPQ", x))
    )

  # Plot cumulative doses for the optimal regimens against all those of the
  # dose regimens from the ascending-dose study in a single pane.
  p_cum <- plot_optimal_regimens_vs_study(optimal_dfs$long, ascending_study)
  p_hb <- plot_optimal_output_intervals(
    tbl_ints |> filter(measure == "Haemoglobin (g/dL)")
  ) +
    scale_y_continuous(
      "Haemoglobin (g/dL)",
      breaks = c(8, 10, 12, 14, 16)
    )
  p_retic <- plot_optimal_output_intervals(
    tbl_ints |> filter(measure == "Reticulocyte (%)")
  ) +
    scale_y_continuous(
      "Reticulocyte (%)",
      breaks = c(0, 5, 10, 15, 20),
      limits = c(0, NA)
    )

  grob_outputs <- gridExtra::arrangeGrob(
    p_cum +
      theme(
        legend.position = "inside",
        legend.position.inside = c(0.15, 0.85)
      ),
    gridExtra::arrangeGrob(p_hb, p_retic, nrow = 2, ncol = 1),
    nrow = 1, ncol = 2,
    widths = c(0.4, 0.6)
  )

  # Save plots to disk.
  ggsave(
    "individuals-exceeding-threshold.png",
    p_cmp,
    width = 6,
    height = 5
  )
  ggsave(
    "individuals-exceeding-threshold-num-optimal.png",
    p_num_optimal,
    width = 6,
    height = 5
  )
  ggsave(
    "individuals-exceeding-threshold-optimal-regimens.png",
    p_optimal,
    width = 6,
    height = 5
  )
  ggsave(
    "individuals-exceeding-threshold-optimal-vs-ascending-study.png",
    p_vs_study,
    width = 6,
    height = 9
  )
  ggsave(
    "individuals-exceeding-threshold-optimal-vs-ascending-study-wide.png",
    p_vs_study_wide,
    width = 9,
    height = 6
  )
  ggsave(
    "individuals-exceeding-threshold-optimal-predicted-results.png",
    grob_outputs,
    width = 10,
    height = 4
  )
}


# Returns a dataframe with columns:
# - threshold (g/dL)
# - fraction_below: fraction of individuals who experience a daily Hb drop
#   that exceeded `threshold`
# - num_optimal: the number of optimal regimens
# - first_optimal_ix: index of the first optimal regimen
# - duration: the regimen duration ("10 days" or "14 days")
load_or_collect_results <- function(results_file) {
  if (! file.exists(results_file)) {
    cat("Collecting results ...\n")

    result_tables <- lapply(
      settings$scenarios,
      function(scenario) {
        evaluate_regimens(
          scenario$chunk_regex,
          scenario$num_chunks,
          thresholds,
          scenario$matrix_file,
          scenario$results_file
        )
      }
    )

    results <- lapply(
      seq_along(plot_tables),
      function(ix) {
        plot_tables[[ix]]$plot |>
          mutate(duration = paste(settings$scenarios[[ix]]$duration, "days"))
      }
    ) |> bind_rows()

    # Save the summary data to disk.
    saveRDS(results, results_file, compress = "xz")
  } else {
    cat("Loading results from", results_file, "...\n")
    results <- readRDS(results_file)
  }

  results
}


fraction_below_threshold <- function(hb_drops, threshold) {
  rowSums(hb_drops < threshold) / ncol(hb_drops)
}


# Returns a matrix of dimensions N_REGIMENS x N_THRESHOLDS.
get_fractions_below_thresholds <- function(regex, num_chunks, thresholds) {
  chunk_files <- list.files(pattern = regex)
  if (length(chunk_files) != num_chunks) {
    stop("Expected ", num_chunks, " files but found ", length(chunk_files))
  }

  chunk_numbers <- as.integer(sub(regex, "\\1", chunk_files))

  chunked_results <- lapply(
    order(chunk_numbers),
    function(ix) {
      cat("Reading", chunk_files[ix], "...\n")
      chunk <- readRDS(chunk_files[ix])
      sapply(
        thresholds,
        function(threshold) {
          fraction_below_threshold(chunk, threshold)
        }
      )
    }
  )

  do.call(rbind, chunked_results)
}


# Returns a data frame of the optimal regimen(s) for each threshold, with
# columns 'regimen', 'threshold', 'fraction_below', 'num_optimal', and
# 'first_optimal_ix'.
fractions_data_frame <- function(fractions_mat, thresholds) {
  column_names <- paste0("threshold_", thresholds)
  colnames(fractions_mat) <- column_names

  as_tibble(fractions_mat) |>
    mutate(regimen = seq_len(nrow(fractions_mat))) |>
    pivot_longer(
      starts_with("threshold_"),
      names_prefix = "threshold_",
      names_to = "threshold",
      values_to = "fraction_below"
    ) |>
    mutate(
      threshold = as.numeric(threshold)
    ) |>
    group_by(threshold) |>
    mutate(regimen_ix = row_number()) |>
    filter(fraction_below == min(fraction_below)) |>
    # NOTE: record the # of optimal regimens for each threshold, and the index
    # of the first optimal regimen.
    mutate(
      num_optimal = n(),
      first_optimal_ix = dplyr::first(regimen_ix)
    ) |>
    ungroup()
}


evaluate_regimens <- function(regex, num_chunks, thresholds, matrix_file, results_file) {
  if (file.exists(matrix_file) && file.exists(results_file)) {
    results_tbl <- readRDS(results_file)
  } else {
    fracs_mat <- get_fractions_below_thresholds(regex, num_chunks, thresholds)
    saveRDS(fracs_mat, matrix_file, compress = "xz")

    results_tbl <- fractions_data_frame(fracs_mat, thresholds)
    saveRDS(results_tbl, results_file, compress = "xz")
  }

  plot_tbl <- results_tbl |>
    select(threshold, fraction_below, num_optimal, first_optimal_ix) |>
    unique() |>
    mutate(threshold = - threshold)

  list(results = results_tbl, plot = plot_tbl)
}


load_optimal_regimens <- function(settings, thresholds) {
  optimal_regimens <- list()

  for (scenario in settings$scenarios) {
    df_results <- readRDS(scenario$results_file) |>
      mutate(threshold = - threshold) |>
      filter(threshold %in% thresholds)

    all_regimens <- read_csv(
      scenario$regimens_file,
      show_col_types = FALSE
    ) |>
      mutate(regimen_ix = row_number())

    df_optimal <- df_results |>
      left_join(all_regimens, by = "regimen_ix") |>
      mutate(duration = paste(scenario$duration, "days"))

    optimal_regimens[[length(optimal_regimens) + 1]] <- df_optimal
  }

  df_wide <- bind_rows(optimal_regimens)

  # Convert to long format and calculate cumulative doses.
  dose_scale <- settings$dose_unit_mg / settings$weight_kg
  df_long <- df_wide |>
    pivot_longer(
      starts_with("day_"),
      names_prefix = "day_",
      names_to = "day",
      names_transform = as.integer,
      values_to = "dose"
    ) |>
    filter(dose > 0) |>
    # NOTE: we assume a 60kg body weight for calculating the cumulative doses
    # for each optimal regimen.
    group_by(regimen_ix, duration, threshold) |>
    mutate(
      cum_dosemgkg = cumsum(dose * dose_scale)
    ) |>
    ungroup()

  list(wide = df_wide, long = df_long)
}


load_ascending_study_regimens <- function() {
  data_env <- new.env()
  load(file.path("..", "Data", "RBC_model_data.RData"), envir = data_env)
  data_env$PQdat |>
    filter(study == "Part1") |>
    mutate(day = Study_Day + 1) |>
    select(ID, day, dose, dosemgkg, weight) |>
    mutate(ID_num = substring(ID, 6) |> as.integer()) |>
    filter(dosemgkg > 0) |>
    # Calculate the cumulative dose over each regimen.
    group_by(ID) |>
    mutate(
      cum_dosemgkg = cumsum(dosemgkg)
    ) |>
    ungroup()
}


predict_outcomes <- function(
  settings, df_regimens, num_draws = 1000, forward_sim = NULL, intervals = TRUE
) {
  calc <- new.env()
  sys.source("calculate-max-daily-Hb-drops.R", envir = calc)

  fit <- settings$load_fit()
  df_draws <- calc$get_draws(fit)
  data_inputs <- settings$get_data_inputs()
  forward_sim_params <- calc$build_forward_sim_args(df_draws, data_inputs)

  draw_ixs <- seq(
    from = 1,
    to = nrow(forward_sim_params),
    length.out = num_draws
  ) |>
    round()
  df_in <- forward_sim_params[draw_ixs, ]

  # Sample body weights from a normal distribution around a reference weight.
  prev_random_state <- .Random.seed
  set.seed(12345)
  body_weights <- rnorm(
    n = num_draws,
    mean = settings$weight_kg,
    sd = settings$weight_sd
  )
  .Random.seed <- prev_random_state

  df_regimens_wide <- df_regimens |>
    select(duration, regimen, day, dose) |>
    pivot_wider(
      id_cols = c(duration, regimen),
      names_from = day,
      names_prefix = "day_",
      values_from = dose
    ) |>
    replace_na()

  nComp_sim <- 29

  if (is.null(forward_sim)) {
    forward_sim <- settings$get_forward_sim_function()
  }

  results <- lapply(
    seq_len(nrow(df_regimens_wide)),
    function(regimen_ix) {
      regimen_units <- df_regimens_wide[regimen_ix, ] |>
        select(starts_with("day_")) |>
        as.numeric() |>
        replace_na(0)
      regimen_zeros <- 0 * seq_len(nComp_sim - length(regimen_units))
      regimen_units <- c(regimen_units, regimen_zeros)

      tbls <- list()

      for (draw_ix in seq_len(nrow(df_in))) {
        # Adjust the dose to account for variations in body weight.
        weight_kg <- body_weights[draw_ix]
        dose_scale <- settings$dose_unit_mg / weight_kg
        regimen_mg_kg <- dose_scale * regimen_units

        draw <- df_in[draw_ix, ] |> as.list()
        draw$dose_weights <- draw$dose_weights[[1]]
        draw$drug_regimen <- regimen_mg_kg

        # NOTE: use a baseline haemoglobin of 15 g/dL.
        draw$Hb_star <- 15

        out <- do.call(forward_sim, draw)

        tbl <- tibble(
          duration = df_regimens_wide$duration[regimen_ix],
          regimen = df_regimens_wide$regimen[regimen_ix],
          draw = rep(draw_ix, nComp_sim),
          day = seq(1, nComp_sim),
          hb_gdl = out[1, ],
          retic_pcnt = out[2, ],
          eff_dose = out[3, ]
        )

        tbls[[length(tbls) + 1]] <- tbl
      }

      bind_rows(tbls)
    }
  )

  tbl_results <- bind_rows(results) |>
    pivot_longer(c(hb_gdl, retic_pcnt, eff_dose), names_to = "measure") |>
    mutate(
      measure = case_when(
        measure == "hb_gdl" ~ "Haemoglobin (g/dL)",
        measure == "retic_pcnt" ~ "Reticulocyte (%)",
        measure == "eff_dose" ~ "Effective dose (mg/kg)"
      )
    )

  # Return mean, median, and 5%-95% intervals for each output.
  if (intervals) {
    tbl_results |>
      # NOTE: convert days from 1..29 to 0..28.
      mutate(day = day - 1) |>
      group_by(duration, regimen, day, measure) |>
      summarise(
        lower = quantile(value, 0.05),
        upper = quantile(value, 0.95),
        mean = mean(value),
        median = median(value),
        .groups = "drop"
      )
  } else {
    tbl_results |>
      # NOTE: convert days from 1..29 to 0..28.
      mutate(day = day - 1)
  }
}


plot_frac_exceeding_vs_threshold <- function(plot_tbl) {
  ggplot(plot_tbl) +
    geom_line(
      aes(threshold, 100 * fraction_below, colour = duration)
    ) +
    geom_point(
      aes(threshold, 100 * fraction_below, colour = duration)
    ) +
    scale_colour_brewer(
      name = "Duration",
      palette = "Dark2"
    ) +
    xlab("Hb drop threshold (g/dL)") +
    ylab("Individuals exceeding threshold (%)") +
    expand_limits(
      x = 0,
      y = c(0, 100)
    ) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.93, 0.93),
      legend.justification = c(1, 1)
    )
}


plot_number_of_optimal_regimens <- function(plot_tbl) {
  # NOTE: have to use factor(threshold) in order to adjust bar widths.
  ggplot(
    plot_tbl |>
      filter(threshold > 0.3, threshold < 1.3)
  ) +
    geom_col(
      aes(factor(threshold), num_optimal, fill = duration),
      position = position_dodge2(),
      width = 0.6
    ) +
    scale_x_discrete(
      breaks = as.character(seq(0.4, 1.4, by = 0.2))
    ) +
    scale_fill_brewer(
      name = "Duration",
      palette = "Dark2"
    ) +
    xlab("Hb drop threshold (g/dL)") +
    ylab("Number of optimal dose regimens") +
    theme(
      legend.position = "top"
    )
}


plot_optimal_regimens_vs_study <- function(df_optimal_long, ascending_study) {
  # Plot cumulative doses for the optimal regimens against all those of the
  # dose regimens from the ascending-dose study in a single pane.
  ggplot() +
    geom_line(
      aes(day, cum_dosemgkg, group = ID),
      ascending_study,
      colour = "#7f7f7f"
    ) +
    geom_line(
      aes(
        day, cum_dosemgkg,
        group = interaction(regimen_ix, duration, threshold),
        colour = duration
      ),
      df_optimal_long,
      linewidth = 1
    ) +
    geom_point(
      aes(day, cum_dosemgkg),
      ascending_study |>
        group_by(ID) |>
        filter(day == max(day)) |>
        ungroup()
    ) +
    geom_point(
      aes(
        day, cum_dosemgkg,
        colour = duration
      ),
      df_optimal_long |>
        group_by(regimen_ix, duration, threshold) |>
        filter(day == max(day)) |>
        ungroup()
    ) +
    scale_colour_brewer(
      name = "Duration",
      palette = "Dark2"
    ) +
    scale_x_continuous(
      "Days since start of primaquine",
      breaks = c(0, 7, 14, 21),
      minor_breaks = NULL,
      limits = c(0, 21)
    ) +
    scale_y_continuous(
      "Total primaquine dose (mg/kg)",
      minor_breaks = NULL,
      limits = c(0, NA)
    ) +
    theme_bw() +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.115, 0.875)
    )
}


plot_optimal_regimen_daily_doses <- function(settings, df_doses) {
  dose_scale <- settings$dose_unit_mg / settings$weight_kg
  ggplot(df_doses) +
    geom_step(
      aes(day, dose * dose_scale, colour = duration, group = regimen_ix)
    ) +
    scale_x_continuous(
      "Day",
      breaks = c(0, 5, 10, 15),
      limits = c(0, 15)
    ) +
    scale_y_continuous(
      "Daily dose (mg/kg)",
      limits = c(0, 0.8)
    ) +
    scale_colour_brewer(
      name = "Duration",
      palette = "Dark2"
    ) +
    facet_wrap(
      ~ threshold,
      labeller = labeller(threshold = ~ paste("<=", .x, "g/dL"))
    ) +
    theme(
      legend.position = "top"
    )
}


plot_optimal_output_intervals <- function(tbl_ints) {
  tbl_durn <- tbl_ints |>
    select(duration) |>
    unique() |>
    mutate(
      Start_Day = 0,
      Final_Day = as.numeric(substring(duration, 1, 2))
    )

  # Use the optimal dose regimen colours for ribbons.
  ribbon_colours <- scales::brewer_pal(palette = "Dark2")(2)

  # Use darker shades for the mean predictions.
  line_colours <- c(
    scales::brewer_pal(palette = "Greens")(9)[9],
    scales::brewer_pal(palette = "Oranges")(9)[9]
  )

  ggplot(tbl_ints) +
    geom_rect(
      aes(xmin = Start_Day, xmax = Final_Day, ymin = -Inf, ymax = Inf),
      tbl_durn,
      fill = "#9f9f9f",
      alpha = 0.3
    ) +
    geom_ribbon(
      aes(day, ymin = lower, ymax = upper, fill = duration),
    ) +
    geom_vline(
      aes(xintercept = Start_Day),
      tbl_durn,
      linetype = "dashed"
    ) +
    geom_vline(
      aes(xintercept = Final_Day),
      tbl_durn,
      linetype = "dashed"
    ) +
    geom_line(
      aes(day, mean, colour = duration),
    ) +
    scale_x_continuous(
      "Days since start of primaquine",
      breaks = 7 * (0:4)
    ) +
    scale_colour_manual(
      values = line_colours,
      guide = "none"
    ) +
    scale_fill_manual(
      values = ribbon_colours,
      guide = "none"
    ) +
    ylab(NULL) +
    facet_wrap(~ duration) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      panel.grid.minor = element_blank()
    )
}


plot_optimal_daily_doses_vs_ascending <- function(
  settings, df_doses, ascending_study
) {
  dose_scale <- settings$dose_unit_mg / settings$weight_kg
  ggplot() +
    geom_step(
      aes(
        day, dose * dose_scale,
        group = interaction(regimen_ix, duration),
        colour = duration
      ),
      df_doses
    ) +
    geom_step(
      aes(day, dosemgkg, group = ID),
      ascending_study
    ) +
    scale_colour_brewer(
      name = "Duration",
      palette = "Dark2"
    ) +
    scale_x_continuous(
      "Day",
      breaks = c(1, 7, 14, 21),
      minor_breaks = 0:21,
      limits = c(0, 21)
    ) +
    scale_y_continuous(
      "Dose (mg/kg)"
    ) +
    facet_wrap(
      ~ ID_num,
      ncol = 4,
      labeller = as_labeller(function(x) paste("ADPQ", x))
    ) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.97, 0.01),
      legend.justification = c(1, 0)
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

call_main("identify-optimal-regimens.R")
