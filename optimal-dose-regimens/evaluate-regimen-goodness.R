#!/usr/bin/env -S Rscript --vanilla

main <- function() {
  library(dplyr, warn.conflicts = FALSE)
  library(ggplot2)
  library(readr)
  library(tibble)
  library(tidyr)

  # Retrieve dose regimen details from `settings.R`.
  settings <- new.env()
  sys.source("settings.R", envir = settings)

  # Identify the matrix column that contains the results for our chosen
  # threshold of 1.0 g/dL.
  chosen_threshold <- 1.0
  dfs <- optimal_for_threshold(chosen_threshold)

  # Calculate daily intervals for the cumulative dose.
  df_intervals <- cumulative_dose_intervals(settings, dfs$all, dfs$best)

  # Plot daily cumulative dose intervals for near-optimal regimens.
  p_intervals <- ggplot() +
    geom_ribbon(
      aes(
        day,
        ymin = cum_dosemgkg_min,
        ymax = cum_dosemgkg_max,
        fill = interval,
        colour = interval,
        group = interaction(interval, duration)
      ),
      df_intervals
    ) +
    scale_fill_brewer("Exceeding 1 g/dL drop:") +
    scale_colour_brewer("Exceeding 1 g/dL drop:") +
    scale_x_continuous(breaks = seq(0, 14, by = 2)) +
    expand_limits(x = 0) +
    xlab("Days since start of primaquine") +
    ylab("Total primaquine dose (mg/kg)") +
    facet_wrap(~ duration, scale = "free_x", nrow = 1) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )

  ggsave(
    "regimens-1gdl-drop-intervals.png", p_intervals, width = 6, height = 5
  )

  # Plot a histogram of dose regimen performance.
  p_hist <- ggplot() +
    geom_histogram(
      aes(pcnt_exceeding, after_stat(density), fill = duration),
      dfs$all,
      binwidth = 1,
      position = "dodge"
    ) +
    geom_vline(
      aes(xintercept = pcnt_exceeding),
      dfs$best,
      linetype = "dashed"
    ) +
    geom_text(
      aes(
        x = pcnt_exceeding,
        y = 0.04,
        label = paste0(pcnt_exceeding, "%")
      ),
      dfs$best,
      hjust = 0,
      nudge_x = 1
    ) +
    scale_fill_brewer(
      "Duration",
      palette = "Dark2"
    ) +
    xlab("Individuals exceeding 1 g/dL drop (%)") +
    ylab(NULL) +
    # facet_wrap(~ duration, ncol = 1) +
    facet_wrap(~ duration, nrow = 1) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.1, 0.89)
    )

  ggsave(
    "regimens-1gdl-drop-histogram.png", p_hist, width = 6, height = 5
  )
}


optimal_for_threshold <- function(chosen_threshold) {
  # Retrieve dose regimen details from `settings.R`.
  settings <- new.env()
  sys.source("settings.R", envir = settings)

  thresholds <- seq(from = 0.1, to = 2, by = 0.1)
  abs_diffs <- abs(thresholds - abs(chosen_threshold))
  threshold_matrix_column <- which(abs_diffs < 1e-3)
  if (length(threshold_matrix_column) != 1) {
    print(chosen_threshold)
    print(thresholds)
    print(abs_diffs)
    stop("Found ", length(threshold_matrix_column), " matches for threshold")
  }

  # Retrieve the fraction of individuals whose greatest daily drop in Hb
  # exceeded our chosen threshold, for each dose regimen.
  df_all <- bind_rows(lapply(
    settings$scenarios,
    function(scenario) {
      scenario_col <- readRDS(scenario$matrix_file)[, threshold_matrix_column]
      tibble(
        duration = paste(scenario$duration, "days"),
        regimen_ix = seq_along(scenario_col),
        frac_exceeding = scenario_col,
        pcnt_exceeding = 100 * scenario_col
      )
    }
  ))

  # Identify the optimal 10-day and 14-day dose regimens.
  df_best <- df_all |>
    group_by(duration) |>
    filter(frac_exceeding == min(frac_exceeding))

  list(all = df_all, best = df_best)
}


load_regimens <- function(scenarios, regimens) {
  bind_rows(lapply(
    scenarios,
    function(scenario) {
      all_regimens <- read_csv(
        scenario$regimens_file,
        show_col_types = FALSE
      ) |>
        mutate(
          regimen_ix = row_number(),
          duration = paste(scenario$duration, "days")
        ) |>
        right_join(
          regimens |>
            select(duration, regimen_ix),
          by = c("duration", "regimen_ix")
        )
    }
  )) |>
    pivot_longer(
      starts_with("day_"),
      names_prefix = "day_",
      names_to = "day",
      names_transform = as.integer,
      values_to = "dose"
    ) |>
    filter(dose > 0)
}


cumulative_dose_intervals <- function(settings, df_all, df_best) {
  interval_spreads_pcnt <- c(1, 2, 5, 10, 20)
  interval_dfs <- list()
  for (spread_pcnt in interval_spreads_pcnt) {
    df_good_enough <- df_all |>
      inner_join(
        df_best |>
          rename(best_pcnt = pcnt_exceeding) |>
          select(duration, best_pcnt),
        by = "duration"
      ) |>
      filter(pcnt_exceeding - best_pcnt <= spread_pcnt)

    df_regimens <- load_regimens(settings$scenarios, df_good_enough)

    scale_mgkg <- settings$dose_unit_mg / settings$weight_kg
    df_interval <- df_regimens |>
      group_by(regimen_ix, duration) |>
      mutate(
        cum_dosemgkg = cumsum(dose * scale_mgkg)
      ) |>
      ungroup() |>
      group_by(duration, day) |>
      summarise(
        cum_dosemgkg_min = min(cum_dosemgkg),
        cum_dosemgkg_max = max(cum_dosemgkg),
        .groups = "drop"
      ) |>
      mutate(
        spread = paste0("+", spread_pcnt, "%")
      )

    interval_dfs[[length(interval_dfs) + 1]] <- df_interval
  }

  df_intervals <- bind_rows(interval_dfs) |>
    mutate(
      interval = factor(
        spread,
        levels = paste0("+", rev(interval_spreads_pcnt), "%"),
        ordered = TRUE
      )
    )
}


main()
