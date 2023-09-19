#!/usr/bin/env -S Rscript --vanilla
#
# Plot the results obtained by fitting the RBC model with independent
# dose-weighting parameters.
#


main <- function() {
  library(readr)
  library(dplyr, warn.conflicts = FALSE)
  library(ggplot2)

  results_file <- "Rout/pop_fit_free_weights.RData"
  load(results_file)

  # Extract the outputs from forward_sim().
  Y_hat1 <- rstan::extract(mod_fit_pop_free_wt, pars = "Y_hat")$Y_hat

  # Produce a density plot of the dose-weighting parameters.
  png("test-free-weights-densities.png",
      width = 8, height = 10, units = "in", res = 150)
  rstan::plot(
    mod_fit_pop_free_wt,
    pars = "dose_weights",
    show_density = TRUE
  )
  invisible(dev.off())

  # Plot model fits against the individual data sets.
  png("test-free-weights-fit-hb.png",
      width = 8, height = 20, units = "in", res = 150)
  plot_haemoglobin_fits(stan_data, Y_hat1, mfrow = c(10, 4))
  invisible(dev.off())

  png("test-free-weights-fit-rc.png",
      width = 8, height = 20, units = "in", res = 150)
  plot_reticulocyte_fits(stan_data, Y_hat1, mfrow = c(10, 4))
  invisible(dev.off())

  # Plot the individual dose response curves.
  png("test-free-weights-dose-response.png",
      width = 8, height = 8, units = "in", res = 150)
  plot_dose_response_curves(stan_data, mod_fit_pop_free_wt)
  invisible(dev.off())
}


#' Plot individual fits for haemoglobin counts.
#'
#' This is adapted from `plot_fits.R`.
#'
#' @param stan_data The input data provided to the model.
#' @param Y_hat1 The `forwardsim()` outputs: haemoglobin and reticulocyte
#'        percentage.
#' @param mfrow The number of columns and rows by which to arrange plots.
plot_haemoglobin_fits <- function(stan_data, Y_hat1, mfrow = c(4, 10)) {
  par(mfrow = mfrow, las = 1, family = "serif",
      cex.lab = 1.3, cex.axis = 1.3, mar = c(2, 2, 2, 2))

  for (id in 1:stan_data$N_experiment) {
    ind_id <- stan_data$ind_start_regimen[id]:stan_data$ind_end_regimen[id]
    Y_hat_hb <- Y_hat1[, 1, ind_id]

    ind_haemocue <- stan_data$t_sim_hemocue_Hb[
      stan_data$ind_start_Hb_hemocue[id]:stan_data$ind_end_Hb_hemocue[id]
    ]
    ind_hb_CBC <- stan_data$t_sim_CBC_Hb[
      stan_data$ind_start_Hb_CBC[id]:stan_data$ind_end_Hb_CBC[id]
    ]

    plot(
      seq_along(ind_id),
      colMeans(Y_hat_hb),
      type = "l",
      lwd = 2,
      panel.first = grid(),
      ylim = c(8.3, 16.3),
      xlab="",
      ylab="",
      # TODO: cannot find where `ID_map` is defined.
      # main = ID_map$ID_combined[id],
      col = "grey"
    )

    dosing_ind <- which(stan_data$drug_regimen[ind_id] > 0)
    dosing_ind <- c(head(dosing_ind, 1), tail(dosing_ind, 1) + 1)
    polygon(
      c(dosing_ind, rev(dosing_ind)),
      c(0, 0, 20 ,20),
      col = adjustcolor("pink", 0.3),
      border = NA
    )
    polygon(
      c(seq_along(ind_id), rev(seq_along(ind_id))),
      c(apply(Y_hat_hb, 2, quantile, probs = 0.1),
        rev(apply(Y_hat_hb, 2, quantile, probs = 0.9))),
      col = adjustcolor("grey", 0.3),
      border = NA
    )
    lines(seq_along(ind_id), colMeans(Y_hat_hb), lwd = 2, col = "grey")
    points(
      seq_along(ind_id)[ind_haemocue],
      stan_data$Hb_Haemocue[
        stan_data$ind_start_Hb_hemocue[id]:stan_data$ind_end_Hb_hemocue[id]
      ],
      pch = 16,
      cex = 1
    )
    points(
      seq_along(ind_id)[ind_hb_CBC],
      stan_data$Hb_CBC[
        stan_data$ind_start_Hb_CBC[id]:stan_data$ind_end_Hb_CBC[id]
      ],
      pch = 16,
      col = "blue",
      cex = 1
    )
  }
}


#' Plot individual fits for reticulocyte counts.
#'
#' This is adapted from `plot_fits.R`.
#'
#' @param stan_data The input data provided to the model.
#' @param Y_hat1 The `forwardsim()` outputs: haemoglobin and reticulocyte
#'        percentage.
#' @param mfrow The number of columns and rows by which to arrange plots.
plot_reticulocyte_fits <- function(stan_data, Y_hat1, mfrow = c(4, 10)) {
  par(mfrow = mfrow, las = 1, family = "serif",
      cex.lab = 1.3, cex.axis = 1.3, mar = c(2, 2, 2, 2))

  for(id in 1:stan_data$N_experiment){
    ind_id <- stan_data$ind_start_regimen[id]:stan_data$ind_end_regimen[id]
    Y_hat_retic <- Y_hat1[, 2, ind_id]

    ind_retic <- stan_data$t_sim_retic[
      stan_data$ind_start_retic[id]:stan_data$ind_end_retic[id]
    ]

    plot(
      seq_along(ind_id),
      colMeans(Y_hat_retic),
      type = "l",
      lwd = 2,
      panel.first = grid(),
      ylim = c(0, 20),
      xlab="",
      ylab="",
      # TODO: cannot find where `ID_map` is defined.
      # main = ID_map$ID_combined[id],
      col = "grey"
    )

    dosing_ind <- which(stan_data$drug_regimen[ind_id] > 0)
    dosing_ind <- c(head(dosing_ind, 1), tail(dosing_ind, 1) + 1)
    polygon(
      c(dosing_ind, rev(dosing_ind)),
      c(0, 0, 20, 20),
      col = adjustcolor("pink", 0.3),
      border = NA
    )
    polygon(
      c(seq_along(ind_id), rev(seq_along(ind_id))),
      c(apply(Y_hat_retic, 2, quantile, probs = 0.1),
        rev(apply(Y_hat_retic, 2, quantile, probs = 0.9))),
      col = adjustcolor("grey", 0.3),
      border = NA
    )
    lines(seq_along(ind_id), colMeans(Y_hat_retic), lwd = 2, col = "grey")
    points(
      seq_along(ind_id)[ind_retic],
      stan_data$Retic_data[
        stan_data$ind_start_retic[id]:stan_data$ind_end_retic[id]
      ],
      pch = 16,
      cex = 1
    )
  }
}

#' Plot individual dose-response curves.
#'
#' This is adapted from `plot_fits.R`.
#'
#' @param stan_data The input data provided to the model.
#' @param stan_fit The output from fitting the RBC model with independent
#'        dose-weighting parameters.
plot_dose_response_curves <- function(stan_data, stan_fit) {
  master <- new.env(parent = baseenv())
  sys.source("master_functions.R", envir = master)

  pars <- c(
    "Hb_star",
    "T_E_star",
    "alpha_diff1",
    "alpha_delta1",
    "alpha_diff2",
    "alpha_delta2",
    "logit_alpha",
    "beta",
    "h",
    "mean_delay",
    "theta_rand"
  )
  my_thetas1 <- rstan::extract(stan_fit, pars = pars)

  par(mfrow = c(1, 1), las = 1, family = "serif",
      cex.lab = 1.3, cex.axis = 1.3)

  xs <- seq(0, 1, length.out = 200)
  effect <- array(dim = c(length(my_thetas1$beta), length(xs)))

  for (i in seq_along(xs)) {
    effect[, i] <- master$emax_f(
      x = xs[i],
      logit_max_drug_effect = my_thetas1$logit_alpha,
      emax_k = my_thetas1$h,
      ed50 = my_thetas1$beta
    )
  }

  effect_mean <- colMeans(effect)

  plot(
    xs,
    100 * effect_mean,
    type = "l",
    lwd = 3,
    xlab = "Primaquine mg/kg",
    ylab = "Reduction in red cell lifespan (%)",
    panel.first = grid(),
    ylim = c(0, 60)
  )
  polygon(
    c(xs, rev(xs)),
    100 * c(apply(effect, 2, quantile, probs = 0.1),
            rev(apply(effect, 2, quantile, probs = 0.9))),
    col = adjustcolor("darkblue", 0.3),
    border = NA
  )

  effect_id_list <- list()
  for (id in seq_len(stan_data$N_subject)) {
    # NOTE: the random effects differ from those in `plot_fits.R`.
    logit_max_drug_effect <- (
      my_thetas1$logit_alpha + my_thetas1$theta_rand[, id, 6]
    )
    emax_k <- my_thetas1$h
    ed50 <- my_thetas1$beta * exp(my_thetas1$theta_rand[, id, 7])

    effect_id <- array(dim = c(length(my_thetas1$beta), length(xs)))
    for (i in seq_along(xs)) {
      effect_id[, i] <- master$emax_f(
        x = xs[i],
        logit_max_drug_effect = logit_max_drug_effect,
        emax_k = emax_k,
        ed50 = ed50
      )
    }

    effect_id_list[[id]] <- colMeans(effect_id)
    lines(xs, 100 * colMeans(effect_id), col = "black")
  }

  lines(xs, 100 * colMeans(effect), lwd = 3, col = "red")
}


main()
