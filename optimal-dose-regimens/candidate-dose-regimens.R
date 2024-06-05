#!/usr/bin/env -S Rscript --vanilla
#
# Identify all candidate drug regimens that satisfy the constraints.
#

main <- function() {
  library(tibble)
  library(readr)

  # Retrieve dose regimen details from `settings.R`.
  settings <- new.env()
  sys.source("settings.R", envir = settings)

  max_daily_dose_units <- floor(
    settings$max_daily_dose_mg / settings$dose_unit_mg
  )
  total_doses <- floor(
    settings$weight_kg * settings$net_dose_mg_kg / settings$dose_unit_mg
  )

  for (scenario in settings$scenarios) {
    num_days <- scenario$duration
    output_file <- scenario$regimens_file

    cat("Collecting dose regimens for", num_days, "days ...")
    regimens <- dose_regimens(
      num_days,
      max_units = max_daily_dose_units,
      total_doses = total_doses
    )

    regimens_tbl <- do.call(rbind, regimens) |>
      as_tibble(.name_repair = function(x) paste0("day_", seq_along(x)))
    cat(" done\n")
    cat("Found", nrow(regimens_tbl), "candidate regimens\n")

    cat("Writing", output_file, "...")
    write_csv(regimens_tbl, output_file)
    cat(" done\n")
  }
}


# How many ways can we distribute `total_doses` over `num_days` such that:
#
# - Each daily dose comprises 1..`max_daily_dose_units` units; and
#
# - The daily dose increases monotonically.
#
# The following uses a recursive approach and is rather slow for 14 days.
#
dose_regimens <- function(num_days, day = 1, min_units = 1, max_units = 18, given_doses = 0, total_doses = 0) {
  combinations <- list()
  for (unit in min_units:max_units) {
    # Create an vector to contain the daily doses.
    regimen <- 0 * (0:(num_days - day))
    regimen[1] <- unit
    new_given_doses <- given_doses + unit

    # Check if we've exceeded the total dose.
    if (new_given_doses > total_doses) {
      break
    }

    if (day < num_days) {
      # This is not the final dose, so generate doses for the remaining days.
      subsequent_combinations <- dose_regimens(
        num_days,
        day = day + 1, min_units = unit, max_units = max_units,
        given_doses = new_given_doses,
        total_doses = total_doses
      )

      # Only retain regimens that deliver the total dose.
      for (ix in seq_along(subsequent_combinations)) {
        regimen[-1] <- subsequent_combinations[[ix]]
        if (total_doses == sum(regimen) + given_doses) {
          combinations[[length(combinations) + 1]] <- regimen
        }
      }
    } else {
      # This is the final dose, only retain the regimen if it delivers the
      # total dose.
      if (total_doses == sum(regimen) + given_doses) {
        combinations[[length(combinations) + 1]] <- regimen
      }
    }
  }
  combinations
}


main()
