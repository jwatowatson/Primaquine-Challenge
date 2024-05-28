#!/usr/bin/env -S Rscript --vanilla
#
# Calculate the number of possible drug regimens that satisfy the constraints.
#

weight_kg <- 60
net_dose_mg_kg <- 5
dose_unit_mg <- 2.5
max_daily_dose_mg <- 45

max_daily_dose_units <- floor(max_daily_dose_mg / dose_unit_mg)
total_doses <- floor(weight_kg * net_dose_mg_kg / dose_unit_mg)

# How many ways can we distribute `total_doses` over `num_days` such that:
#
# - Each daily dose comprises 1..`max_daily_dose_units` units; and
#
# - The daily dose increases monotonically.
#
# The following uses a recursive approach and is rather slow for 14 days.

count_combinations <- function(day = 1, min_units = 1, max_units = 18,
                               given_doses = 0) {
  combinations <- list()
  for (unit in min_units:max_units) {
    # Create an vector to contain the daily doses.
    regimen <- 0 * (0:(num_days - day))
    regimen[1] <- unit
    new_given_doses <- given_doses + unit
    if (new_given_doses > total_doses) {
      break
    }
    if (day < num_days) {
      subsequent_combinations <- count_combinations(
        day = day + 1, min_units = unit, max_units = max_units,
        given_doses = new_given_doses
      )
      for (ix in seq_along(subsequent_combinations)) {
        regimen[-1] <- subsequent_combinations[[ix]]
        if (total_doses == sum(regimen) + given_doses) {
          combinations[[length(combinations) + 1]] <- regimen
        }
      }
    } else {
      if (total_doses == sum(regimen) + given_doses) {
        combinations[[length(combinations) + 1]] <- regimen
      }
    }
  }
  combinations
}

# NOTE: 3,648,057 possible regimens.
num_days <- 14
combs_14_days <- count_combinations()

# NOTE: 78,796 possible regimens.
num_days <- 10
combs_10_days <- count_combinations()

# Convert these regimens to data frames and save them to CSV files.
library(tibble)
library(readr)

df_10_days <- do.call(rbind, combs_10_days) |>
  as_tibble(.name_repair = function(x) paste0("day_", seq_along(x)))
write_csv(df_10_days, "combinations-10-days.csv")

df_14_days <- do.call(rbind, combs_14_days) |>
  as_tibble(.name_repair = function(x) paste0("day_", seq_along(x)))
write_csv(df_14_days, "combinations-14-days.csv")
