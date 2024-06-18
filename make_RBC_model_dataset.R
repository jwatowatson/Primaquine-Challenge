#!/usr/bin/env -S Rscript --vanilla
#
# Extract the study data, and adjust the initial Hb measurements
# for individual "ADPQ 2" in the ascending dose study.
#

suppressPackageStartupMessages(library(tidyverse))

load_data <- function(fix_ADPQ2 = TRUE) {
  load("Data/PQdata.rds")
  PQdat$Day <- as.numeric(PQdat$Day)

  # Individual "ADPQ 2" has clearly incorrect Hb measurements on day 0.
  # These should be replaced by the pre-study values.
  ADPQ2_initial_hb <- PQdat |>
    filter(label == "ADPQ 2", Day < 0) |>
    group_by(label) |>
    summarise(
      Haemocue_hb = mean(Haemocue_hb, na.rm = TRUE),
      CBC_hb = mean(CBC_hb, na.rm = TRUE)
    ) |>
    select(! label)

  cols <- c(
    "ID", "ID2", "study", "Study_Day", "G6PD_variant", "Haemocue_hb",
    "CBC_hb", "CBC_retic", "Manual_retic", "Mean_hb", "Mean_retic",
    "dose", "dosemgkg", "weight"
  )

  study_data <- PQdat |>
    filter(
      # Subject 18 is excluded from the analysis.
      label != "ADPQ 18",
      Day >= 0,
      Day <= 30
    ) |>
    group_by(label, Day) |>
    mutate(
      Study_Day = unique(Day),
      Haemocue_hb = mean(Haemocue_hb, na.rm = TRUE),
      CBC_hb = mean(CBC_hb, na.rm = TRUE),
      CBC_retic = mean(CBC_retic, na.rm = TRUE),
      Manual_retic = mean(Manual_retic, na.rm = TRUE),
      dose = sum(PQdose, na.rm = TRUE)
    ) |>
    ungroup() |>
    distinct(
      label,
      Day,
      .keep_all = TRUE
    ) |>
    arrange(label, Day) |>
    group_by(label, Day) |>
    mutate(
      dosemgkg = dose / weight,
      Mean_retic = mean(c(Manual_retic, CBC_retic), na.rm = TRUE),
      Mean_hb = mean(c(Haemocue_hb, CBC_hb), na.rm = TRUE)
    ) |>
    ungroup() |>
    mutate(
      ID = ifelse(! is.na(Ascending_ID), Ascending_ID, label),
      ID2 = label
    ) |>
    mutate_all(function(x) ifelse(is.nan(x), NA, x)) |>
    select(all_of(cols))

  if (fix_ADPQ2) {
    mask <- study_data$ID == "ADPQ 2" & study_data$Study_Day == 0

    for (column in names(ADPQ2_initial_hb)) {
      study_data[[column]][mask] <- ADPQ2_initial_hb[[column]]
    }
  }

  writeLines("Unique individuals:")
  print(unique(sort(study_data$ID)))

  writeLines("Unique IDs with repeats:")
  print(unique(sort(study_data$ID2)))

  output_file <- file.path("Data", "RBC_model_data.RData")
  cat("Saving to", output_file, "...\n")
  PQdat <- study_data
  save(PQdat, file = output_file)

  study_data
}

PQdat <- load_data(fix_ADPQ2 = TRUE)
