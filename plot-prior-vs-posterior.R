#!/usr/bin/env -S Rscript --vanilla

utils <- new.env()
sys.source("cmdstan_utils.R", envir = utils)
utils$load_packages(plot_libs = TRUE)

fit_file <- file.path(
  "Rout", "pop_fit_free_weights_cmdstan_max_delay_9_job3.rds"
)
fit <- readRDS(fit_file)

re_match <- "pop_fit_free_weights_cmdstan_max_delay_9_job(.*)\\.rds"

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
)

# Replace the "logit_alpha" values with "alpha" values (%).
all_draws <- bind_rows(
  draws |>
    filter(name != "logit_alpha"),
  draws |>
    filter(name == "logit_alpha") |>
    mutate(
      value = 100 * exp(value) / (1 + exp(value)),
      name = "alpha"
    )
) |>
  arrange(.draw, job, name)

# Create descriptive labels for each facet.
#
# NOTE: use plotmath expressions and `label_parsed` to display the
# descriptive label and the notation used in the manuscript.
labels <- c(
  "alpha" = "paste('Max RBC lifespan reduction (', alpha, ', %)')",
  "beta" = "paste('Half-maximal effect dose (', beta, ', mg/kg)')",
  "h" = "paste('Dose-response slope (', h, ')')",
  "Hb_star" = "paste('Steady-state haemoglobin (', 'Hb'*'*', ', g/dl)')",
  "T_E_star" = "paste('Steady-state RBC lifespan (', T[E]*'*', ', days)')"
)

labelled_draws <- all_draws |>
  mutate(
    label = factor(
      name,
      levels = names(labels),
      labels = unname(labels),
      ordered = TRUE
    )
  )

inv_logit <- function(x) {
  exp(x) / (1 + exp(x))
}

tbl_prior <- tibble(
  alpha = 100 * inv_logit(rnorm(n = 1000, mean = 0, sd = 2)),
  beta = rnorm(n = 1000, mean = 0.75, sd = 0.2),
  h = rexp(n = 1000, rate = 1),
  Hb_star =  rnorm(n = 1000, mean = 15, sd = 1),
  T_E_star = rnorm(n = 1000, mean = 90, sd = 20)
) |>
  pivot_longer(everything()) |>
  mutate(
    label = factor(
      name,
      levels = names(labels),
      labels = unname(labels),
      ordered = TRUE
    )
  )

p_dists <- ggplot() +
  geom_density(
    aes(value, colour = job),
    labelled_draws
  ) +
  geom_density(
    aes(value),
    tbl_prior,
    colour = "black",
    linetype = "dashed",
    linewidth = 1
  ) +
  scale_colour_brewer(
    NULL,
    palette = "Set2"
  ) +
  xlab(NULL) +
  ylab("Density") +
  facet_wrap(~ label, scales = "free", labeller = label_parsed) +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.95, 0.05),
    legend.justification = c(1, 0),
    strip.background = element_blank()
  )

plot_file <- "prior_vs_posterior.png"
cat("Writing", plot_file, "...")
png(plot_file, width = 8, height = 4, unit = "in", res = 150)
print(p_dists)
invisible(dev.off())
cat(" done\n")
