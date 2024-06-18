#
# Values and functions that are used by multiple scripts.
#

# The net dose is 5 mg/kg.
net_dose_mg_kg <- 5

# The dose unit is 2.5 mg (half of a 5 mg tablet).
dose_unit_mg <- 2.5

# The maximum daily dose is 45 mg (9 x 5 mg tablets).
max_daily_dose_mg <- 45

# For simplicity we assume a body weight of 60kg, so that each dose regimen
# comprises 300 mg, which can be divided into 120 x 2.5mg doses.
weight_kg <- 60

# For sampling body weights for each simulated individual, we assume a normal
# distribution with mean `weight_kg` and the following standard deviation.
weight_sd <- 5

# For 10-day and 14-day dose regimens, define the various file names and how
# many chunks into which the maximum daily drop calculations are divided.
scenarios <- list(
  list(
    duration = 10,
    num_chunks = 1,
    regimens_file = "combinations-10-days.csv",
    evaluation_file_prefix = "evaluation-10-days-part-",
    chunk_regex = "^evaluation-10-days-part-([0-9]+)\\.rds$",
    matrix_file = "optimal-10-day-regimens-matrix.rds",
    results_file = "optimal-10-day-regimens-table.rds"
  ),
  list(
    duration = 14,
    num_chunks = 20,
    regimens_file = "combinations-14-days.csv",
    evaluation_file_prefix = "evaluation-14-days-part-",
    chunk_regex = "^evaluation-14-days-part-([0-9]+)\\.rds$",
    matrix_file = "optimal-14-day-regimens-matrix.rds",
    results_file = "optimal-14-day-regimens-table.rds"
  )
)


get_stan_model_functions <- function() {
  utils <- new.env()
  sys.source(file.path("..", "cmdstan_utils.R"), envir = utils)
  utils$load_packages(plot_libs = TRUE)

  model_file <- file.path(
    "..",
    "Stan_models", "RBC_model_master_pop_free_weights_cmdstan.stan"
  )
  model <- utils$compile_model_with_exposed_functions(model_file)

  model$functions
}


get_forward_sim_function <- function() {
  functions <- get_stan_model_functions()
  functions$forwardsim
}


get_effective_dose_function <- function() {
  functions <- get_stan_model_functions()
  functions$compute_effective_dose
}


load_fit <- function() {
  fit_file <- file.path(
    "..",
    "Rout",
    "pop_fit_free_weights_cmdstan_max_delay_9_job3.rds"
  )
  readRDS(fit_file)
}


get_data_inputs <- function() {
  # NOTE: need to adjust the working directory for `create_job_data()`.
  curr_wd <- getwd()
  setwd("..")

  utils <- new.env()
  sys.source(file.path("cmdstan_utils.R"), envir = utils)
  utils$load_packages(plot_libs = TRUE)

  # NOTE: the following arguments are data:
  data_args <- c(
    "T_nmblast", "T_retic", "T_RBC_max", "T_transit_steady_state", "K_weights",
    "sigma"
  )
  job_data <- utils$create_job_data(
    job_number = 3,
    max_dose_delay = 9,
    quiet = TRUE
  )
  args <- job_data[data_args]
  args$nComp_sim <- 29

  # NOTE: restore the original working directory.
  setwd(curr_wd)

  args
}
