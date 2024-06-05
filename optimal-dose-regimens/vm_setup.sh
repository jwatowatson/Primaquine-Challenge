#!/bin/sh

sudo apt install git screen rsync r-cran-tidyverse r-cran-furrr r-cran-posterior r-cran-bayesplot

R --vanilla --quiet << END
if (! require("cmdstanr", quietly = TRUE)) {
  install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", "https://mc-stan.org/r-packages/", getOption("repos")))
  library(cmdstanr)
  install_cmdstan()
}
END
