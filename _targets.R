# _targets.R

# setup -------------------------------------------------------------------

library(targets)
library(tarchetypes)

# set target options
tar_option_set(
  packages = c("tidyverse"),
  format = "qs"
)

# source R scripts
tar_source()

# parallel processing
options(clustermq.scheduler = "multicore")
future::plan(future.callr::callr(workers = future::availableCores() - 1))


# targets -----------------------------------------------------------------

list(
  tar_target(
    name = data,
    command = tibble(x = rnorm(100), y = rnorm(100))
    #   format = "feather" # efficient storage of large data frames # nolint
  ),
  tar_target(
    name = model,
    command = coefficients(lm(y ~ x, data = data))
  )
)
