if(! "nfidd" %in% rownames(installed.packages())) {
  pak::pak("nfidd/nfidd", dependencies = "all", upgrade = TRUE)
}

library("nfidd")
library("ggplot2")
library("dplyr")
library("tidyr")
library("lubridate")
library("here")
library("cmdstanr")

set.seed(123)
options(cmdstanr_print_line_numbers = TRUE)
