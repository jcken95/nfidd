# packages ----

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

# useful for my excersise solns
library("patchwork")
library("purrr")
library("tidybayes")

# environment ----

set.seed(123)
options(cmdstanr_print_line_numbers = TRUE)

# helper functions ----


fit_model <- function(
    stan_model,
    true_params = tibble("meanlog" = 1.75, "sdlog" = 0.5)
) {
  # model outline
  # y = time to hospitalisation after symptom first present
  # likelihood:
  # y ~ lognormal(mu, sig)
  # prior:
  # mu ~ normal(0, 10)
  # sig ~ half-normal(0, 10)

  data(infection_times)

  true_params_long <- pivot_longer(true_params, cols = everything(), names_to = "param", values_to = "value")

  # adapted from nfidd/snippets/onset-hosp.r
  df <- infection_times |>
    mutate(
      # Gamma distribution for onset times
      onset_time = infection_time + rgamma(n(), shape = 5, rate = 1),
      # Lognormal distribution for hospitalisation times
      hosp_time = onset_time + rlnorm(n(),
                                      meanlog = true_params$meanlog,
                                      sdlog = true_params$sdlog)
    ) |>
    mutate(
      hosp_time = if_else(
        # use the binomial distribution for random binary outcomes
        rbinom(n = n(), size = 1, prob = 0.3) == 1,
        hosp_time,
        NA_real_
      )
    )

  df_dates <- df |>
    mutate(across(everything(), floor)) |>
    mutate(
      incubation_period = onset_time - infection_time,
      onset_to_hosp = hosp_time - onset_time
    )


  stan_data <- drop_na(df_dates, onset_to_hosp)
  ## Use the data to sample from the model posterior
  res <- stan_model$sample(
    data = list(
      n = nrow(stan_data),
      y = stan_data$onset_to_hosp
    )
  )

  posterior_samples <- res$draws() |>
    tidybayes::tidy_draws() |>
    select(-lp__) |>
    pivot_longer(
      cols = c(meanlog, sdlog),
      values_to = "value",
      names_to = "param"
    ) |>
    left_join(true_params_long,
              suffix = c("", "_true"),
              by = c("param"))

  # perhaps some bias, but the true values are within the (marginal) regions of high posterior density

  naive_scatter <- res$draws() |>
    tidybayes::tidy_draws() |>
    ggplot(aes(x  = meanlog, y = sdlog)) +
    geom_point() +
    geom_density_2d() +
    geom_point(data = true_params, colour = "red", pch = 15)

  results <- list()
  results$plot <- naive_scatter
  results$posterior <- res
  results$input_data <- list(
    n = nrow(stan_data),
    y = stan_data$onset_to_hosp
  )

  results
}

overall_error <- function(observed, expected) {
  top <- (observed - expected)^2

  sum(top / expected)
}
