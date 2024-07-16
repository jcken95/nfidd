# Bias in delay distributions
# https://nfidd.github.io/nfidd/sessions/biases-in-delay-distributions.html

# setup ----

source(here::here("R", "00-setup.R"))

# load data ----

# simulate infections

source(here("nfidd", "snippets", "onset-hosp.r"))

head(df)

# dates not days ----

# commonly, we don't have the exact infection time. We normally have dates
# this sensors the observation

# we make the data in df more realistic by rounding the date down to integer
# in practice, we would work out the number of days between two key dates

df_dates <- mutate(df, across(everything(), floor))
head(df_dates)

# our data is now double interval sensored
#  * double: delays represent the time between two events that are both censored
#  * interval: we only know the time intercal at which the event happened (between 00:00 and 23:59)


# estimating delay distributions - accounting for sensoring

# naive approach: ignore the fact that data are censored
# just compute difference between censored times

df_dates <- df_dates |>
  mutate(
    incubation_period = onset_time - infection_time,
    onset_to_hosp = hosp_time - onset_time
  )

# Excersise 1 ----
# Fit the lognormal distribution - are the parameters recovered?


mod <- cmdstan_model(here("nfidd", "stan", "lognormal.stan"))
mod
# model outline
# y = time to hospitalisation after symptom first present
# likelihood:
# y ~ lognormal(mu, sig)
# prior:
# mu ~ normal(0, 10)
# sig ~ half-normal(0, 10)

stan_data <- drop_na(df_dates, onset_to_hosp)
## Use the data to sample from the model posterior
res <- mod$sample(
  data = list(
    n = nrow(stan_data),
    y = stan_data$onset_to_hosp
  )
)

# true params:
# mu  = 1.75 (meanlog)
# sig = 0.5 (sdlog)
true_params <- tibble(param = c("meanlog", "sdlog"), value = c(1.75, 0.5))
true_params_wide <- pivot_wider(true_params,
                                names_from = "param")
res$summary()

posterior_samples <- res$draws() |>
  tidybayes::tidy_draws() |>
  select(-lp__) |>
  pivot_longer(
    cols = c(meanlog, sdlog),
    values_to = "value",
    names_to = "param"
  ) |>
  left_join(true_params,
            suffix = c("", "_true"),
            by = c("param"))

posterior_samples |>
  ggplot() +
  geom_histogram(aes(x = value)) +
  geom_vline(aes(xintercept = value_true)) +
  facet_wrap(~param, scales = "free")

# perhaps some bias, but the true values are within the (marginal) regions of high posterior density

naive_scatter <- res$draws() |>
  tidybayes::tidy_draws() |>
  ggplot(aes(x  = meanlog, y = sdlog)) +
  geom_point() +
  geom_density_2d() +
  geom_point(data = true_params_wide, colour = "red", pch = 15)

# estimation with double interval censoring ----

cmod <- cmdstan_model(here("nfidd", "stan", "censored-delay-model.stan"))
cmod
# same log normal model as before, but now the onset times are
# true_time = integer time of hospital admission + hospital time of day - onset time of day

cres <- cmod$sample(
  data = list(
    n = nrow(na.omit(df_dates)),
    onset_to_hosp = na.omit(df_dates)$onset_to_hosp
  )
)

cres$summary()

uniform_scatter <- cres$draws() |>
  tidybayes::tidy_draws() |>
  ggplot(aes(x  = meanlog, y = sdlog)) +
  geom_point() +
  geom_density_2d() +
  geom_point(data = true_params_wide, colour = "red", pch = 15)

patchwork::wrap_plots(
  naive_scatter + ggtitle("naive"),
  uniform_scatter + ggtitle("uniform")
)


# Exercise 2 ----
# Try re-simulating the delays using different parameters of the delay distribution.
# Can you establish under which conditions the bias in estimation gets worse?

params <- expand_grid(
  meanlog = seq(2, 4,  by = 0.5),
  sdlog = seq(0.3, 0.7, by = 0.1)
)

simulation <- params |>
  mutate(
    aux_posterior = map2(
      params$meanlog,
      params$sdlog,

      \(m, s) {
        withr::with_seed(4321, fit_model(mod, true_params = tibble(meanlog = m, sdlog = s)))
      }
    ),
    plot = map(aux_posterior, \(.x) .x$plot),
    posterior = map(aux_posterior, \(.x) .x$posterior),
    input_data = map(aux_posterior, \(.x) .x$input_data)
  ) |>
  select(-aux_posterior) |>
  mutate(
    e_meanlog = map(posterior, \(post) mean_from_summary(post, "meanlog")),
    e_sdlog  = map(posterior, \(post)  mean_from_summary(post, "sdlog"))
  ) |>
  unnest(c(e_meanlog, e_sdlog))

simulation <- simulation |>
  mutate(
    mean = exp(meanlog + 0.5 * sdlog^2),
    variance = (exp(sdlog^2) - 1) * exp(sdlog^2 + 2 * meanlog)
  )

simulation <- simulation |>
  rowwise() |>
  mutate(overall_bias = overall_error(c(e_meanlog, e_sdlog), c(meanlog, sdlog))) |>
  ungroup()

simulation |>
  ggplot(aes(x = meanlog, y = sdlog, fill = overall_bias)) +
  geom_tile()

simulation |>
  ggplot(aes(x = meanlog, y = overall_bias, colour = sdlog)) +
  geom_point()

# reduced relative error as meanlog increases

simulation |>
  ggplot(aes(x = sdlog, y = overall_bias, colour = meanlog)) +
  geom_point()

# reduced relative error as sd log increases (maybe?)

simulation |>
  ggplot(aes(x = mean, y = overall_bias, colour = variance)) +
  geom_point()

# larger mean -> reduced relative error

simulation |>
  ggplot(aes(x = variance, y = overall_bias, colour = mean)) +
  geom_point()

# larger var -> reduced relative error (just about ...)
