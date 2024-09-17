# Using delay distributions to model the data generating process
# https://nfidd.github.io/nfidd/sessions/using-delay-distributions-to-model-the-data-generating-process-of-an-epidemic.html

# setup ----

source(here::here("R", "00-setup.R"))

# Delay distributions & convolutions ----

source(
  here::here("nfidd", "snippets", "load-ts.r"),
  echo = TRUE, max.deparse.length = 1000
)

head(inf_ts)

## discretising a delay distn ----

# load in censored_delay_pmf()
source(here::here("nfidd", "functions", "censored-delay-pmf.r"))
censored_delay_pmf

# generate the density of delay (in days) from infection to symptom onset
# assumptions:
##  1. maximum delay is 14 days
##  2. the exact delay distribution is gamma(5, 1) - trunctated above at 14
gamma_pmf <- censored_delay_pmf(rgamma, max = 14, shape = 5, rate = 1)
gamma_pmf
plot(gamma_pmf)

# small probab of being > 14 days with a longer delay

gamma_pmf_28 <- censored_delay_pmf(rgamma, max = 28, shape = 5, rate = 1)
gamma_pmf_28
plot(gamma_pmf_28)


lognormal_pdf <- censored_delay_pmf(rlnorm, max = 14, meanlog = 0.5, sdlog = 1)
lognormal_pdf
plot(lognormal_pdf)

## Applying a convolution

source(here::here("nfidd", "functions", "convolve-with-delay.r"))
convolve_with_delay

onsets <- convolve_with_delay(inf_ts$infections, gamma_pmf)

combined <- inf_ts |>
  rename(time = infection_day) |>
  mutate(onsets = onsets)
ggplot(combined, aes(x = time, y = onsets)) +
  geom_bar(stat = "identity")

# observation uncertainty ----

# convolution is deterministic
# add some randomness - we have counts so a poisson is a sensible choice

combined <- combined |>
  mutate(observed = rpois(n(), onsets))
ggplot(combined, aes(x = time, y = observed)) +
  geom_bar(stat = "identity")

# Estimating a time series of infections

# don't usually have data on number of infections
# use symptom onsets (what we normall have) to estimate number of infections over time

mod <- cmdstan_model(here("nfidd", "stan", "estimate-infections.stan"))
mod
data <- list(
  n = nrow(combined),
  obs = combined$observed,
  ip_max = length(gamma_pmf) - 1,
  ip_pmf = gamma_pmf
)

inf_fit <- mod$sample(data = data, parallel_chains = 4)

inf_fit

# Extract posterior draws
inf_posterior <- inf_fit |>
  gather_draws(infections[infection_day]) |>
  ungroup() |>
  mutate(infection_day = infection_day - 1) |>
  filter(.draw %in% sample(.draw, 100))

ggplot(mapping = aes(x = infection_day)) +
  geom_line(
    data = inf_posterior, mapping = aes(y = .value, group = .draw), alpha = 0.1
  ) +
  geom_line(data = inf_ts, mapping = aes(y = infections), colour = "red") +
  labs(title = "Infections per day",
       subtitle = "True data (red), and sample of estimated trajectories (black)")
