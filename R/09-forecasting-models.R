# forecasting models
# https://nfidd.github.io/nfidd/sessions/forecasting-models.html

library("nfidd")
library("dplyr")
library("tidyr")
library("ggplot2")
library("here")
library("cmdstanr")
library("tidybayes")
library("scoringutils")
set.seed(123)
options(cmdstanr_print_line_numbers = TRUE)

# Forecastin as a spectrum

# idea: use a mixture / combintino of mechanical and statistical models
# to make forecasts

# hopin for best of both worlds

# identifying best model can be hard


# add mech structure to renewal model

# log r_t = log s_{t-1} + log r_0 - log N

# this term is proportinoal to the number of susceptibles in the population
# basic idea behind the SIR model

# simulation

n <- 100
N <- 1000
R0 <- 1.5
S <- rep(NA, n)
S[1] <- N
Rt <- rep(NA, n) ## reproduction number
Rt[1] <- R0
I <- rep(NA, n)
I[1] <- 1
for (i in 2:n) {
  Rt[i] <- (S[i-1]) / N * R0
  I[i] <- I[i-1] * Rt[i]
  S[i] <- S[i-1] - I[i]
}

data <- tibble(t = 1:n, Rt = Rt)

ggplot(data, aes(x = t, y = Rt)) +
  geom_line() +
  labs(title = "Simulated data from an SIR model",
       x = "Time",
       y = "Rt")
#fit in stan
mech_mod <- cmdstan_model(here("nfidd", "stan", "mechanistic-r.stan"))
mech_mod

# forecasting with mech/stat models

source(here("nfidd", "snippets", "simulate-onsets.r"))
onset_df

# we'll make a forecast on day non day 41, pretending we haven't seen the later data
cutoff <- 41

filtered_onset_df <- onset_df |>
  filter(day <= cutoff)

# fit mostly mech model
horizon <- 28

data <- list(
  n =nrow(filtered_onset_df),
  I0 = 1,
  obs = filtered_onset_df$onsets,
  gen_time_max = length(gen_time_pmf),
  gen_time_pmf = gen_time_pmf,
  ip_max = length(ip_pmf) - 1,
  ip_pmf = ip_pmf,
  h = horizon, # Here we set the number of days to forecast into the future
  N_prior = c(10000, 2000) # the prior for the population size
)
mech_forecast_fit <- mech_mod$sample(
  data = data, parallel_chains = 4
)
mech_forecast_fit
# more stat model
data <- list(
  n =nrow(filtered_onset_df),
  I0 = 1,
  obs = filtered_onset_df$onsets,
  gen_time_max = length(gen_time_pmf),
  gen_time_pmf = gen_time_pmf,
  ip_max = length(ip_pmf) - 1,
  ip_pmf = ip_pmf,
  h = horizon # Here we set the number of days to forecast into the future
)

stat_mod <- cmdstan_model(here("nfidd", "stan/statistical-r.stan"))
stat_mod

stat_forecast_fit <- stat_mod$sample(
  data = data, parallel_chains = 4,
  init = \() list(init_R = 0, rw_sd = 0.01) # again set the initial values to make fitting more numerically stable
)

# extract & plt forecasts
mech_forecast <- mech_forecast_fit |>
  gather_draws(forecast[day]) |>
  ungroup() |>
  mutate(day = day + cutoff)

stat_forecast <- stat_forecast_fit |>
  gather_draws(forecast[day]) |>
  ungroup() |>
  mutate(day = day + cutoff)

forecast <- bind_rows(
  mutate(mech_forecast, model = "more mechanistic"),
  mutate(stat_forecast, model = "more statistical")
) |>
  ungroup()

target_onsets <- onset_df |>
  filter(day > cutoff) |>
  filter(day <= cutoff + horizon)
forecast |>
  filter(.draw %in% sample(.draw, 100)) |>
  ggplot(aes(x = day)) +
  geom_line(alpha = 0.1, aes(y = .value, group = interaction(.draw, model), colour = model)) +
  geom_point(data = target_onsets, aes(x = day, y = onsets), color = "black") +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  lims(y = c(0, 500))
## on log scale
forecast |>
  filter(.draw %in% sample(.draw, 100)) |>
  ggplot(aes(x = day)) +
  geom_line(alpha = 0.1, aes(y = .value, group = interaction(.draw, model), colour = model)) +
  geom_point(data = target_onsets, aes(x = day, y = onsets), color = "black") +
  scale_y_log10() +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))

## mechanical model seems better, but did say we were cheating

data(rw_forecasts, stat_forecasts, mech_forecasts)
forecasts <- bind_rows(
  mutate(rw_forecasts, model = "Random walk"),
  mutate(stat_forecasts, model = "More statistical"),
  mutate(mech_forecasts, model = "More mechanistic")
) |>
  ungroup()

head(onset_df)


forecasts |>
  filter(.draw %in% sample(.draw, 100)) |>
  ggplot(aes(x = day)) +
  geom_line(aes(y = .value, group = interaction(.draw, target_day), col = target_day), alpha = 0.1) +
  geom_point(data = onset_df |>
               filter(day >= 21),
             aes(x = day, y = onsets), color = "black") +
  scale_color_binned(type = "viridis") +
  facet_wrap(~model) +
  lims(y = c(0, 500))


forecasts |>
  filter(.draw %in% sample(.draw, 100)) |>
  ggplot(aes(x = day)) +
  geom_line(aes(y = .value, group = interaction(.draw, target_day), col = target_day), alpha = 0.1) +
  geom_point(data = onset_df, aes(x = day, y = onsets), color = "black") +
  scale_y_log10(limits = c(NA, 500)) +
  scale_color_binned(type = "viridis") +
  facet_wrap(~model)


# mechanistic model captures downturn really well
# uncertainty: rw > stat > mech
# which is best would probably depend on our prior understanding of the epidemic?

sc_forecasts <- forecasts |>
  left_join(onset_df, by = "day") |>
  filter(!is.na(.value)) |>
  as_forecast(
    forecast_unit = c("target_day", "horizon", "model"),
    forecast_type = "sample",
    observed = "onsets",
    predicted = ".value",
    model = "model",
    sample_id = ".draw"
  )
sc_forecasts


sc_scores <- sc_forecasts |>
  score()

sc_scores

sc_scores |>
  summarise_scores(by = "model")

sc_scores |>
  summarise_scores(by = c("model", "horizon")) |>
  ggplot(aes(x = horizon, y = crps, col = model)) +
  geom_point()

sc_scores |>
  summarise_scores(by = c("target_day", "model")) |>
  ggplot(aes(x = target_day, y = crps, col = model)) +
  geom_point()

sc_forecasts |>
  get_pit(by = "model") |>
  plot_pit() +
  facet_wrap(~model)
