# Forecasting concepts
# https://nfidd.github.io/nfidd/sessions/forecasting-concepts.html

source(here::here("R", "00-setup.R"))

library("dplyr")
library("tidyr")
library("ggplot2")
library("here")
library("cmdstanr")
library("tidybayes")
library("scoringutils")

set.seed(123)
options(cmdstanr_print_line_numbers = TRUE)

# Good forecast

# .1 Calibration - forecasted probabilities match observed frequencies
# .2 unbiased - should not consistently over or under predict
# .3 accuracy - forecasted value ~= observed value
# .4 sharpness - forecast intervals should be as narrow as possible
#                excessively wide intervals may not be useful

mod <- cmdstan_model(here("nfidd", "stan", "estimate-inf-and-r-rw-forecast.stan"))
mod

# model modifications:

## `h` allows us to predict `h` time periods into the future
## llh has `n` obs => only use first `n` elts of `onsets`
## generated quantities block allows us to forecast

# potential limitations

## including `h` in model and param blocks increases copmutational load
## vs a separate forecasting model

## bodge because my structure means i can't source things with subprojects

source(here::here("nfidd", "snippets", "load-ts.r"))
source(here::here("nfidd", "functions", "censored-delay-pmf.r"))
source(here::here("nfidd", "functions", "convolve-with-delay.r"))

gen_time_pmf <- censored_delay_pmf(rgamma, max = 14, shape = 4, rate = 1)
gen_time_pmf <- gen_time_pmf[-1] ## remove first element
gen_time_pmf <- gen_time_pmf / sum(gen_time_pmf) ## renormalise

ip_pmf <- censored_delay_pmf(rgamma, max = 14, shape = 5, rate = 1)
onsets <- convolve_with_delay(inf_ts$infections, ip_pmf)
onsets <- rpois(n = length(onsets), lambda = onsets)
onset_df <- tibble(day = seq_along(onsets), onsets = onsets) |>
  left_join(
    inf_ts |> select(day = infection_day, infections),
    by = "day"
  ) |>
  replace_na(list(infections = 0))
## end bodge

cutoff <- 71
filtered_onset_df <- onset_df |>
  filter(day <= cutoff)
tail(onset_df)


## forecast into future

horizon <- 28

data <- list(
  n = nrow(filtered_onset_df),
  I0 = 1,
  obs = filtered_onset_df$onsets,
  gen_time_max = length(gen_time_pmf),
  gen_time_pmf = gen_time_pmf,
  ip_max = length(ip_pmf) - 1,
  ip_pmf = ip_pmf,
  h = horizon # Here we set the number of days to forecast into the future
)
rw_forecast <- mod$sample(
  data = data, parallel_chains = 4, adapt_delta = 0.95,
  init = \() list(init_R = 0, rw_sd = 0.01)
)


rw_forecast

## plot forecast

forecast <- rw_forecast |>
  gather_draws(forecast[day]) |>
  ungroup() |>
  mutate(day = day + cutoff)

target_onsets <- onset_df |>
  filter(day > cutoff) |>
  filter(day <= cutoff + horizon)

forecast |>
  filter(.draw %in% sample(.draw, 100)) |>
  ggplot(aes(x = day)) +
  geom_line(alpha = 0.1, aes(y = .value, group = .draw)) +
  geom_point(data = target_onsets, aes(x = day, y = onsets), color = "black") +
  labs(title = "Symptom onsets", subtitle = "Forecast (trajectories) and observed (points)")

# quite high uncertainty
# often outbreaks are exponential => plot on log scale

forecast |>
  filter(.draw %in% sample(.draw, 100)) |>
  ggplot(aes(x = day)) +
  geom_line(alpha = 0.1, aes(y = .value, group = .draw)) +
  geom_point(data = target_onsets, aes(x = day, y = onsets), color = "black") +
  scale_y_log10() +
  labs(title = "Symptom onsets, log scale", subtitle = "Forecast and observed")


## forecast ok in short term - after ~ 80 days starts to strugle
## .e. systematic underprediction

## we have forecastied based on reproduction number - how is this doing?

long_onset_df <- onset_df |>
  filter(day <= cutoff + horizon)

long_data <- list(
  n =nrow(long_onset_df),
  I0 = 1,
  obs = long_onset_df$onsets,
  gen_time_max = length(gen_time_pmf),
  gen_time_pmf = gen_time_pmf,
  ip_max = length(ip_pmf) - 1,
  ip_pmf = ip_pmf,
  h = 0
)


rw_long <- mod$sample(
  data = long_data, parallel_chains = 4, adapt_delta = 0.95,
  init = \() list(init_R = 0, rw_sd = 0.01)
)

## grab forecast & repro numbers

forecast_r <- rw_forecast |>
  gather_draws(R[day]) |>
  ungroup() |>
  mutate(type = "forecast")

long_r <- rw_long |>
  gather_draws(R[day]) |>
  ungroup() |>
  mutate(type = "estimate")

# plot

forecast_r |>
  bind_rows(long_r) |>
  filter(.draw %in% sample(.draw, 100)) |>
  ggplot(aes(x = day)) +
  geom_vline(xintercept = cutoff, linetype = "dashed") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_line(aes(y = .value, group = interaction(.draw, type), color = type), alpha = 0.1)+
  labs(title = "Estimated R",
       subtitle = "Estimated over whole time series (red), and forecast (blue)") +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))


## before vertical line (the point at which we start forecasting)
## estiamtes / forecast are similar

## estimate has a downward trend whereas forecast has approx constant mean

## makes sense because randow walks (by defn) have constant mean
# but increasing variance


# maybe we need rw + drift??

# Forecast evaluation

data(rw_forecasts)
rw_forecasts |>
  ungroup()
# visualist

rw_forecasts |>
  filter(.draw %in% sample(.draw, 100)) |>
  ggplot(aes(x = day)) +
  geom_line(aes(y = .value, group = interaction(.draw, target_day), col = target_day), alpha = 0.1) +
  geom_point(data = onset_df |>
               filter(day >= 21),
             aes(x = day, y = onsets), color = "black") +
  scale_color_binned(type = "viridis") +
  labs(title = "Weekly forecasts of symptom onsets over an outbreak",
       col = "Forecast start day")
### and on log scale

rw_forecasts |>
  filter(.draw %in% sample(.draw, 100)) |>
  ggplot(aes(x = day)) +
  geom_line(aes(y = .value, group = interaction(.draw, target_day), col = target_day), alpha = 0.1) +
  geom_point(data = onset_df, aes(x = day, y = onsets), color = "black") +
  scale_y_log10() +
  scale_color_binned(type = "viridis") +
  labs(title = "Weekly symptom onset forecasts: log scale",
       col = "Forecast start day")

## forecasts + uncertainty look good - perhaps a little too much uncertainty?
## better in the short term than they are in the long term - not surprising
## some suggestion of overprediction on the raw scale - less clear on log

# proper scoring rules for forecasts
sc_forecasts <- rw_forecasts |>
  left_join(onset_df, by = "day") |>
  filter(!is.na(.value)) |>
  as_forecast_sample(
    forecast_unit = c("target_day", "horizon", "model"),
   # forecast_type = "sample",
    observed = "onsets",
    predicted = ".value",
    model = "model",
    sample_id = ".draw"
  )
sc_forecasts


sc_scores <- sc_forecasts |>
  score()

sc_scores


summarise_scores(sc_scores, by = "model")


# CPRS

## enrealiseation of MAE to a distribution

sc_scores |>
  summarise_scores(by = "horizon") |>
  ggplot(aes(x = horizon, y = crps)) +
  geom_point() +
  labs(title = "CRPS by daily forecast horizon",
       subtitle = "Summarised across all forecasts")


sc_scores |>
  summarise_scores(by = "target_day") |>
  ggplot(aes(x = target_day, y = crps)) +
  geom_point() +
  labs(title = "CRPS by forecast start date",
       subtitle = "Summarised across all forecasts", x = "forecast date")

# probability integral transofrm

# distribution should be U(0, 1) if asumptions/fit are okay

sc_forecasts |>
  get_pit(by = "model") |>
  plot_pit() +
  labs(title = "PIT histogram")

sc_forecasts |>
  mutate(group_horizon = case_when(
    horizon <= 3 ~ "1-3",
    horizon <= 7 ~ "4-7",
    horizon <= 14 ~ "8-14"
  )) |>
  get_pit(by = "group_horizon") |>
  plot_pit() +
  facet_wrap(~group_horizon) +
  labs(title = "PIT by forecast horizon (days)")

sc_forecasts |>
  get_pit(by = "target_day") |>
  plot_pit() +
  facet_wrap(~target_day) +
  labs(title = "PIT by forecast date")
# PIT indicate that the model is overpredicitn
# bias increased at longer oforecast lags

## score on log cale

log_sc_forecasts <- sc_forecasts |>
  transform_forecasts(
    fun = log_shift,
    offset = 1,
    append = FALSE
  )

log_scores <- log_sc_forecasts |>
  score()

log_scores |>
  summarise_scores(by = "model")

log_scores |>
  summarise_scores(by = "horizon") |>
  ggplot(aes(x = horizon, y = crps)) +
  geom_point() +
  labs(title = "CRPS by daily forecast horizon, scored on the log scale")
log_scores |>
  summarise_scores(by = "target_day") |>
  ggplot(aes(x = target_day, y = crps)) +
  geom_point() +
  labs(title = "CRPS by forecast date, scored on the log scale")
# pit for log scale

log_sc_forecasts |>
  get_pit(by = "model") |>
  plot_pit() +
  labs(title = "PIT histogram, scored on the log scale")


log_sc_forecasts |>
  mutate(group_horizon = case_when(
    horizon <= 3 ~ "1-3",
    horizon <= 7 ~ "4-7",
    horizon <= 14 ~ "8-14"
  )) |>
  get_pit(by = "group_horizon") |>
  plot_pit() +
  facet_wrap(~group_horizon) +
  labs(title = "PIT by forecast horizon, scored on the log scale")


log_sc_forecasts |>
  get_pit(by = "target_day") |>
  plot_pit() +
  facet_wrap(~target_day) +
  labs(title = "PIT by forecast date, scored on the log scale")
