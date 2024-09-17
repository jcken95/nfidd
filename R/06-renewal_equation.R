# R_t estimateion & the renewal equation
# https://nfidd.github.io/nfidd/sessions/R-estimation-and-the-renewal-equation.html

# setup ----

source(here::here("R", "00-setup.R"))

# The renewal equation as a process model for infectious diseases ----

# Simulating an epidemic using the renewal equation  ----

source(here("nfidd", "functions", "renewal.r"))
renewal

# Estimating R_t from a time series of infections ----

source(here("nfidd", "snippets", "load-ts.r"))
head(inf_ts)
source(here("nfidd", "functions", "censored-delay-pmf.r"))
gen_time_pmf <- censored_delay_pmf(rgamma, max = 14, shape = 4, rate = 1)

gen_time_pmf <- gen_time_pmf[-1] ## remove first element
gen_time_pmf <- gen_time_pmf / sum(gen_time_pmf) ## renormalise

r_mod <- cmdstan_model(here("nfidd", "stan", "estimate-r.stan"))
r_mod

data <- list(
  n = nrow(inf_ts) - 1,
  obs = inf_ts$infections[-1],
  I0 = inf_ts$infections[1],
  gen_time_max = length(gen_time_pmf),
  gen_time_pmf = gen_time_pmf
)
r_fit <- r_mod$sample(data = data, parallel_chains = 4)

r_fit

# Extract posterior draws
r_posterior <- r_fit |>
  gather_draws(R[infection_day]) |>
  ungroup() |>
  mutate(infection_day = infection_day - 1) |>
  filter(.draw %in% sample(.draw, 100))

ggplot(
  data = r_posterior,
  aes(x = infection_day, y = .value, group = .draw))  +
  geom_line(alpha =  0.1) +
  labs(title = "Estimated Rt",
       subtitle = "Model 1: renewal equation from infections")
