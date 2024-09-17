# Delay distributions
# https://nfidd.github.io/nfidd/sessions/delay-distributions.html

# setup ----

source(here::here("R", "00-setup.R"))


# load data ----
data(infection_times)
head(infection_times)

### visualise the infection curve

# infection_times are (decimal) values indicating the time at which an individual was
# infected, relative to the first time of infection

ggplot(infection_times, aes(x = infection_time)) +
  geom_histogram(binwidth = 1) +
  scale_x_continuous(n.breaks = 10) +
  labs(x = "Infection time (in days)", y = "Number of infections",
       title = "Infections during an outbreak")

# in this case, approx symmetric
dplyr::summarise(
  infection_times,
  mean = mean(infection_time),
  sd = sd(infection_time),
  median = median(infection_time),
  iqr = IQR(infection_time),
  infection_time |>
    stats::quantile(probs = c(5, 25, 75, 95)  / 100) |>
    as.list() |>
    tibble::as_tibble() |>
    dplyr::rename_with(
      \(str) {
        str_name <- str |>
          stringr::str_remove("%") |>
          stringr::str_pad(2, "left", "0")
        glue::glue(".q{str_name}")
      }
    )

)

# simulation of delays after an infection ----

## goal: simulation of hospital arrivals after an outbreak

## assumptions:
##  * every infection -> symtoms
##  * incubation period  ~ gamma(shape = 5, rate = 1)
##  * p(hospitalised | symptoms) = 0.3
##  * time from symptom onset -> admission ~ lognormal(mu = 1.75, sigma = 0.5)

# Simulate times for symptom onset and hospitalisation for each infection

# remember: onset ~ gamma
# time to hosp | symptoms present ~ I(hopsitalised) * lognormal + (1 - I(hosptialised)) * NULL

source(
  here::here("nfidd", "snippets", "onset-hosp.r"),
  echo = TRUE,
  max.deparse.length = 1000
)

# convert our data frame to long format
dfl <- df |>
  pivot_longer(
    cols = c(infection_time, onset_time, hosp_time),
    names_to = "type", values_to = "time"
  )
# plot
ggplot(dfl, aes(x = time)) +
  geom_histogram(position = "dodge", binwidth = 1) +
  facet_wrap(~ type, ncol = 1) +
  xlab("Time (in days)") +
  ylab("Count")

# as expected, around 30% of patients are hospitalised

1 - mean(is.na(df$hosp_time))

# estimation of delay distribution from data ----

## lognormal model ----

# aim: recover model parameters (time to hospitalisation)

mod <- cmdstan_model(here("nfidd", "stan", "lognormal.stan"))
mod
# model outline
# y = time to hospitalisation after symptom first present
# likelihood:
# y ~ lognormal(mu, sig)
# prior:
# mu ~ normal(0, 10)
# sig ~ half-normal(0, 10)

df_onset_to_hosp <- df |>
  mutate(onset_to_hosp = hosp_time - onset_time) |>
  # exclude infections that didn't result in hospitalisation
  drop_na(onset_to_hosp)
## Use the data to sample from the model posterior
res <- mod$sample(
  data = list(
    n = nrow(df_onset_to_hosp),
    y = df_onset_to_hosp$onset_to_hosp
  )
)

# true params:
# mu  = 1.75 (meanlog)
# sig = 0.5 (sdlog)

res$summary()

# posterior means of these params close to theoretical values
# E(mu, sig | data) = (1.73, 0.494)

## get shape and rate samples from the posterior
res_df <- res |>
  as_draws_df() |>
  filter(.draw %in% sample(.draw, 100)) # sample 100 draws

## find the value (x) that includes 99% of the cumulative density
max_x <- max(qlnorm(0.99, meanlog = res_df$meanlog, sdlog = res_df$sdlog))

## calculate density on grid of x values
x <- seq(0, max_x, length.out = 100)
res_df <- res_df |>
  crossing(x = x) |> ## add grid to data frame
  mutate(density = dlnorm(x, meanlog, sdlog))

## plot
lognormal_density <- ggplot(res_df, aes(x = x, y = density, group = .draw)) +
  geom_line(alpha = 0.3) +
  ggtitle("lognormal")

# Exercise:
#  In this session we were in the enviable situation of knowing which
# distribution was used to generate the data. With real data, of course,
# we don’t have this information available. Try using a different distribution
# for inference (e.g. normal, or gamma). Do you get a good fit?

##  gamma model ----

gamma_mod <- cmdstan_model(here("stan", "03-gamma.stan"))
gamma_mod


df_onset_to_hosp_gamma <- df |>
  mutate(onset_to_hosp = hosp_time - onset_time) |>
  # exclude infections that didn't result in hospitalisation
  drop_na(onset_to_hosp)
## Use the data to sample from the model posterior
res_gamma <- gamma_mod$sample(
  data = list(
    n = nrow(df_onset_to_hosp),
    y = df_onset_to_hosp$onset_to_hosp
  )
)

res_gamma$summary()

# posterior means of these params close to theoretical values
# E(mu, sig | data) = (1.73, 0.494)

## get shape and rate samples from the posterior
res_df_gamma <- res_gamma |>
  as_draws_df() |>
  filter(.draw %in% sample(.draw, 100)) # sample 100 draws

## find the value (x) that includes 99% of the cumulative density
#max_x_gamma <- max(qgamma(0.99999, shape = res_df_gamma$alpha, scale = res_df_gamma$beta))

## calculate density on grid of x values
x <- seq(0, max_x, length.out = 100)
res_df_gamma <- res_df_gamma |>
  crossing(x = x) |> ## add grid to data frame
  mutate(density = dgamma(x, alpha, beta))

## plot
gamma_density <- ggplot(res_df_gamma, aes(x = x, y = density, group = .draw)) +
  geom_line(alpha = 0.3) +
  ggtitle("gamma")
# not too different from the log normal
# note if priors are (for gamma model) N(0, 10), as with lognormal
# results are very different


## normal model ----

normal_mod <- cmdstan_model(here("stan", "03-normal.stan"))
normal_mod


df_onset_to_hosp_normal <- df |>
  mutate(onset_to_hosp = hosp_time - onset_time) |>
  # exclude infections that didn't result in hospitalisation
  drop_na(onset_to_hosp)
## Use the data to sample from the model posterior
res_normal <- normal_mod$sample(
  data = list(
    n = nrow(df_onset_to_hosp),
    y = df_onset_to_hosp$onset_to_hosp
  )
)

res_normal$summary()

# posterior means of these params close to theoretical values
# E(mu, sig | data) = (1.73, 0.494)

## get shape and rate samples from the posterior
res_df_normal <- res_normal |>
  as_draws_df() |>
  filter(.draw %in% sample(.draw, 100)) # sample 100 draws

## find the value (x) that includes 99% of the cumulative density
#max_x_normal <- max(qnormal(0.99999, shape = res_df_normal$alpha, scale = res_df_normal$beta))

## calculate density on grid of x values
x <- seq(0, max_x, length.out = 100)
res_df_normal <- res_df_normal |>
  crossing(x = x) |> ## add grid to data frame
  mutate(density = dnorm(x, mu, sigma))

## plot
normal_density <- ggplot(res_df_normal, aes(x = x, y = density, group = .draw)) +
  geom_line(alpha = 0.3) +
  ggtitle("normal")
# not too different from the log normal
# note if priors are (for normal model) N(0, 10), as with lognormal
# results are very different

wrap_plots(lognormal_density, gamma_density, normal_density)

# gamma and lognormal are reasonably simlar - hard to tell apart by eye
# normal isn't as good
#  * non negligible density below zero
#  * normal symmetric whereas truth (and gamma) are skewed
