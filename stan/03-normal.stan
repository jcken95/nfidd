// 03-normal.stan
data {
  int<lower=0> n; // number of data points
  array[n] real y; // data
}
parameters {
  real mu;
  real<lower=0> sigma;
}

model {
  mu ~ normal(0, 1);  // prior distribution
  sigma ~ normal(0, 1); // prior distribution

  y ~ normal(mu, sigma);
}
