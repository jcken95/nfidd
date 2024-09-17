// 03-gamma.stan
data {
  int<lower=0> n; // number of data points
  array[n] real y; // data
}
parameters {
  real<lower=0> alpha;
  real<lower=0> beta;
}

model {
  alpha ~ normal(0, 1);  // prior distribution
  beta ~ normal(0, 1); // prior distribution

  y ~ gamma(alpha, beta);
}
