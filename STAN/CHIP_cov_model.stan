data {
    int N;
    int Y[N];
}
parameters {
    vector[N] r;
    real<lower=0> s_r;
}
model {
    target += normal_lpdf(r[3:N] | 2*r[2:(N-1)] - r[1:(N-2)], s_r);
    Y ~ poisson_log(r);
}
generated quantities {
  vector[N] Y_mean;
  Y_mean = exp(r);
}
