data {
    int Np;
    int Nt;
    real T[Nt];
    real Y[Np, Nt];
    int new_T;
    real new_Time[new_T];
}
parameters {
    real m_a;
    real m_b;
    real log_a[Np];
    real log_b[Np];
    real<lower=0> s_a;
    real<lower=0> s_b;
    real <lower=0> s_Y;
}
transformed parameters {
    real a[Np];
    real b[Np];
    for (i in 1:Np) {
        a[i] = exp(log_a[i]);
        b[i] = exp(log_b[i]);
    }

}
model {
    for (i in 1:Np) {
        log_a[i] ~ normal(m_a, s_a);
        log_b[i] ~ normal(m_b, s_b);
    }

    for (i in 1:Np)
        for (j in 1:Nt)
            Y[i, j] ~ normal(a[i]*(1-exp(-b[i]*T[j])), s_Y);
}
generated quantities {
    real new_Y[Np, new_T];
    for (i in 1:Np)
        for (j in 1:new_T)
            new_Y[i, j] = normal_rng(a[i]*(1-exp(-b[i]*new_Time[j])), s_Y);
}
